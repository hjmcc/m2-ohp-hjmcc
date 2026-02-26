"""Lightweight FITS directory scanner for at-telescope use.

Walks a data directory, reads FITS headers, classifies frames (bias/flat/science),
and builds an archive JSON compatible with the calibration pipeline.

Usage:
    python -m pipeline scan --year 2025 --telescope T120
"""

import logging
import re
from collections import defaultdict
from pathlib import Path

from astropy.io import fits

from . import config
from .utils import write_json

log = logging.getLogger("pipeline")

FITS_EXT = {".fits", ".fit"}

# Pipeline output directories to skip when scanning for raw frames
_SKIP_DIRS = {"master_bias", "master_flat", "reduced", "astrometry",
              "photometry", "psf", "qc"}

# ── Filter normalization (subset from scan_archive.py) ────────────────────

_FILTER_RULES = [
    (["R", "R_Cousins", "r_Cousins", "Rc", "RC", "Rouge", "rouge", "Red",
      "R Cousins"], "R"),
    (["V", "V_Johnson", "V_Cousins", "V_cousins", "V Cousins", "v_Gunn",
      "Vert", "vert", "Green"], "V"),
    (["B", "B_Johnson", "B_Cousins", "B_cousins", "B Cousins",
      "Bleu", "bleu", "Blue"], "B"),
    (["I", "Ic", "IC", "I_Cousins"], "I"),
    (["U", "U'_Cousins", "U''_Cousins"], "U"),
    (["g", "G", "g'", "g''", "g_Gunn"], "g"),
    (["r", "r'", "r''", "r_Gunn"], "r"),
    (["i'", "i''", "i_Gunn", "i Gunn"], "i"),
    (["z'", "z''", "z_Gunn"], "z"),
    (["u'", "u''", "u_Gunn"], "u"),
    (["SDSS g", "SDSS_g"], "SDSS_g"),
    (["SDSS r", "SDSS_r"], "SDSS_r"),
    (["SDSS i", "SDSS_i"], "SDSS_i"),
    (["SDSS z", "SDSS_z"], "SDSS_z"),
    (["SDSS u", "SDSS_u"], "SDSS_u"),
    (["H_alpha", "Halpha", "Ha", "H-alpha", "H_Alpha", "H alpha",
      "Halpha_OHP", "H alpha OHP"], "Halpha"),
    (["OIII", "O III", "[OIII]", "O_III", "OIII 5nm"], "OIII"),
    (["SII", "[SII]", "S_II", "S II"], "SII"),
    (["L", "Luminance", "luminance", "Lum"], "L"),
    (["CLEAR", "clear", "Clear", "C", "EMPTY"], "Clear"),
    (["Grism", "grism", "GRISM"], "Grism"),
]

FILTER_MAP = {}
for _variants, _canonical in _FILTER_RULES:
    for _v in _variants:
        FILTER_MAP[_v] = _canonical


def normalize_filter(raw):
    """Normalize filter name to canonical form."""
    if not raw:
        return ""
    raw = raw.strip()
    if raw in FILTER_MAP:
        return FILTER_MAP[raw]
    raw_lower = raw.lower()
    for key, val in FILTER_MAP.items():
        if raw_lower == key.lower():
            return val
    return raw


# ── Frame classification ──────────────────────────────────────────────────

def _is_flat_from_context(filename, obj_name):
    """Detect flat frames that have IMAGETYP='Light Frame' but are actually flats."""
    fn = filename.lower()
    obj = obj_name.lower().strip() if obj_name else ""
    if any(fn.startswith(p) for p in ("flatdome", "flatsky", "domeflat",
                                       "skyflat", "flat_", "flat-")):
        return True
    if obj and "flat" in re.split(r"[\s_\-]+", obj):
        return True
    return False


def classify_frame(header, filename):
    """Classify a FITS frame as bias, flat, dark, or science.

    Returns (imagetyp, object_name, filter_canonical, exptime, date_obs, binning_key).
    """
    imagetyp_raw = header.get("IMAGETYP", "")
    obj_raw = header.get("OBJECT", "")
    filter_raw = header.get("FILTER", "")
    exptime = header.get("EXPTIME", header.get("EXPOSURE", 0))
    date_obs = header.get("DATE-OBS", header.get("DATE", ""))
    xbin = header.get("XBINNING", "")
    ybin = header.get("YBINNING", "")

    try:
        exptime = float(exptime)
    except (ValueError, TypeError):
        exptime = 0.0

    # Binning key
    if xbin and ybin:
        bkey = f"{int(xbin)}x{int(ybin)}"
    else:
        bkey = "unknown"

    # Filter normalization
    filt = normalize_filter(filter_raw)

    # Object name cleanup
    obj_name = obj_raw.strip().replace("_", " ")
    obj_name = re.sub(r"\s+", " ", obj_name)

    fn = filename.lower()

    # Flat detection override (catches flats with IMAGETYP='Light Frame')
    if _is_flat_from_context(filename, obj_raw):
        return "flat", obj_name, filt, exptime, date_obs, bkey

    # Standard IMAGETYP parsing
    if imagetyp_raw:
        r = imagetyp_raw.lower().strip()
        if r in ("flat field", "flat", "flatfield"):
            return "flat", obj_name, filt, exptime, date_obs, bkey
        if r in ("bias frame", "bias"):
            return "bias", obj_name, filt, exptime, date_obs, bkey
        if r in ("dark frame", "dark"):
            return "dark", obj_name, filt, exptime, date_obs, bkey
        if r in ("light frame", "light", "object"):
            return "science", obj_name, filt, exptime, date_obs, bkey
    else:
        # No IMAGETYP — guess from filename
        if any(fn.startswith(p) for p in ("bias", "biais", "offset")):
            return "bias", obj_name, filt, exptime, date_obs, bkey
        if any(fn.startswith(p) for p in ("flat", "ff")):
            return "flat", obj_name, filt, exptime, date_obs, bkey
        if fn.startswith("dark"):
            return "dark", obj_name, filt, exptime, date_obs, bkey

    return "science", obj_name, filt, exptime, date_obs, bkey


# ── Directory scanner ─────────────────────────────────────────────────────

def scan_directory(data_dir, year, telescope):
    """Scan a telescope data directory for FITS files and classify them.

    Expected layout:
        data_dir/year/telescope/YYYYMMDD/*.fits

    Returns list of record dicts.
    """
    base = Path(data_dir) / year / telescope
    if not base.exists():
        log.error("Directory not found: %s", base)
        return []

    records = []
    fits_files = sorted(
        p for p in base.rglob("*")
        if p.suffix.lower() in FITS_EXT and p.is_file()
        and not any(part in _SKIP_DIRS for part in p.relative_to(base).parts)
    )
    log.info("Found %d FITS files in %s", len(fits_files), base)

    for i, fpath in enumerate(fits_files, 1):
        try:
            hdr = fits.getheader(str(fpath))
        except Exception as exc:
            log.warning("Cannot read header: %s (%s)", fpath.name, exc)
            continue

        imagetyp, obj_name, filt, exptime, date_obs, bkey = classify_frame(
            hdr, fpath.name
        )

        records.append({
            "filename": fpath.name,
            "path": str(fpath),
            "telescope": telescope,
            "year": year,
            "date": str(date_obs)[:10] if date_obs else "",
            "object": obj_name,
            "imagetyp": imagetyp,
            "filter": filt,
            "exptime": exptime,
            "binning": bkey,
        })

        if i % 100 == 0:
            log.info("  scanned %d/%d files...", i, len(fits_files))

    log.info("Classification: %s",
             {k: sum(1 for r in records if r["imagetyp"] == k)
              for k in ("bias", "flat", "dark", "science")})
    return records


def build_archive(records, year, telescope):
    """Build archive JSON (cal_index + file_index) from scan records.

    Returns dict matching the format consumed by calibration.py.
    """
    cal_index = {}
    file_index = defaultdict(list)

    for r in records:
        itype = r["imagetyp"]
        bkey = r["binning"]

        if itype in ("bias", "dark", "flat"):
            cal_index.setdefault(year, {})
            cal_index[year].setdefault(telescope, {})
            cal_index[year][telescope].setdefault(
                bkey, {"bias": [], "dark": [], "flat": {}}
            )
            bucket = cal_index[year][telescope][bkey]
            entry = {
                "filename": r["filename"],
                "date": r["date"],
                "exptime": r["exptime"],
                "path": r["path"],
            }
            if itype == "flat":
                filt = r["filter"] or "unknown"
                bucket["flat"].setdefault(filt, [])
                bucket["flat"][filt].append(entry)
            else:
                bucket[itype].append(entry)

        elif itype == "science" and r["object"]:
            file_index[r["object"]].append({
                "filename": r["filename"],
                "date": r["date"],
                "telescope": r["telescope"],
                "filter": r["filter"],
                "exptime": r["exptime"],
                "path": r["path"],
            })

    # Sort entries
    for y in cal_index:
        for tel in cal_index[y]:
            for bkey in cal_index[y][tel]:
                bucket = cal_index[y][tel][bkey]
                bucket["bias"].sort(key=lambda f: (f["date"], f["filename"]))
                bucket["dark"].sort(key=lambda f: (f["date"], f["filename"]))
                for filt in bucket["flat"]:
                    bucket["flat"][filt].sort(
                        key=lambda f: (f["date"], f["filename"])
                    )

    for obj in file_index:
        file_index[obj].sort(key=lambda f: (f["date"], f["filename"]))

    return {
        "cal_index": cal_index,
        "file_index": dict(file_index),
    }


def run_scan(year, telescope, data_root=None):
    """Scan a telescope directory and write archive JSON.

    Returns path to the written archive JSON.
    """
    data_dir = data_root or config.DATA_ROOT
    records = scan_directory(data_dir, year, telescope)
    if not records:
        log.error("No FITS files found — nothing to write")
        return None

    archive = build_archive(records, year, telescope)

    # Write alongside the data
    out_path = Path(data_dir) / year / telescope / "archive.json"
    write_json(archive, out_path)

    n_bias = sum(1 for r in records if r["imagetyp"] == "bias")
    n_flat = sum(1 for r in records if r["imagetyp"] == "flat")
    n_sci = sum(1 for r in records if r["imagetyp"] == "science")
    log.info("Archive written: %s", out_path)
    log.info("  %d bias, %d flat, %d science frames", n_bias, n_flat, n_sci)
    return out_path
