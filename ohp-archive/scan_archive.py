#!/usr/bin/env python3
"""scan_archive.py — OHP M2 Archive Scanner

Scans FITS files across 8 years of OHP student observing data (2018-2025),
builds a CSV database, aggregated JSON, and self-contained HTML summary.

Pure stdlib Python — reads FITS headers as raw 80-byte ASCII cards.
No astropy or other dependencies required.

Usage:
    python3 scan_archive.py
"""

import copy
import csv
import html as html_mod
import json
import os
import re
import shutil
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

ARCHIVE_ROOT = Path("/Users/hjmcc/Archive/ohp")
OUTPUT_DIR = Path(__file__).resolve().parent
CSV_PATH = OUTPUT_DIR / "ohp_archive.csv"
JSON_PATH = OUTPUT_DIR / "ohp_archive.json"
HTML_PATH = OUTPUT_DIR / "archive_summary.html"
WEB_HTML_PATH = OUTPUT_DIR / "ohp-m2-archive.html"
DOCS_DIR = OUTPUT_DIR.parent / "docs"
WEB_BASE_URL = "https://ohp.ias.universite-paris-saclay.fr/archive/"
# Local directory names that differ from the remote archive server
RUN_NAME_MAP = {
    "202503OHP1_M2": "202503OHP_M2",
}
SIMBAD_CACHE = OUTPUT_DIR / "simbad_cache.json"

FITS_EXT = {".fits", ".fit"}
MAX_HEADER_BYTES = 14400  # 5 × 2880-byte FITS blocks

TELESCOPES = {
    "T120": {
        "name": "1.2m Newton",
        "instrument": "Andor Tech CCD",
        "type": "Imaging",
        "desc": "1.2m f/6 Newton reflector with Andor 1k×1k CCD (2×2 bin)",
    },
    "T080": {
        "name": "80cm Cassegrain",
        "instrument": "SBIG STXL-6303",
        "type": "Imaging",
        "desc": "80cm f/15 Cassegrain with SBIG STXL-6303 CCD (2×2 bin)",
    },
    "IRIS": {
        "name": "50cm robotic (IRIS)",
        "instrument": "FLI ProLine",
        "type": "Imaging",
        "desc": "50cm f/8.14 robotic telescope with FLI ProLine 2k×2k CCD",
    },
    "T152": {
        "name": "1.52m Coudé",
        "instrument": "Andor DU940P_BV",
        "type": "Spectroscopy",
        "desc": "1.52m coudé spectrograph with Andor DU940P_BV echelle detector",
    },
}

CSV_COLUMNS = [
    "file_path", "filename", "year", "run_label", "telescope",
    "date_obs", "object_raw", "object_clean", "imagetyp_raw", "imagetyp",
    "filter_raw", "filter_canonical", "exptime", "ra", "dec", "airmass",
    "naxis1", "naxis2", "xbinning", "ybinning", "ccd_temp",
    "observer", "instrume",
    "simbad_name", "simbad_ra", "simbad_dec", "simbad_type",
    "error",
]

# ═══════════════════════════════════════════════════════════════════════════════
# Filter normalization
# ═══════════════════════════════════════════════════════════════════════════════

_FILTER_RULES = [
    # Cousins/Johnson broadband
    (["R", "R_Cousins", "r_Cousins", "Rc", "RC", "Rouge", "rouge", "Red",
      "R Cousins"], "R"),
    (["V", "V_Johnson", "V_Cousins", "V_cousins", "V Cousins", "v_Gunn",
      "Vert", "vert", "Green"], "V"),
    (["B", "B_Johnson", "B_Cousins", "B_cousins", "B Cousins",
      "Bleu", "bleu", "Blue"], "B"),
    (["I", "Ic", "IC", "I_Cousins"], "I"),
    (["U", "U'_Cousins", "U''_Cousins"], "U"),
    # Sloan primed (unqualified single-char handled below)
    (["g", "G", "g'", "g''", "g_Gunn"], "g"),
    (["r", "r'", "r''", "r_Gunn"], "r"),
    (["i'", "i''", "i_Gunn", "i Gunn"], "i"),
    (["z'", "z''", "z_Gunn"], "z"),
    (["u'", "u''", "u_Gunn"], "u"),
    # SDSS explicit prefix
    (["SDSS g", "SDSS_g"], "SDSS_g"),
    (["SDSS r", "SDSS_r"], "SDSS_r"),
    (["SDSS i", "SDSS_i"], "SDSS_i"),
    (["SDSS z", "SDSS_z"], "SDSS_z"),
    (["SDSS u", "SDSS_u"], "SDSS_u"),
    # Narrowband
    (["H_alpha", "Halpha", "Ha", "H-alpha", "H_Alpha", "H alpha",
      "Halpha_OHP", "H alpha OHP"], "Halpha"),
    (["OIII", "O III", "[OIII]", "O_III", "OIII 5nm"], "OIII"),
    (["SII", "[SII]", "S_II", "S II"], "SII"),
    # Luminance / clear
    (["L", "Luminance", "luminance", "Lum"], "L"),
    (["CLEAR", "clear", "Clear", "C", "EMPTY"], "Clear"),
    # Grism
    (["Grism", "grism", "GRISM"], "Grism"),
]

FILTER_MAP = {}
for _variants, _canonical in _FILTER_RULES:
    for _v in _variants:
        FILTER_MAP[_v] = _canonical

# Simbad object-type codes -> human labels
OTYPE_MAP = {
    "G": "Galaxy", "GiG": "Galaxy", "GiC": "Galaxy", "GiP": "Galaxy",
    "IG": "Interacting Galaxy", "PaG": "Galaxy Pair",
    "ClG": "Galaxy Cluster", "GrG": "Galaxy Group", "SCG": "Supercluster",
    "AGN": "AGN", "Sy1": "Seyfert 1", "Sy2": "Seyfert 2",
    "QSO": "Quasar", "BLL": "BL Lac", "LIN": "LINER",
    "*": "Star", "**": "Double Star", "V*": "Variable Star",
    "PM*": "High-PM Star", "HV*": "High-Vel Star",
    "SB*": "Spectroscopic Binary", "EB*": "Eclipsing Binary",
    "Ce*": "Cepheid", "RR*": "RR Lyrae", "LP*": "Long-Period Var",
    "WR*": "Wolf-Rayet", "Be*": "Be Star", "em*": "Emission-line Star",
    "Pe*": "Peculiar Star", "HB*": "Horiz-Branch Star",
    "RGB*": "Red Giant", "AGB*": "AGB Star", "C*": "Carbon Star",
    "WD*": "White Dwarf", "BD*": "Brown Dwarf", "TT*": "T Tauri",
    "Or*": "Orion Variable", "Ae*": "Herbig Ae/Be",
    "PN": "Planetary Nebula", "SNR": "SNR", "HII": "HII Region",
    "OpC": "Open Cluster", "GlC": "Globular Cluster", "Cl*": "Cluster",
    "As*": "Association", "Neb": "Nebula", "RNe": "Refl. Nebula",
    "GNe": "Nebula", "DNe": "Dark Nebula", "MoC": "Mol. Cloud",
    "Pl": "Planet", "Cme": "Comet",
}


def otype_label(code):
    """Convert Simbad otype short code to human-readable label."""
    if not code:
        return ""
    return OTYPE_MAP.get(code.strip(), code.strip())


# ═══════════════════════════════════════════════════════════════════════════════
# FITS header parsing
# ═══════════════════════════════════════════════════════════════════════════════

def parse_fits_header(filepath):
    """Read FITS header as raw 80-byte ASCII cards. Returns dict of keyword->value."""
    headers = {}
    try:
        with open(filepath, "rb") as f:
            raw = f.read(MAX_HEADER_BYTES)
    except Exception as e:
        return {"_error": str(e)}

    for i in range(0, len(raw), 80):
        card = raw[i : i + 80]
        if len(card) < 80:
            break
        try:
            line = card.decode("ascii", errors="replace")
        except Exception:
            continue

        # END card terminates header
        if line[:3] == "END" and line[3:].strip() == "":
            break

        key = line[:8].strip()
        if not key or key in ("HISTORY", "COMMENT"):
            continue

        # Value indicator: '= ' at columns 9-10 (standard) or '=' at column 9
        if line[8:10] == "= ":
            value_part = line[10:]
        elif len(line) > 8 and line[8] == "=":
            value_part = line[9:]
        else:
            continue

        value_part = value_part.strip()

        # String value (single-quoted)
        if value_part.startswith("'"):
            end = 1
            while end < len(value_part):
                if value_part[end] == "'":
                    if end + 1 < len(value_part) and value_part[end + 1] == "'":
                        end += 2  # escaped quote
                    else:
                        break
                else:
                    end += 1
            val = value_part[1:end].strip()
        else:
            # Numeric / logical — take everything before inline comment
            val = value_part.split("/")[0].strip()

        headers[key] = val

    return headers


# ═══════════════════════════════════════════════════════════════════════════════
# Telescope identification
# ═══════════════════════════════════════════════════════════════════════════════

def identify_telescope(filepath_str, headers):
    """Identify telescope from directory path; fallback to INSTRUME header."""
    fp = filepath_str.upper()
    if "/T120/" in fp or "/OHP_T120_" in fp:
        return "T120"
    if "/T080/" in fp or "/OHP_T080_" in fp:
        return "T080"
    if "/IRIS/" in fp:
        return "IRIS"
    if "/T152/" in fp or "/OHP_T152_" in fp:
        return "T152"

    # Fallback to header keywords
    instrume = headers.get("INSTRUME", "").lower()
    head = headers.get("HEAD", "").upper()
    if "du940" in instrume or "DU940" in head:
        return "T152"
    if "andor" in instrume:
        return "T120"
    if "sbig" in instrume or "stxl" in instrume:
        return "T080"
    if "fli" in instrume or "proline" in instrume or "flipro" in instrume:
        return "IRIS"
    return "unknown"


# ═══════════════════════════════════════════════════════════════════════════════
# Normalization helpers
# ═══════════════════════════════════════════════════════════════════════════════

def classify_t152(filename):
    """Classify T152 file type from filename prefix."""
    fn = os.path.basename(filename).lower()
    if fn.startswith(("ff", "ff_", "flat", "flatfield", "thungstene", "tungsten")):
        return "flat"
    if fn.startswith(("bias", "biais", "offset")):
        return "bias"
    if fn.startswith(("thar", "th-ar", "th_ar", "cl_", "cl ", "cathode")):
        return "arc"
    if fn.startswith("dark"):
        return "dark"
    return "science"


def _is_flat_from_context(filename, obj_name):
    """Detect flat frames that have IMAGETYP='Light Frame' but are actually flats."""
    fn = filename.lower()
    obj = obj_name.lower().strip() if obj_name else ""
    # Filename-based detection
    if any(fn.startswith(p) for p in ("flatdome", "flatsky", "domeflat", "skyflat", "flat_", "flat-")):
        return True
    # OBJECT-based: any token is "flat" (catches "flat Halpha", "T080_flat", etc.)
    if obj and "flat" in re.split(r"[\s_\-]+", obj):
        return True
    return False


def normalize_imagetyp(raw, filename="", telescope="", obj_name=""):
    """Normalize IMAGETYP to: science, flat, bias, dark, arc."""
    if telescope == "T152":
        return classify_t152(filename)

    fn = filename.lower()

    # Override: detect flats masquerading as Light Frame
    if _is_flat_from_context(filename, obj_name):
        return "flat"

    if not raw:
        if any(fn.startswith(p) for p in ("bias", "biais", "offset")):
            return "bias"
        if any(fn.startswith(p) for p in ("flat", "ff")):
            return "flat"
        if fn.startswith("dark"):
            return "dark"
        return "science"
    r = raw.lower().strip()
    if r in ("light frame", "light", "object"):
        return "science"
    if r in ("flat field", "flat", "flatfield"):
        return "flat"
    if r in ("bias frame", "bias"):
        return "bias"
    if r in ("dark frame", "dark"):
        return "dark"
    return raw


def normalize_filter(raw):
    """Normalize filter name to canonical form."""
    if not raw:
        return ""
    raw = raw.strip()
    # Direct match
    if raw in FILTER_MAP:
        return FILTER_MAP[raw]
    # Case-insensitive fallback
    raw_lower = raw.lower()
    for key, val in FILTER_MAP.items():
        if raw_lower == key.lower():
            return val
    return raw


def clean_object_name(raw):
    """Clean object name: strip whitespace, normalize underscores to spaces."""
    if not raw:
        return ""
    name = raw.strip()
    name = name.replace("_", " ")
    name = re.sub(r"\s+", " ", name)
    return name


def extract_t152_object(filename):
    """Extract target name from T152 filename.

    2021: HD51530_Mon Mar  1 2021_20.28.40_00022.fits -> HD51530
    2022: HD127614_1200_00039.fits -> HD127614
    """
    fn = os.path.splitext(os.path.basename(filename))[0]
    parts = fn.split("_")
    if parts:
        return parts[0]
    return ""


def extract_iris_object(filename):
    """Extract target from IRIS/ACP filename.

    RAW-M_67-S001-R001-C001-SDSS_g.fits -> M 67
    NGC_5846-S001-R001-C003-SDSS_r.fits -> NGC 5846
    """
    fn = os.path.splitext(os.path.basename(filename))[0]
    # Strip RAW- prefix
    if fn.upper().startswith("RAW-"):
        fn = fn[4:]
    # Match ACP pattern: TARGET-S###-R###-C###-FILTER
    m = re.match(r"^(.+?)-S\d{3}-R\d{3}-C\d{3}", fn)
    if m:
        return m.group(1).replace("_", " ")
    return ""


def safe_float(val, default=""):
    """Convert header value to float string; return default on failure."""
    if not val:
        return default
    try:
        return str(float(val))
    except (ValueError, TypeError):
        return default


def extract_run_info(filepath):
    """Extract (run_label, year) from filepath.

    .../201801OHP_M2/DATA/... -> ('201801OHP_M2', '2018')
    """
    for part in Path(filepath).parts:
        if re.match(r"^\d{6}OHP\d?_M2$", part):
            return part, part[:4]
    return "", ""


# ═══════════════════════════════════════════════════════════════════════════════
# File scanner
# ═══════════════════════════════════════════════════════════════════════════════

def scan_file(filepath):
    """Scan a single FITS file and return a record dict."""
    filepath_str = str(filepath)
    filename = os.path.basename(filepath_str)
    run_label, year = extract_run_info(filepath_str)

    headers = parse_fits_header(filepath_str)
    error = headers.pop("_error", "")

    telescope = identify_telescope(filepath_str, headers)

    # Date
    if telescope == "T152":
        date_obs = headers.get("DATE", "")
    else:
        date_obs = headers.get("DATE-OBS", headers.get("DATE", ""))

    # Exposure time
    if telescope == "T152":
        exptime_raw = headers.get("EXPOSURE", "")
    else:
        exptime_raw = headers.get("EXPTIME", headers.get("EXPOSURE", ""))

    # Object name (extract before imagetyp so we can use it for flat detection)
    obj_raw_header = headers.get("OBJECT", "")
    if telescope == "T152":
        # T152 has no OBJECT header; extract from filename
        obj_raw = extract_t152_object(filename)
    elif telescope == "IRIS" and not obj_raw_header.strip():
        obj_raw = extract_iris_object(filename)
    else:
        obj_raw = obj_raw_header
    object_clean = clean_object_name(obj_raw)

    # Image type (uses obj_raw to detect misclassified flats)
    imagetyp_raw = headers.get("IMAGETYP", "")
    imagetyp = normalize_imagetyp(imagetyp_raw, filename, telescope, obj_name=obj_raw_header)

    # For T152 non-science, clear the object name
    if telescope == "T152" and imagetyp != "science":
        obj_raw = ""
        object_clean = ""
    # For flats masquerading as science, clear the object name
    if imagetyp == "flat" and imagetyp_raw.lower().strip() in ("light frame", "light", "object", ""):
        object_clean = ""  # don't count "FlatDome" as a science target

    # Filter
    filter_raw = headers.get("FILTER", "")
    filter_canonical = normalize_filter(filter_raw)

    # RA / Dec
    ra = headers.get("OBJCTRA", headers.get("RA", ""))
    dec = headers.get("OBJCTDEC", headers.get("DEC", ""))

    return {
        "file_path": filepath_str,
        "filename": filename,
        "year": year,
        "run_label": run_label,
        "telescope": telescope,
        "date_obs": date_obs,
        "object_raw": obj_raw,
        "object_clean": object_clean,
        "imagetyp_raw": imagetyp_raw,
        "imagetyp": imagetyp,
        "filter_raw": filter_raw,
        "filter_canonical": filter_canonical,
        "exptime": safe_float(exptime_raw),
        "ra": ra,
        "dec": dec,
        "airmass": safe_float(headers.get("AIRMASS", "")),
        "naxis1": headers.get("NAXIS1", ""),
        "naxis2": headers.get("NAXIS2", ""),
        "xbinning": headers.get("XBINNING", ""),
        "ybinning": headers.get("YBINNING", ""),
        "ccd_temp": safe_float(headers.get("CCD-TEMP", headers.get("SET-TEMP", ""))),
        "observer": headers.get("OBSERVER", ""),
        "instrume": headers.get("INSTRUME", headers.get("HEAD", "")),
        "simbad_name": "",
        "simbad_ra": "",
        "simbad_dec": "",
        "simbad_type": "",
        "error": error,
    }


def scan_archive():
    """Walk the archive and scan all FITS files. Returns list of record dicts."""
    fits_files = []

    print(f"Scanning archive: {ARCHIVE_ROOT}")
    for root, _dirs, files in os.walk(ARCHIVE_ROOT):
        for f in files:
            if os.path.splitext(f)[1].lower() in FITS_EXT:
                fits_files.append(os.path.join(root, f))

    total = len(fits_files)
    print(f"Found {total} FITS files")

    records = []
    t0 = time.time()
    for i, fp in enumerate(sorted(fits_files)):
        if (i + 1) % 1000 == 0 or i == 0:
            print(f"  Scanning {i + 1}/{total} ({time.time() - t0:.1f}s)...")
        try:
            records.append(scan_file(fp))
        except Exception as e:
            records.append({
                "file_path": fp,
                "filename": os.path.basename(fp),
                "error": str(e),
                **{k: "" for k in CSV_COLUMNS if k not in ("file_path", "filename", "error")},
            })

    elapsed = time.time() - t0
    print(f"Scanned {total} files in {elapsed:.1f}s")
    return records


# ═══════════════════════════════════════════════════════════════════════════════
# Simbad resolution
# ═══════════════════════════════════════════════════════════════════════════════

def load_simbad_cache():
    if SIMBAD_CACHE.exists():
        with open(SIMBAD_CACHE) as f:
            return json.load(f)
    return {}


def save_simbad_cache(cache):
    with open(SIMBAD_CACHE, "w") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)


def query_simbad_tap(name):
    """Query Simbad TAP for a single object name. Returns dict or None."""
    escaped = name.replace("'", "''")
    adql = (
        f"SELECT TOP 1 main_id, ra, dec, otype "
        f"FROM basic JOIN ident ON basic.oid = ident.oidref "
        f"WHERE ident.id = '{escaped}'"
    )
    params = urllib.parse.urlencode({
        "request": "doQuery",
        "lang": "adql",
        "format": "json",
        "query": adql,
    })
    url = f"https://simbad.cds.unistra.fr/simbad/sim-tap/sync?{params}"

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "OHP-Scanner/1.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        if "data" in data and data["data"]:
            row = data["data"][0]
            return {
                "main_id": str(row[0] or ""),
                "ra": round(row[1], 6) if row[1] is not None else "",
                "dec": round(row[2], 6) if row[2] is not None else "",
                "otype": str(row[3] or ""),
            }
    except Exception:
        pass
    return None


def _simbad_variants(name):
    """Generate plausible Simbad name variants for an object."""
    variants = [name]

    # Add/remove space after catalogue prefix:  NGC5846 -> NGC 5846,  NGC 5846 -> NGC5846
    m = re.match(
        r"^(NGC|IC|HD|HIP|HR|SAO|UGC|PGC|Abell|ACO|Mrk|Arp|VV|WASP|HAT-P|TOI|TYC)\s*[-]?\s*(\S+)$",
        name, re.IGNORECASE,
    )
    if m:
        prefix, num = m.group(1), m.group(2)
        variants.append(f"{prefix} {num}")
        variants.append(f"{prefix}{num}")
        variants.append(f"{prefix}-{num}")

    # Messier: M67 -> M 67, M 67 -> M67
    m = re.match(r"^M\s*(\d+)$", name, re.IGNORECASE)
    if m:
        variants.append(f"M {m.group(1)}")
        variants.append(f"M{m.group(1)}")

    # ACO -> Abell
    if name.upper().startswith("ACO"):
        num = name[3:].strip().lstrip("-").strip()
        variants.append(f"Abell {num}")
        variants.append(f"ACO {num}")

    # Hyphens <-> spaces
    if "-" in name:
        variants.append(name.replace("-", " "))
    if " " in name:
        variants.append(name.replace(" ", "-"))

    # HAT-P targets
    m = re.match(r"^HAT[\s-]*P[\s-]*(\d+)", name, re.IGNORECASE)
    if m:
        variants.append(f"HAT-P-{m.group(1)}")

    # TOI targets
    m = re.match(r"^TOI[\s-]*(\d+)", name, re.IGNORECASE)
    if m:
        variants.append(f"TOI-{m.group(1)}")

    # Remove duplicate/original while preserving order
    seen = set()
    unique = []
    for v in variants:
        if v not in seen:
            seen.add(v)
            unique.append(v)
    return unique


def resolve_objects(records):
    """Resolve unique science object names via Simbad TAP; update records in-place."""
    cache = load_simbad_cache()

    # Collect unique science objects
    names = set()
    for r in records:
        if r["imagetyp"] == "science" and r["object_clean"]:
            names.add(r["object_clean"])

    to_query = [n for n in sorted(names) if n not in cache]
    print(f"\nSimbad: {len(names)} unique targets, {len(names) - len(to_query)} cached, {len(to_query)} to query")

    resolved = 0
    failed = 0
    for i, name in enumerate(to_query):
        result = None
        for variant in _simbad_variants(name):
            result = query_simbad_tap(variant)
            if result:
                break
            time.sleep(0.15)

        if result:
            cache[name] = result
            resolved += 1
        else:
            cache[name] = None
            failed += 1

        if (i + 1) % 10 == 0:
            print(f"  Queried {i + 1}/{len(to_query)} ({resolved} resolved, {failed} failed)")
            save_simbad_cache(cache)

        time.sleep(0.15)  # rate limit

    save_simbad_cache(cache)
    print(f"  Done: {resolved} resolved, {failed} unresolved")

    # Apply to records
    for r in records:
        name = r["object_clean"]
        info = cache.get(name)
        if info:
            r["simbad_name"] = info.get("main_id", "")
            r["simbad_ra"] = str(info.get("ra", ""))
            r["simbad_dec"] = str(info.get("dec", ""))
            r["simbad_type"] = otype_label(info.get("otype", ""))


# ═══════════════════════════════════════════════════════════════════════════════
# CSV output
# ═══════════════════════════════════════════════════════════════════════════════

def write_csv(records):
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for r in records:
            writer.writerow(r)
    print(f"Wrote {len(records)} rows → {CSV_PATH}")


# ═══════════════════════════════════════════════════════════════════════════════
# JSON aggregation
# ═══════════════════════════════════════════════════════════════════════════════

def _aggregate_targets(recs):
    """Aggregate science targets from a list of records."""
    targets = defaultdict(lambda: {
        "object": "", "simbad_name": "", "simbad_type": "",
        "simbad_ra": "", "simbad_dec": "",
        "telescopes": set(), "filters": set(),
        "frames": 0, "total_exp": 0.0, "dates": set(),
    })
    for r in recs:
        if r["imagetyp"] != "science" or not r["object_clean"]:
            continue
        obj = r["object_clean"]
        t = targets[obj]
        t["object"] = obj
        t["simbad_name"] = r["simbad_name"] or t["simbad_name"]
        t["simbad_type"] = r["simbad_type"] or t["simbad_type"]
        t["simbad_ra"] = r["simbad_ra"] or t["simbad_ra"]
        t["simbad_dec"] = r["simbad_dec"] or t["simbad_dec"]
        t["telescopes"].add(r["telescope"])
        if r["filter_canonical"]:
            t["filters"].add(r["filter_canonical"])
        t["frames"] += 1
        try:
            t["total_exp"] += float(r["exptime"]) if r["exptime"] else 0
        except ValueError:
            pass
        if r["date_obs"]:
            t["dates"].add(r["date_obs"][:10])
    return targets


def _targets_to_list(targets, include_coords=False):
    """Convert targets dict to sorted JSON-safe list."""
    result = []
    for _obj, t in sorted(targets.items(), key=lambda x: -x[1]["total_exp"]):
        entry = {
            "object": t["object"],
            "simbad_name": t["simbad_name"],
            "simbad_type": t["simbad_type"],
            "telescopes": sorted(t["telescopes"]),
            "filters": sorted(t["filters"]),
            "frames": t["frames"],
            "total_exp": round(t["total_exp"], 1),
            "dates": sorted(t["dates"]),
        }
        if include_coords:
            entry["simbad_ra"] = t["simbad_ra"]
            entry["simbad_dec"] = t["simbad_dec"]
        result.append(entry)
    return result


def build_json(records):
    """Build aggregated JSON for the HTML page."""
    data = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "archive_root": str(ARCHIVE_ROOT),
        "total_files": len(records),
        "total_science": sum(1 for r in records if r["imagetyp"] == "science"),
        "years": sorted(set(r["year"] for r in records if r["year"])),
        "telescopes": {k: v for k, v in TELESCOPES.items()},
        "per_year": {},
        "master_targets": [],
        "filter_summary": [],
        "t152_log": [],
    }

    # ── Per-year ──
    by_year = defaultdict(list)
    for r in records:
        if r["year"]:
            by_year[r["year"]].append(r)

    for year in sorted(by_year):
        recs = by_year[year]
        tel_counts = {}
        for tel in ["T120", "T080", "IRIS", "T152"]:
            tel_recs = [r for r in recs if r["telescope"] == tel]
            if not tel_recs:
                continue
            counts = defaultdict(int)
            for r in tel_recs:
                counts[r["imagetyp"]] += 1
            tel_counts[tel] = dict(counts)
            tel_counts[tel]["total"] = len(tel_recs)

        # Unknown telescope files
        unk_recs = [r for r in recs if r["telescope"] not in TELESCOPES]
        if unk_recs:
            tel_counts["other"] = {"total": len(unk_recs)}

        targets = _aggregate_targets(recs)
        run_labels = sorted(set(r["run_label"] for r in recs if r["run_label"]))

        data["per_year"][year] = {
            "run_labels": run_labels,
            "total_files": len(recs),
            "telescope_counts": tel_counts,
            "targets": _targets_to_list(targets),
        }

    # ── Master target catalogue ──
    all_targets = _aggregate_targets(records)
    # Add years to each target
    for r in records:
        if r["imagetyp"] == "science" and r["object_clean"] and r["year"]:
            all_targets[r["object_clean"]].setdefault("years", set())
            all_targets[r["object_clean"]]["years"].add(r["year"])
    master = []
    for _obj, t in sorted(all_targets.items(), key=lambda x: -x[1]["total_exp"]):
        master.append({
            "object": t["object"],
            "simbad_name": t["simbad_name"],
            "simbad_type": t["simbad_type"],
            "simbad_ra": t["simbad_ra"],
            "simbad_dec": t["simbad_dec"],
            "telescopes": sorted(t["telescopes"]),
            "filters": sorted(t["filters"]),
            "frames": t["frames"],
            "total_exp": round(t["total_exp"], 1),
            "years": sorted(t.get("years", set())),
            "dates": sorted(t["dates"]),
        })
    data["master_targets"] = master

    # ── File index (per-object list of individual files) ──
    file_index = defaultdict(list)
    for r in records:
        if r["imagetyp"] == "science" and r["object_clean"]:
            file_index[r["object_clean"]].append({
                "filename": r["filename"],
                "date": r["date_obs"][:10] if r["date_obs"] else "",
                "telescope": r["telescope"],
                "filter": r["filter_canonical"],
                "exptime": r["exptime"],
                "path": r["file_path"],
            })
    # Sort each object's files by date then filename
    for obj in file_index:
        file_index[obj].sort(key=lambda f: (f["date"], f["filename"]))
    data["file_index"] = dict(file_index)

    # ── Calibration index (T080/T120 only) ──
    cal_index = {}  # year → telescope → binning → {bias, dark, flat}
    for r in records:
        tel = r["telescope"]
        if tel not in ("T080", "T120"):
            continue
        itype = r["imagetyp"]
        if itype not in ("bias", "dark", "flat"):
            continue
        year = r["year"]
        if not year:
            continue
        xb = r.get("xbinning", "")
        yb = r.get("ybinning", "")
        if xb and yb:
            bkey = f"{xb}x{yb}"
        else:
            bkey = "unknown"
        cal_index.setdefault(year, {})
        cal_index[year].setdefault(tel, {})
        cal_index[year][tel].setdefault(bkey, {"bias": [], "dark": [], "flat": {}})
        bucket = cal_index[year][tel][bkey]
        entry = {
            "filename": r["filename"],
            "date": r["date_obs"][:10] if r["date_obs"] else "",
            "exptime": r["exptime"],
            "path": r["file_path"],
        }
        if itype == "flat":
            filt = r["filter_canonical"] or r["filter_raw"] or "unknown"
            bucket["flat"].setdefault(filt, [])
            bucket["flat"][filt].append(entry)
        else:
            bucket[itype].append(entry)
    # Sort each group by date
    for year in cal_index:
        for tel in cal_index[year]:
            for bkey in cal_index[year][tel]:
                bucket = cal_index[year][tel][bkey]
                bucket["bias"].sort(key=lambda f: (f["date"], f["filename"]))
                bucket["dark"].sort(key=lambda f: (f["date"], f["filename"]))
                for filt in bucket["flat"]:
                    bucket["flat"][filt].sort(key=lambda f: (f["date"], f["filename"]))
    data["cal_index"] = cal_index

    # ── Filter summary ──
    fstats = defaultdict(lambda: {"frames": 0, "total_exp": 0.0})
    for r in records:
        filt = r["filter_canonical"] or r["filter_raw"] or "(none)"
        fstats[filt]["frames"] += 1
        try:
            fstats[filt]["total_exp"] += float(r["exptime"]) if r["exptime"] else 0
        except ValueError:
            pass
    data["filter_summary"] = [
        {"filter": f, "frames": s["frames"], "total_exp": round(s["total_exp"], 1)}
        for f, s in sorted(fstats.items(), key=lambda x: -x[1]["frames"])
    ]

    # ── T152 spectroscopy log ──
    t152 = [r for r in records if r["telescope"] == "T152" and r["imagetyp"] == "science"]
    data["t152_log"] = [
        {
            "object": r["object_clean"],
            "simbad_name": r["simbad_name"],
            "date": r["date_obs"],
            "exptime": r["exptime"],
            "filename": r["filename"],
        }
        for r in sorted(t152, key=lambda r: r.get("date_obs", ""))
    ]

    return data


def write_json(data):
    with open(JSON_PATH, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Wrote JSON → {JSON_PATH}")


# ═══════════════════════════════════════════════════════════════════════════════
# HTML generation
# ═══════════════════════════════════════════════════════════════════════════════

def _esc(text):
    """HTML-escape a string."""
    return html_mod.escape(str(text)) if text else ""


def _fmt_exp(seconds):
    """Format exposure time for display."""
    if not seconds:
        return ""
    try:
        s = float(seconds)
    except (ValueError, TypeError):
        return str(seconds)
    if s < 60:
        return f"{s:.0f}s"
    if s < 3600:
        return f"{s / 60:.1f}min"
    return f"{s / 3600:.1f}h"


def generate_html(data, web_mode=False):
    """Generate self-contained dark-themed HTML summary.

    If web_mode is True, file paths in file_index are relative URLs and
    the JS renders them as clickable download links.
    """
    h = _esc  # alias
    out_path = WEB_HTML_PATH if web_mode else HTML_PATH
    title = "OHP M2 Observing Archive" if web_mode else "OHP M2 Archive Summary"

    lines = []
    w = lines.append  # writer shorthand

    # ── Document head ──
    w("<!DOCTYPE html>")
    w('<html lang="en">')
    w("<head>")
    w('<meta charset="UTF-8">')
    w('<meta name="viewport" content="width=device-width, initial-scale=1.0">')
    w(f"<title>{title}</title>")
    w("<style>")
    w(_CSS)
    w("</style>")
    w("</head>")
    w("<body>")

    # ── Header ──
    w('<header class="page-header">')
    w(f"<h1>{h(title)}</h1>")
    w('<p class="subtitle">Observatoire de Haute-Provence &middot; Student Observing Campaigns 2018&ndash;2025</p>')
    w("</header>")

    # ── Dashboard cards ──
    w('<section class="dashboard">')
    years = data["years"]
    tels_used = set()
    for ydata in data["per_year"].values():
        tels_used.update(ydata["telescope_counts"].keys())
    tels_used.discard("other")
    cards = [
        (f"{data['total_files']:,}", "Total Files"),
        (f"{data['total_science']:,}", "Science Frames"),
        (str(len(years)), "Observation Years"),
        (str(len(tels_used)), "Telescopes"),
        (str(len(data['master_targets'])), "Unique Targets"),
    ]
    for val, label in cards:
        w(f'<div class="card"><div class="card-val">{val}</div><div class="card-lbl">{h(label)}</div></div>')
    w("</section>")

    # ── Telescope legend ──
    w('<section class="section">')
    w("<h2>Telescopes</h2>")
    w('<table class="tbl"><thead><tr><th>ID</th><th>Name</th><th>Instrument</th><th>Type</th><th>Description</th></tr></thead><tbody>')
    for tid, tinfo in TELESCOPES.items():
        w(f'<tr><td><strong>{h(tid)}</strong></td><td>{h(tinfo["name"])}</td>'
          f'<td>{h(tinfo["instrument"])}</td><td>{h(tinfo["type"])}</td><td>{h(tinfo["desc"])}</td></tr>')
    w("</tbody></table>")
    w("</section>")

    # ── Per-year controls ──
    w('<section class="section">')
    w('<div class="section-head"><h2>Per-Year Summary</h2>')
    w('<span class="controls"><button onclick="toggleAll(true)">Expand All</button>'
      '<button onclick="toggleAll(false)">Collapse All</button></span></div>')

    for year in sorted(data["per_year"], reverse=True):
        ydata = data["per_year"][year]
        run_str = ", ".join(ydata["run_labels"]) if ydata["run_labels"] else year
        w(f'<details class="year-details"><summary><strong>{h(year)}</strong> &mdash; {h(run_str)} '
          f'({ydata["total_files"]:,} files)</summary>')
        w('<div class="year-body">')

        # Telescope breakdown table
        w("<h3>Telescope Breakdown</h3>")
        types = ["science", "flat", "bias", "dark", "arc"]
        w('<table class="tbl"><thead><tr><th>Telescope</th>')
        for t in types:
            w(f"<th>{t.title()}</th>")
        w("<th>Total</th></tr></thead><tbody>")
        for tel, counts in sorted(ydata["telescope_counts"].items()):
            if tel == "other":
                continue
            w(f"<tr><td><strong>{h(tel)}</strong></td>")
            for t in types:
                v = counts.get(t, 0)
                w(f'<td>{v if v else "&mdash;"}</td>')
            w(f'<td><strong>{counts.get("total", 0)}</strong></td></tr>')
        w("</tbody></table>")

        # ── Calibration Inventory (collapsible, right after telescope breakdown) ──
        targets = ydata["targets"]
        cal_year = data.get("cal_index", {}).get(year, {})
        # Collect science filters per telescope+binning for gap detection
        sci_filters = defaultdict(lambda: defaultdict(set))  # tel -> bin -> set(filter)
        for t in targets:
            for tel_name in t["telescopes"]:
                if tel_name in ("T080", "T120"):
                    for filt in t["filters"]:
                        for bkey in cal_year.get(tel_name, {}):
                            sci_filters[tel_name][bkey].add(filt)
        if cal_year:
            # Count missing flats for summary
            n_missing = 0
            for tel_name in ("T080", "T120"):
                if tel_name not in cal_year:
                    continue
                for bkey in cal_year[tel_name]:
                    bucket = cal_year[tel_name][bkey]
                    if not bucket["bias"]:
                        n_missing += 1
                    sci_for_tb = sci_filters.get(tel_name, {}).get(bkey, set())
                    n_missing += len(sci_for_tb - set(bucket["flat"].keys()))
            warn_badge = f' <span class="cal-warn">{n_missing} missing</span>' if n_missing else ""
            w(f'<details class="cal-details"><summary>Calibration Inventory{warn_badge}</summary>')
            w('<div class="cal-section">')
            for tel_name in ("T080", "T120"):
                if tel_name not in cal_year:
                    continue
                for bkey in sorted(cal_year[tel_name]):
                    bucket = cal_year[tel_name][bkey]
                    n_bias = len(bucket["bias"])
                    n_dark = len(bucket["dark"])
                    flat_filters = sorted(bucket["flat"].keys())
                    w(f"<h4>{h(tel_name)} &mdash; {h(bkey)} binning</h4>")
                    w('<table class="tbl"><thead><tr>'
                      "<th>Type</th><th>Filter</th><th>Count</th>"
                      "<th>Date Range</th><th>Exp Range</th></tr></thead><tbody>")
                    # Bias row
                    if n_bias:
                        dates_b = [f["date"] for f in bucket["bias"] if f["date"]]
                        d_range = f'{dates_b[0]} &ndash; {dates_b[-1]}' if len(dates_b) > 1 else (dates_b[0] if dates_b else "&mdash;")
                        exps_b = sorted(set(f["exptime"] for f in bucket["bias"] if f["exptime"]))
                        e_range = ", ".join(_fmt_exp(e) for e in exps_b[:3]) if exps_b else "&mdash;"
                        cal_key_b = f"{year}_{tel_name}_{bkey}_bias"
                        w(f'<tr><td><strong>Bias</strong></td><td>&mdash;</td>'
                          f'<td><a href="#" class="cal-link" data-cal-type="bias" '
                          f'data-cal-year="{h(year)}" data-cal-tel="{h(tel_name)}" '
                          f'data-cal-bin="{h(bkey)}" data-cal-key="{h(cal_key_b)}" '
                          f'onclick="showCalFiles(this);return false">{n_bias}</a></td>'
                          f'<td class="mono">{d_range}</td><td>{e_range}</td></tr>')
                    else:
                        w('<tr><td><strong>Bias</strong></td><td>&mdash;</td>'
                          '<td>0 <span class="cal-warn severe">missing</span></td>'
                          '<td>&mdash;</td><td>&mdash;</td></tr>')
                    # Flat rows per filter
                    sci_for_tel_bin = sci_filters.get(tel_name, {}).get(bkey, set())
                    for filt in flat_filters:
                        ffiles = bucket["flat"][filt]
                        n_flat = len(ffiles)
                        dates_f = [f["date"] for f in ffiles if f["date"]]
                        d_range = f'{dates_f[0]} &ndash; {dates_f[-1]}' if len(dates_f) > 1 else (dates_f[0] if dates_f else "&mdash;")
                        exps_f = sorted(set(float(f["exptime"]) for f in ffiles if f["exptime"]))
                        e_range = ", ".join(_fmt_exp(e) for e in exps_f[:5]) if exps_f else "&mdash;"
                        cal_key_f = f"{year}_{tel_name}_{bkey}_flat_{filt}"
                        w(f'<tr><td><strong>Flat</strong></td><td>{h(filt)}</td>'
                          f'<td><a href="#" class="cal-link" data-cal-type="flat" '
                          f'data-cal-year="{h(year)}" data-cal-tel="{h(tel_name)}" '
                          f'data-cal-bin="{h(bkey)}" data-cal-filter="{h(filt)}" '
                          f'data-cal-key="{html_mod.escape(cal_key_f, quote=True)}" '
                          f'onclick="showCalFiles(this);return false">{n_flat}</a></td>'
                          f'<td class="mono">{d_range}</td><td>{e_range}</td></tr>')
                    # Gap warnings: science filters with no matching flat
                    missing_flats = sorted(sci_for_tel_bin - set(flat_filters))
                    for mf in missing_flats:
                        w(f'<tr><td><strong>Flat</strong></td><td>{h(mf)}</td>'
                          f'<td>0 <span class="cal-warn">no flat</span></td>'
                          f'<td>&mdash;</td><td>&mdash;</td></tr>')
                    w("</tbody></table>")
            w("</div></details>")

        # Science targets table
        if targets:
            w("<h3>Science Targets</h3>")
            w('<table class="tbl"><thead><tr>'
              "<th>Object</th><th>Simbad Name</th><th>Type</th>"
              "<th>Telescope</th><th>Filters</th><th>Frames</th>"
              "<th>Total Exp</th><th>Dates</th></tr></thead><tbody>")
            for t in targets:
                obj_esc = html_mod.escape(t["object"], quote=True)
                w(f'<tr><td><a href="#" class="obj-link" data-obj="{obj_esc}" onclick="showFiles(this);return false">{h(t["object"])}</a></td>'
                  f'<td>{h(t["simbad_name"])}</td>'
                  f'<td>{h(t["simbad_type"])}</td>'
                  f'<td>{", ".join(t["telescopes"])}</td>'
                  f'<td>{", ".join(t["filters"])}</td>'
                  f'<td>{t["frames"]}</td>'
                  f'<td>{_fmt_exp(t["total_exp"])}</td>'
                  f'<td class="mono">{", ".join(t["dates"])}</td></tr>')
            w("</tbody></table>")

        w("</div></details>")

    w("</section>")

    # ── Master target catalogue ──
    w('<section class="section" id="master-section">')
    w("<h2>Master Target Catalogue</h2>")
    # Sticky filter bar
    all_tels = sorted(tels_used)
    all_years = sorted(data["per_year"].keys(), reverse=True)
    w('<div class="filter-bar">')
    w('<input type="text" id="master-search" class="search-box" '
      'placeholder="Search targets... ( / )" oninput="filterMaster()">')
    w('<select id="filter-telescope" onchange="filterMaster()">')
    w('<option value="">All telescopes</option>')
    for tel in all_tels:
        w(f'<option value="{h(tel)}">{h(tel)}</option>')
    w('</select>')
    w('<select id="filter-year" onchange="filterMaster()">')
    w('<option value="">All years</option>')
    for yr in all_years:
        w(f'<option value="{h(yr)}">{h(yr)}</option>')
    w('</select>')
    w('<span class="filter-count" id="filter-count"></span>')
    w('</div>')
    w('<table class="tbl sortable" id="master-table">')
    w("<thead><tr>")
    cols = ["Object", "Simbad Name", "Type", "RA", "Dec",
            "Telescopes", "Filters", "Frames", "Total Exp", "Years"]
    for i, col in enumerate(cols):
        w(f'<th onclick="sortTable(\'master-table\',{i})" class="sortable-th">{col}</th>')
    w("</tr></thead><tbody>")
    for t in data["master_targets"]:
        ra_str = f'{t["simbad_ra"]:.4f}' if isinstance(t.get("simbad_ra"), (int, float)) and t["simbad_ra"] != "" else h(t.get("simbad_ra", ""))
        dec_str = f'{t["simbad_dec"]:.4f}' if isinstance(t.get("simbad_dec"), (int, float)) and t["simbad_dec"] != "" else h(t.get("simbad_dec", ""))
        obj_esc = html_mod.escape(t["object"], quote=True)
        w(f'<tr><td><a href="#" class="obj-link" data-obj="{obj_esc}" onclick="showFiles(this);return false">{h(t["object"])}</a></td>'
          f'<td>{h(t["simbad_name"])}</td>'
          f'<td>{h(t["simbad_type"])}</td>'
          f'<td class="mono">{ra_str}</td>'
          f'<td class="mono">{dec_str}</td>'
          f'<td>{", ".join(t["telescopes"])}</td>'
          f'<td>{", ".join(t["filters"])}</td>'
          f'<td>{t["frames"]}</td>'
          f'<td data-sort="{t["total_exp"]}">{_fmt_exp(t["total_exp"])}</td>'
          f'<td>{", ".join(t.get("years", []))}</td></tr>')
    w("</tbody></table>")
    w("</section>")

    # ── Filter usage summary ──
    w('<section class="section">')
    w("<h2>Filter Usage</h2>")
    w('<table class="tbl"><thead><tr><th>Filter</th><th>Frames</th><th>Total Exposure</th></tr></thead><tbody>')
    for f in data["filter_summary"]:
        w(f'<tr><td><strong>{h(f["filter"])}</strong></td>'
          f'<td>{f["frames"]:,}</td>'
          f'<td>{_fmt_exp(f["total_exp"])}</td></tr>')
    w("</tbody></table>")
    w("</section>")

    # ── T152 spectroscopy log ──
    if data["t152_log"]:
        w('<section class="section">')
        w("<h2>T152 Spectroscopy Log</h2>")
        w('<table class="tbl"><thead><tr><th>Object</th><th>Simbad Name</th>'
          "<th>Date</th><th>Exposure</th><th>Filename</th></tr></thead><tbody>")
        for r in data["t152_log"]:
            w(f'<tr><td>{h(r["object"])}</td>'
              f'<td>{h(r["simbad_name"])}</td>'
              f'<td class="mono">{h(r["date"])}</td>'
              f'<td>{_fmt_exp(r["exptime"])}</td>'
              f'<td class="mono">{h(r["filename"])}</td></tr>')
        w("</tbody></table>")
        w("</section>")

    # ── Footer ──
    w("<footer>")
    if web_mode:
        w(f'<p>Generated: {h(data["generated"])} &middot; '
          f'{data["total_files"]:,} files scanned</p>')
    else:
        w(f'<p>Archive: <code>{h(str(ARCHIVE_ROOT))}</code> &middot; '
          f'Generated: {h(data["generated"])} &middot; '
          f'{data["total_files"]:,} files scanned</p>')
    w('<p class="kbd-hint">Keyboard: '
      '<kbd>/</kbd> search &middot; '
      '<kbd>1</kbd>&ndash;<kbd>8</kbd> jump to year &middot; '
      '<kbd>Esc</kbd> close panels</p>')
    w("</footer>")

    # ── Embedded JSON + JavaScript ──
    w("<script>")
    w(f"const DATA = {json.dumps(data, ensure_ascii=False)};")
    w(f"const WEB_MODE = {'true' if web_mode else 'false'};")
    w(_JS)
    w("</script>")

    w("</body></html>")

    html_str = "\n".join(lines)
    with open(out_path, "w") as f:
        f.write(html_str)
    print(f"Wrote HTML → {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# CSS
# ═══════════════════════════════════════════════════════════════════════════════

_CSS = """\
:root {
  --bg: #0d1117; --surface: #161b22; --border: #30363d;
  --text: #c9d1d9; --text-dim: #8b949e; --heading: #f0f6fc;
  --accent: #58a6ff; --green: #3fb950; --orange: #d29922; --red: #f85149;
  --font: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
  --mono: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
}
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
  background: var(--bg); color: var(--text); font-family: var(--font);
  font-size: 14px; line-height: 1.6; max-width: 1400px;
  margin: 0 auto; padding: 20px 24px 60px;
}
a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }
h1, h2, h3 { color: var(--heading); font-weight: 600; }
h1 { font-size: 1.8em; }
h2 { font-size: 1.3em; margin-bottom: 12px; border-bottom: 1px solid var(--border); padding-bottom: 6px; }
h3 { font-size: 1.05em; margin: 16px 0 8px; color: var(--text); }
code { font-family: var(--mono); font-size: 0.92em; background: var(--surface); padding: 2px 6px; border-radius: 4px; }
.mono { font-family: var(--mono); font-size: 0.88em; }

/* Header */
.page-header { text-align: center; margin-bottom: 32px; padding: 24px 0 16px; border-bottom: 1px solid var(--border); }
.subtitle { color: var(--text-dim); margin-top: 4px; }

/* Dashboard */
.dashboard {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 16px; margin-bottom: 32px;
}
.card {
  background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
  padding: 20px 16px; text-align: center;
}
.card-val { font-size: 2em; font-weight: 700; color: var(--accent); }
.card-lbl { color: var(--text-dim); font-size: 0.85em; margin-top: 4px; }

/* Sections */
.section { margin-bottom: 32px; }
.section-head { display: flex; justify-content: space-between; align-items: baseline; flex-wrap: wrap; gap: 8px; }
.controls button {
  background: var(--surface); color: var(--text); border: 1px solid var(--border);
  border-radius: 6px; padding: 4px 12px; cursor: pointer; font-size: 0.85em;
}
.controls button:hover { border-color: var(--accent); color: var(--accent); }

/* Tables */
.tbl {
  width: 100%; border-collapse: collapse; margin-bottom: 12px;
  font-size: 0.9em;
}
.tbl th, .tbl td {
  padding: 6px 10px; text-align: left; border-bottom: 1px solid var(--border);
  white-space: nowrap;
}
.tbl td:last-child, .tbl th:last-child { white-space: normal; }
.tbl thead th { background: var(--surface); color: var(--heading); font-weight: 600; position: sticky; top: 0; }
.tbl tbody tr:hover { background: rgba(88,166,255,0.06); }
.sortable-th { cursor: pointer; user-select: none; }
.sortable-th:hover { color: var(--accent); }
.sortable-th::after { content: ' \\2195'; color: var(--text-dim); font-size: 0.8em; }

/* Year accordion */
.year-details { margin-bottom: 8px; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }
.year-details summary {
  padding: 10px 16px; background: var(--surface); cursor: pointer;
  list-style: none; display: flex; align-items: center; gap: 8px;
}
.year-details summary::-webkit-details-marker { display: none; }
.year-details summary::before { content: '\\25B6'; font-size: 0.7em; color: var(--text-dim); transition: transform 0.2s; }
.year-details[open] summary::before { transform: rotate(90deg); }
.year-details summary:hover { background: rgba(88,166,255,0.08); }
.year-body { padding: 12px 16px; }
.year-body .tbl { font-size: 0.85em; }

/* Search */
.search-box {
  width: 100%; max-width: 400px; padding: 8px 12px; margin-bottom: 12px;
  background: var(--surface); color: var(--text); border: 1px solid var(--border);
  border-radius: 6px; font-size: 0.9em; outline: none;
}
.search-box:focus { border-color: var(--accent); }

/* Footer */
footer {
  margin-top: 48px; padding-top: 16px; border-top: 1px solid var(--border);
  color: var(--text-dim); font-size: 0.85em; text-align: center;
}

/* Object file panel */
.obj-link { cursor: pointer; }
.file-panel td { padding: 0 !important; border-bottom: none !important; }
.file-panel-inner {
  background: var(--surface); border: 1px solid var(--border); border-radius: 6px;
  margin: 4px 8px 8px; padding: 10px 14px; overflow-x: auto;
}
.file-panel-inner table {
  width: 100%; border-collapse: collapse; font-family: var(--mono); font-size: 0.82em;
}
.file-panel-inner th {
  text-align: left; padding: 3px 8px; color: var(--text-dim);
  border-bottom: 1px solid var(--border); font-weight: 600;
}
.file-panel-inner td { padding: 2px 8px; color: var(--text); }
.file-panel-inner tr:hover td { background: rgba(88,166,255,0.06); }
.file-panel-inner .fp-path { user-select: all; }
.fp-toolbar {
  display: flex; align-items: center; gap: 10px; margin-bottom: 6px; padding: 4px 0;
}
.fp-toolbar label { font-size: 0.85em; color: var(--text-dim); cursor: pointer; user-select: none; }
.fp-toolbar label input { margin-right: 4px; vertical-align: middle; }
.fp-toolbar .fp-sel-count { font-size: 0.82em; color: var(--text-dim); }
.fp-toolbar button {
  background: var(--accent); color: #fff; border: none; border-radius: 4px;
  padding: 4px 12px; font-size: 0.82em; font-weight: 600; cursor: pointer;
}
.fp-toolbar button:disabled { opacity: 0.4; cursor: default; }
.fp-toolbar button:hover:not(:disabled) { opacity: 0.85; }
.fp-cb { width: 15px; text-align: center; }
.fp-cb input { cursor: pointer; }
.dl-modal-overlay {
  position: fixed; top: 0; left: 0; width: 100%; height: 100%;
  background: rgba(0,0,0,0.6); z-index: 1000; display: flex; align-items: center; justify-content: center;
}
.dl-modal {
  background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
  padding: 20px; max-width: 700px; width: 90%; max-height: 80vh; display: flex; flex-direction: column;
}
.dl-modal h3 { margin: 0 0 8px; color: var(--text); font-size: 1em; }
.dl-modal textarea {
  width: 100%; height: 300px; background: var(--bg); color: var(--green); border: 1px solid var(--border);
  border-radius: 4px; padding: 10px; font-family: var(--mono); font-size: 0.82em;
  resize: vertical; white-space: pre; overflow-wrap: normal; overflow-x: auto;
}
.dl-modal .dl-actions { display: flex; gap: 8px; margin-top: 10px; justify-content: flex-end; }
.dl-modal .dl-actions button {
  background: var(--accent); color: #fff; border: none; border-radius: 4px;
  padding: 6px 16px; font-size: 0.85em; font-weight: 600; cursor: pointer;
}
.dl-modal .dl-actions button.dl-close { background: var(--surface); color: var(--text); border: 1px solid var(--border); }

/* Calibration inventory */
.cal-details { margin-top: 12px; margin-bottom: 8px; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }
.cal-details summary { cursor: pointer; padding: 8px 14px; font-weight: 600; font-size: 0.95em; background: rgba(22,27,34,0.5); user-select: none; list-style: none; }
.cal-details summary::-webkit-details-marker { display: none; }
.cal-details summary::before { content: '\25B6'; font-size: 0.7em; color: var(--text-dim); margin-right: 8px; display: inline-block; transition: transform 0.2s; }
.cal-details[open] summary::before { transform: rotate(90deg); }
.cal-details summary:hover { background: rgba(88,166,255,0.08); }
.cal-section { padding: 8px 16px 12px; }
.cal-section .tbl { margin-bottom: 4px; }
.cal-warn { display: inline-block; background: var(--orange); color: #000; font-size: 0.75em; font-weight: 700; padding: 1px 7px; border-radius: 10px; margin-left: 4px; vertical-align: middle; }
.cal-warn.severe { background: var(--red); color: #fff; }
.cal-link { cursor: pointer; color: var(--accent); }
.cal-link:hover { text-decoration: underline; }
.cal-panel td { padding: 0 !important; border-bottom: none !important; }
.cal-panel-inner {
  background: var(--bg); border: 1px solid var(--border); border-radius: 6px;
  margin: 4px 8px 8px; padding: 10px 14px; overflow-x: auto;
}
.cal-panel-inner table { width: 100%; border-collapse: collapse; font-family: var(--mono); font-size: 0.82em; }
.cal-panel-inner th { text-align: left; padding: 3px 8px; color: var(--text-dim); border-bottom: 1px solid var(--border); font-weight: 600; }
.cal-panel-inner td { padding: 2px 8px; color: var(--text); }
.cal-panel-inner tr:hover td { background: rgba(88,166,255,0.06); }
.cal-summary { font-size: 0.85em; color: var(--text-dim); margin-bottom: 6px; padding: 6px 8px; background: rgba(22,27,34,0.6); border-radius: 4px; border-left: 3px solid var(--accent); }
.cal-summary .cal-ok { color: var(--green); }
.cal-summary .cal-missing { color: var(--orange); font-weight: 600; }

/* Sticky filter bar */
.filter-bar {
  position: sticky; top: 0; z-index: 100;
  background: var(--bg); border-bottom: 1px solid var(--border);
  padding: 8px 0; margin: -8px 0 12px;
  display: flex; align-items: center; gap: 10px; flex-wrap: wrap;
}
.filter-bar .search-box { margin-bottom: 0; max-width: 280px; }
.filter-bar select {
  background: var(--surface); color: var(--text); border: 1px solid var(--border);
  border-radius: 6px; padding: 7px 10px; font-size: 0.9em; outline: none; cursor: pointer;
}
.filter-bar select:focus { border-color: var(--accent); }
.filter-bar .filter-count { color: var(--text-dim); font-size: 0.85em; margin-left: auto; }
.kbd-hint {
  color: var(--text-dim); font-size: 0.8em; margin-top: 4px;
}
.kbd-hint kbd {
  display: inline-block; background: var(--surface); border: 1px solid var(--border);
  border-radius: 3px; padding: 1px 5px; font-family: var(--mono); font-size: 0.9em;
  color: var(--text); margin: 0 1px;
}

/* Responsive */
@media (max-width: 768px) {
  body { padding: 12px; }
  .tbl { display: block; overflow-x: auto; }
  .dashboard { grid-template-columns: repeat(2, 1fr); }
  .filter-bar { flex-direction: column; align-items: stretch; }
  .filter-bar .search-box { max-width: 100%; }
  .filter-bar .filter-count { margin-left: 0; }
}
"""

# ═══════════════════════════════════════════════════════════════════════════════
# JavaScript
# ═══════════════════════════════════════════════════════════════════════════════

_JS = """\
function toggleAll(expand) {
  document.querySelectorAll('details.year-details').forEach(d => d.open = expand);
}

function sortTable(tableId, colIdx) {
  const table = document.getElementById(tableId);
  const tbody = table.querySelector('tbody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  const th = table.querySelectorAll('thead th')[colIdx];

  // Determine sort direction
  const curDir = th.dataset.sortDir || 'none';
  // Reset all headers
  table.querySelectorAll('thead th').forEach(h => h.dataset.sortDir = 'none');
  const newDir = curDir === 'asc' ? 'desc' : 'asc';
  th.dataset.sortDir = newDir;

  rows.sort((a, b) => {
    let aCell = a.cells[colIdx];
    let bCell = b.cells[colIdx];
    let aVal = aCell.dataset.sort !== undefined ? aCell.dataset.sort : aCell.textContent.trim();
    let bVal = bCell.dataset.sort !== undefined ? bCell.dataset.sort : bCell.textContent.trim();

    // Try numeric comparison
    let aNum = parseFloat(aVal);
    let bNum = parseFloat(bVal);
    if (!isNaN(aNum) && !isNaN(bNum)) {
      return newDir === 'asc' ? aNum - bNum : bNum - aNum;
    }
    // String comparison
    let cmp = aVal.localeCompare(bVal, undefined, {sensitivity: 'base'});
    return newDir === 'asc' ? cmp : -cmp;
  });

  rows.forEach(r => tbody.appendChild(r));
}

function filterMaster() {
  const query = document.getElementById('master-search').value.toLowerCase();
  const telFilter = document.getElementById('filter-telescope').value;
  const yearFilter = document.getElementById('filter-year').value;
  const tbody = document.querySelector('#master-table tbody');
  let shown = 0, total = 0;
  tbody.querySelectorAll('tr').forEach(row => {
    if (row.classList.contains('file-panel')) { row.remove(); return; }
    total++;
    const text = row.textContent.toLowerCase();
    const telCol = row.cells[5] ? row.cells[5].textContent : '';
    const yearCol = row.cells[9] ? row.cells[9].textContent : '';
    const matchText = !query || text.includes(query);
    const matchTel = !telFilter || telCol.includes(telFilter);
    const matchYear = !yearFilter || yearCol.includes(yearFilter);
    const vis = matchText && matchTel && matchYear;
    row.style.display = vis ? '' : 'none';
    if (vis) shown++;
  });
  const countEl = document.getElementById('filter-count');
  if (countEl) {
    countEl.textContent = (query || telFilter || yearFilter) ? shown + ' / ' + total : total + ' targets';
  }
}

// Keyboard shortcuts
document.addEventListener('keydown', function(e) {
  // Ignore when typing in an input
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT' || e.target.tagName === 'TEXTAREA') {
    if (e.key === 'Escape') { e.target.blur(); return; }
    return;
  }

  // / → focus search
  if (e.key === '/') {
    e.preventDefault();
    const search = document.getElementById('master-search');
    if (search) { search.focus(); search.scrollIntoView({behavior: 'smooth', block: 'center'}); }
    return;
  }

  // Esc → close any open file/cal panels
  if (e.key === 'Escape') {
    document.querySelectorAll('.file-panel, .cal-panel').forEach(p => p.remove());
    return;
  }

  // 1-8 → jump to year accordion (1 = most recent)
  const digit = parseInt(e.key);
  if (digit >= 1 && digit <= 8) {
    const details = document.querySelectorAll('details.year-details');
    const idx = digit - 1;
    if (idx < details.length) {
      details[idx].open = true;
      details[idx].scrollIntoView({behavior: 'smooth', block: 'start'});
    }
    return;
  }
});

function getCalSummary(obj, year) {
  // Build calibration summary for a science object in a given year (or all years)
  const files = (DATA.file_index || {})[obj];
  if (!files || !DATA.cal_index) return '';
  // Determine telescope+binning combos used and science filters
  const combos = {};  // "tel|bin" -> Set of filters
  for (const f of files) {
    if (year && f.date && !f.date.startsWith(year)) continue;
    const tel = f.telescope;
    if (tel !== 'T080' && tel !== 'T120') continue;
    // Look up binning from cal_index structure (files don't carry binning)
    // We'll check all binnings for this tel/year
    const y = f.date ? f.date.substring(0, 4) : '';
    if (!y) continue;
    const key = tel + '|' + y;
    if (!combos[key]) combos[key] = new Set();
    if (f.filter) combos[key].add(f.filter);
  }
  if (Object.keys(combos).length === 0) return '';

  const parts = [];
  for (const key of Object.keys(combos).sort()) {
    const [tel, y] = key.split('|');
    const sciFilters = combos[key];
    const calYear = (DATA.cal_index || {})[y];
    if (!calYear || !calYear[tel]) {
      parts.push('<span class="cal-missing">No calibrations for ' + esc(tel) + ' ' + esc(y) + '</span>');
      continue;
    }
    // Check each binning available
    for (const [bkey, bucket] of Object.entries(calYear[tel])) {
      const flatFilters = Object.keys(bucket.flat || {});
      const flatParts = flatFilters.map(f => bucket.flat[f].length + '\\u00d7 ' + esc(f)).join(', ');
      let line = '<b>' + esc(tel) + ' ' + esc(bkey) + ' (' + esc(y) + ')</b>: ';
      line += '<span class="cal-ok">' + bucket.bias.length + ' bias</span>';
      if (flatParts) line += ', ' + flatParts;
      // Gap detection
      const missing = [];
      for (const sf of sciFilters) {
        if (!bucket.flat[sf] || bucket.flat[sf].length === 0) missing.push(sf);
      }
      if (missing.length > 0) {
        line += ' <span class="cal-missing">[no flat: ' + missing.map(esc).join(', ') + ']</span>';
      }
      parts.push(line);
    }
  }
  return parts.length ? '<div class="cal-summary">' + parts.join('<br>') + '</div>' : '';
}

function showFiles(el) {
  const obj = el.getAttribute('data-obj');
  const row = el.closest('tr');
  const next = row.nextElementSibling;

  // Toggle off if already open for this object
  if (next && next.classList.contains('file-panel') && next.dataset.obj === obj) {
    next.remove();
    return;
  }
  // Remove any other open panel in the same table
  const tbody = row.closest('tbody');
  const existing = tbody.querySelector('tr.file-panel');
  if (existing) existing.remove();

  const files = (DATA.file_index || {})[obj];
  if (!files || files.length === 0) return;

  // Determine year context from enclosing year-details if any
  const yearDetails = row.closest('.year-details');
  let year = null;
  if (yearDetails) {
    const summary = yearDetails.querySelector('summary strong');
    if (summary) year = summary.textContent.trim();
  }

  const cols = row.closest('table').querySelector('thead tr').cells.length;
  const panel = document.createElement('tr');
  panel.className = 'file-panel';
  panel.dataset.obj = obj;
  const td = document.createElement('td');
  td.colSpan = cols;

  // Calibration summary at top
  const calSummary = getCalSummary(obj, year);

  // Filter files by year
  const filtered = files.filter(f => !(year && f.date && !f.date.startsWith(year)));

  const cbCol = WEB_MODE ? '<th class="fp-cb"><input type="checkbox" title="Select all" onchange="toggleFileCheckboxes(this)"></th>' : '';
  const toolbar = WEB_MODE
    ? '<div class="fp-toolbar"><label><input type="checkbox" onchange="toggleFileCheckboxes(this)"> Select all</label>' +
      '<span class="fp-sel-count"></span>' +
      '<button onclick="generateScript(this)" disabled>Download script</button></div>'
    : '';

  let html = calSummary + '<div class="file-panel-inner">' + toolbar + '<table><thead><tr>' +
    cbCol + '<th>Filename</th><th>Date</th><th>Telescope</th><th>Filter</th><th>Exposure</th><th>Path</th>' +
    '</tr></thead><tbody>';
  for (const f of filtered) {
    const exp = parseFloat(f.exptime);
    const expStr = isNaN(exp) ? f.exptime : (exp >= 3600 ? (exp/3600).toFixed(1)+'h' : exp >= 60 ? (exp/60).toFixed(1)+'m' : exp+'s');
    const pathCell = WEB_MODE
      ? '<td class="fp-path"><a href="' + esc(f.path) + '" download>' + esc(f.path) + '</a></td>'
      : '<td class="fp-path">' + esc(f.path) + '</td>';
    const cb = WEB_MODE ? '<td class="fp-cb"><input type="checkbox" data-url="' + esc(f.path) + '" onchange="updateSelCount(this)"></td>' : '';
    html += '<tr>' + cb +
      '<td>' + esc(f.filename) + '</td>' +
      '<td>' + esc(f.date) + '</td>' +
      '<td>' + esc(f.telescope) + '</td>' +
      '<td>' + esc(f.filter) + '</td>' +
      '<td>' + expStr + '</td>' +
      pathCell + '</tr>';
  }
  html += '</tbody></table></div>';
  td.innerHTML = html;
  panel.appendChild(td);
  row.after(panel);
}

function showCalFiles(el) {
  const calType = el.getAttribute('data-cal-type');
  const calYear = el.getAttribute('data-cal-year');
  const calTel = el.getAttribute('data-cal-tel');
  const calBin = el.getAttribute('data-cal-bin');
  const calFilter = el.getAttribute('data-cal-filter') || null;
  const row = el.closest('tr');
  const next = row.nextElementSibling;

  // Toggle off
  if (next && next.classList.contains('cal-panel') && next.dataset.calKey === el.dataset.calKey) {
    next.remove();
    return;
  }
  // Remove any other open cal panel in the same table
  const tbody = row.closest('tbody');
  const existing = tbody.querySelector('tr.cal-panel');
  if (existing) existing.remove();

  const bucket = ((DATA.cal_index || {})[calYear] || {})[calTel];
  if (!bucket || !bucket[calBin]) return;
  const data = bucket[calBin];
  let files;
  if (calFilter) {
    files = (data.flat || {})[calFilter] || [];
  } else {
    files = data[calType] || [];
  }
  if (files.length === 0) return;

  const cols = row.closest('table').querySelector('thead tr').cells.length;
  const panel = document.createElement('tr');
  panel.className = 'cal-panel';
  panel.dataset.calKey = el.dataset.calKey;
  const td = document.createElement('td');
  td.colSpan = cols;

  const cbCol = WEB_MODE ? '<th class="fp-cb"><input type="checkbox" title="Select all" onchange="toggleFileCheckboxes(this)"></th>' : '';
  const toolbar = WEB_MODE
    ? '<div class="fp-toolbar"><label><input type="checkbox" onchange="toggleFileCheckboxes(this)"> Select all</label>' +
      '<span class="fp-sel-count"></span>' +
      '<button onclick="generateScript(this)" disabled>Download script</button></div>'
    : '';

  let html = '<div class="cal-panel-inner">' + toolbar + '<table><thead><tr>' +
    cbCol + '<th>Filename</th><th>Date</th><th>Exposure</th><th>Path</th>' +
    '</tr></thead><tbody>';
  for (const f of files) {
    const exp = parseFloat(f.exptime);
    const expStr = isNaN(exp) ? (f.exptime || '') : (exp >= 3600 ? (exp/3600).toFixed(1)+'h' : exp >= 60 ? (exp/60).toFixed(1)+'m' : exp+'s');
    const pathCell = WEB_MODE
      ? '<td class="fp-path"><a href="' + esc(f.path) + '" download>' + esc(f.path) + '</a></td>'
      : '<td class="fp-path">' + esc(f.path) + '</td>';
    const cb = WEB_MODE ? '<td class="fp-cb"><input type="checkbox" data-url="' + esc(f.path) + '" onchange="updateSelCount(this)"></td>' : '';
    html += '<tr>' + cb + '<td>' + esc(f.filename) + '</td><td>' + esc(f.date) + '</td><td>' + expStr + '</td>' + pathCell + '</tr>';
  }
  html += '</tbody></table></div>';
  td.innerHTML = html;
  panel.appendChild(td);
  row.after(panel);
}

function toggleFileCheckboxes(src) {
  const panel = src.closest('.file-panel-inner, .cal-panel-inner');
  if (!panel) return;
  const checked = src.checked;
  // Sync all checkboxes (header th, toolbar label, and row checkboxes)
  panel.querySelectorAll('input[type="checkbox"]').forEach(cb => cb.checked = checked);
  updateSelCount(src);
}

function updateSelCount(el) {
  const panel = el.closest('.file-panel-inner, .cal-panel-inner');
  if (!panel) return;
  const boxes = panel.querySelectorAll('tbody input[type="checkbox"]');
  const checked = panel.querySelectorAll('tbody input[type="checkbox"]:checked');
  const n = checked.length;
  const countEl = panel.querySelector('.fp-sel-count');
  if (countEl) countEl.textContent = n > 0 ? n + ' of ' + boxes.length + ' selected' : '';
  const btn = panel.querySelector('.fp-toolbar button');
  if (btn) btn.disabled = (n === 0);
  // Sync select-all checkboxes
  const allChecked = n === boxes.length && n > 0;
  panel.querySelectorAll('.fp-toolbar input[type="checkbox"], thead .fp-cb input[type="checkbox"]')
    .forEach(cb => cb.checked = allChecked);
}

function generateScript(btn) {
  const panel = btn.closest('.file-panel-inner, .cal-panel-inner');
  if (!panel) return;
  const urls = [];
  panel.querySelectorAll('tbody input[type="checkbox"]:checked').forEach(cb => {
    const url = cb.getAttribute('data-url');
    if (url) urls.push(url);
  });
  if (urls.length === 0) return;

  const script = '#!/bin/bash\\n# Download ' + urls.length + ' file' + (urls.length > 1 ? 's' : '') +
    ' from OHP archive\\n# Generated by ohp-m2-archive\\n\\n' +
    urls.map(u => 'wget -nc "' + u + '"').join('\\n') + '\\n';

  // Show modal
  const overlay = document.createElement('div');
  overlay.className = 'dl-modal-overlay';
  overlay.innerHTML =
    '<div class="dl-modal">' +
    '<h3>Download script (' + urls.length + ' file' + (urls.length > 1 ? 's' : '') + ')</h3>' +
    '<textarea readonly>' + script.replace(/&/g,'&amp;').replace(/</g,'&lt;') + '</textarea>' +
    '<div class="dl-actions">' +
    '<button class="dl-close" onclick="this.closest(\\'.dl-modal-overlay\\').remove()">Close</button>' +
    '<button onclick="copyScript(this)">Copy to clipboard</button>' +
    '<button onclick="saveScript(this)">Save as .sh</button>' +
    '</div></div>';
  overlay.addEventListener('click', function(e) { if (e.target === overlay) overlay.remove(); });
  document.body.appendChild(overlay);
}

function copyScript(btn) {
  const ta = btn.closest('.dl-modal').querySelector('textarea');
  navigator.clipboard.writeText(ta.value).then(() => {
    btn.textContent = 'Copied!';
    setTimeout(() => btn.textContent = 'Copy to clipboard', 1500);
  });
}

function saveScript(btn) {
  const ta = btn.closest('.dl-modal').querySelector('textarea');
  const blob = new Blob([ta.value], {type: 'text/x-shellscript'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'download_ohp.sh';
  a.click();
  URL.revokeObjectURL(a.href);
}

function esc(s) {
  if (!s) return '';
  const d = document.createElement('div');
  d.textContent = s;
  return d.innerHTML;
}

// Init: show target count
document.addEventListener('DOMContentLoaded', function() { filterMaster(); });
"""


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()

    # 1. Scan all FITS files
    records = scan_archive()

    # 2. Resolve science targets via Simbad
    resolve_objects(records)

    # 3. Write CSV
    write_csv(records)

    # 4. Build and write JSON
    data = build_json(records)
    write_json(data)

    # 5. Generate HTML (local version)
    generate_html(data)

    # 6. Generate web HTML (absolute URLs for external hosting)
    web_data = copy.deepcopy(data)
    web_data["archive_root"] = WEB_BASE_URL
    archive_prefix = str(ARCHIVE_ROOT) + "/"

    def _to_web_url(local_path):
        """Convert local archive path to remote URL, applying run name fixes."""
        if not local_path.startswith(archive_prefix):
            return local_path
        rel = local_path[len(archive_prefix):]
        for local_name, remote_name in RUN_NAME_MAP.items():
            if rel.startswith(local_name + "/") or rel == local_name:
                rel = remote_name + rel[len(local_name):]
                break
        return WEB_BASE_URL + rel

    for obj, files in web_data.get("file_index", {}).items():
        for f in files:
            f["path"] = _to_web_url(f["path"])
    for year_ci in web_data.get("cal_index", {}).values():
        for tel_ci in year_ci.values():
            for bin_ci in tel_ci.values():
                for cal_type in ("bias", "dark"):
                    for f in bin_ci.get(cal_type, []):
                        f["path"] = _to_web_url(f["path"])
                for filt_files in bin_ci.get("flat", {}).values():
                    for f in filt_files:
                        f["path"] = _to_web_url(f["path"])
    generate_html(web_data, web_mode=True)

    # 7. Copy web HTML to docs/ for GitHub Pages
    DOCS_DIR.mkdir(exist_ok=True)
    docs_index = DOCS_DIR / "index.html"
    shutil.copy2(WEB_HTML_PATH, docs_index)
    print(f"Copied → {docs_index}")

    elapsed = time.time() - t_start
    print(f"\nDone in {elapsed:.1f}s — {len(records):,} files processed")
    print(f"  CSV:  {CSV_PATH}")
    print(f"  JSON: {JSON_PATH}")
    print(f"  HTML: {HTML_PATH}")
    print(f"  Web:  {WEB_HTML_PATH}")
    print(f"  Docs: {docs_index}")


if __name__ == "__main__":
    main()
