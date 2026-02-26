"""Calibration: master bias, master flat, science frame reduction."""

import logging
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import ccdproc
from astropy.nddata import CCDData
import astropy.units as u

from . import config
from .utils import read_ccd, write_ccd, output_dir, write_json, get_science_frames
from .quality import check_bias_frames, check_flat_frames

log = logging.getLogger("pipeline")


def make_master_bias(bias_entries: list, year: str, telescope: str,
                     force: bool = False) -> tuple[CCDData | None, dict]:
    """Build master bias from QC-checked frames.

    Returns (master_bias_ccd, qc_report).
    """
    tel = config.TELESCOPES[telescope]
    outdir = output_dir(year, telescope, "master_bias")
    outfile = outdir / f"master_bias_{year}_{telescope}_{tel['binning']}.fits"
    qc_file = outdir / "bias_qc.json"

    if outfile.exists() and not force:
        log.info("Master bias exists: %s (use --force to rebuild)", outfile)
        master = read_ccd(str(outfile))
        # Load existing QC
        import json
        with open(qc_file) as f:
            qc = json.load(f)
        return master, qc

    # QC
    qc = check_bias_frames(bias_entries, telescope)
    if qc["n_accepted"] == 0:
        log.error("No bias frames survived QC for %s %s", year, telescope)
        write_json(qc, qc_file)
        return None, qc

    # Load accepted frames
    log.info("Combining %d bias frames for %s...", qc["n_accepted"], telescope)
    paths = [e["path"] for e in bias_entries
             if e["filename"] in {a["filename"] for a in qc["accepted"]}]
    ccds = [read_ccd(p) for p in paths]

    # Combine
    master = ccdproc.combine(
        ccds, method="average",
        sigma_clip=True,
        sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
        sigma_clip_func=np.ma.median, sigma_clip_devfunc=np.ma.std,
    )

    # Add headers
    master.header["NCOMBINE"] = (qc["n_accepted"], "Number of frames combined")
    master.header["COMBMETH"] = ("average", "Combination method")
    master.header["BIASQUAL"] = (qc["bias_qual"], "Bias quality flag")
    master.header["PIPELINE"] = (config.PIPELINE_NAME, "Pipeline version")
    master.header["PROCDATE"] = (datetime.now(timezone.utc).isoformat(),
                                  "Processing date (UTC)")
    if qc["bias_qual"] == "TEMP_MISMATCH":
        master.header["COMMENT"] = (
            "WARNING: Bias CCD temperature does not match science temperature. "
            "Readout pattern is valid but absolute bias level may differ."
        )

    write_ccd(master, outfile)
    write_json(qc, qc_file)
    log.info("Master bias saved: %s (%d frames, qual=%s)",
             outfile.name, qc["n_accepted"], qc["bias_qual"])
    return master, qc


def make_master_flats(flat_entries_by_filter: dict, master_bias: CCDData,
                      year: str, telescope: str,
                      force: bool = False) -> tuple[dict, dict]:
    """Build master flats per filter.

    Returns (dict of filter→CCDData, flat_qc_report).
    """
    tel = config.TELESCOPES[telescope]
    outdir = output_dir(year, telescope, "master_flat")
    qc_file = outdir / "flat_qc.json"

    # QC all flats
    qc = check_flat_frames(flat_entries_by_filter, telescope,
                           master_bias_data=master_bias.data if master_bias else None)

    master_flats = {}
    for filt, filt_qc in qc.items():
        outfile = outdir / f"master_flat_{year}_{telescope}_{tel['binning']}_{filt}.fits"

        if outfile.exists() and not force:
            log.info("Master flat [%s] exists: %s", filt, outfile)
            master_flats[filt] = read_ccd(str(outfile))
            continue

        if filt_qc["n_accepted"] == 0:
            log.warning("No flat frames survived QC for %s %s [%s]",
                        year, telescope, filt)
            continue

        # Load and bias-subtract accepted frames
        accepted_names = {a["filename"] for a in filt_qc["accepted"]}
        entries = [e for e in flat_entries_by_filter[filt]
                   if e["filename"] in accepted_names]
        log.info("Combining %d flat frames for %s [%s]...",
                 len(entries), telescope, filt)

        ccds = []
        for e in entries:
            ccd = read_ccd(e["path"])
            if master_bias is not None:
                ccd = ccdproc.subtract_bias(ccd, master_bias)
            ccds.append(ccd)

        # Combine with median, scaling each to median=1
        master = ccdproc.combine(
            ccds, method="median",
            scale=lambda arr: 1.0 / np.ma.median(arr),
            sigma_clip=True,
            sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
        )

        # Ensure median is exactly 1.0
        master.data = master.data / np.median(master.data)

        # Headers
        master.header["NCOMBINE"] = (len(entries), "Number of frames combined")
        master.header["COMBMETH"] = ("median", "Combination method")
        master.header["FILTER"] = (filt, "Filter name")
        master.header["PIPELINE"] = (config.PIPELINE_NAME, "Pipeline version")
        master.header["PROCDATE"] = (datetime.now(timezone.utc).isoformat(),
                                      "Processing date (UTC)")

        write_ccd(master, outfile)
        master_flats[filt] = master
        log.info("Master flat [%s] saved: %s (%d frames)",
                 filt, outfile.name, len(entries))

    write_json(qc, qc_file)
    return master_flats, qc


def reduce_science_frame(entry: dict, master_bias: CCDData,
                         master_flats: dict, year: str, telescope: str,
                         force: bool = False) -> Path | None:
    """Bias-subtract and flat-correct a single science frame.

    Returns output path, or None if skipped.
    """
    tel = config.TELESCOPES[telescope]
    # Extract night directory from source path to avoid filename collisions
    # Source: .../T120/YYYYMMDD/file.fits or .../OHP_T120_YYYYMMDD/file.fits
    src = Path(entry["path"])
    night = src.parent.name  # e.g. "20250318" or "Calibration" subfolder
    if night == "Calibration":
        night = src.parent.parent.name
    outdir = output_dir(year, telescope, f"reduced/{night}")
    safe_name = entry["filename"].replace(" ", "_")
    outfile = outdir / safe_name

    if outfile.exists() and not force:
        return outfile

    # Quick QC: skip junk frames
    try:
        ccd = read_ccd(entry["path"])
    except Exception as exc:
        log.error("Cannot read %s: %s", entry["filename"], exc)
        return None

    hdr = ccd.header
    exptime = hdr.get("EXPTIME", 0)
    obj = hdr.get("OBJECT", "")

    # Skip bias/test frames that leaked into science index
    if exptime == 0:
        log.debug("Skip %s: EXPTIME=0", entry["filename"])
        return None
    obj_lower = obj.lower()
    if any(tok in obj_lower for tok in ("bias", "dark", "flat", "test")):
        log.warning("Skip %s: OBJECT='%s' is cal/test, not science",
                    entry["filename"], obj)
        return None

    # Check dimensions match telescope expectation
    if ccd.data.shape != tel["expected_shape"]:
        log.warning("Skip %s: shape %s != expected %s",
                    entry["filename"], ccd.data.shape, tel["expected_shape"])
        return None

    # Bias subtraction
    if master_bias is not None:
        ccd = ccdproc.subtract_bias(ccd, master_bias)
        ccd.header["BIASCORR"] = (True, "Bias subtracted")
    else:
        ccd.header["BIASCORR"] = (False, "No master bias available")

    # Flat correction
    filt = entry.get("filter", hdr.get("FILTER", ""))
    if filt in master_flats:
        ccd = ccdproc.flat_correct(ccd, master_flats[filt])
        ccd.header["FLATCORR"] = (filt, "Flat-corrected with this filter")
    else:
        ccd.header["FLATCORR"] = ("NONE", f"No master flat for filter {filt}")
        if filt:
            log.debug("No flat for %s [%s] — bias-only correction",
                      entry["filename"], filt)

    # Pipeline metadata
    ccd.header["PIPELINE"] = (config.PIPELINE_NAME, "Pipeline version")
    ccd.header["PROCDATE"] = (datetime.now(timezone.utc).isoformat(),
                               "Processing date (UTC)")

    write_ccd(ccd, outfile)
    return outfile


def run_calibration(archive: dict, year: str, telescope: str,
                    force: bool = False):
    """Full calibration pipeline for one year/telescope."""
    from .utils import get_cal_frames

    tel = config.TELESCOPES[telescope]
    cal = get_cal_frames(archive, year, telescope, tel["binning"])

    log.info("=" * 60)
    log.info("Calibrating %s %s", year, telescope)
    log.info("=" * 60)

    # ── Master bias ───────────────────────────────────────────────────
    bias_entries = cal.get("bias", [])
    if not bias_entries:
        log.warning("No bias frames found for %s %s", year, telescope)
    master_bias, bias_qc = make_master_bias(bias_entries, year, telescope,
                                             force=force)

    # ── Master flats ──────────────────────────────────────────────────
    flat_entries = cal.get("flat", {})
    if not flat_entries:
        log.warning("No flat frames found for %s %s", year, telescope)
    master_flats, flat_qc = make_master_flats(flat_entries, master_bias,
                                               year, telescope, force=force)

    log.info("Master flats available: %s", list(master_flats.keys()))

    # ── Science reduction ─────────────────────────────────────────────
    science = get_science_frames(archive, year, telescope)
    log.info("Reducing %d science frames for %s %s...", len(science), year, telescope)

    n_ok = 0
    n_skip = 0
    n_noflat = 0
    for i, entry in enumerate(science, 1):
        result = reduce_science_frame(entry, master_bias, master_flats,
                                       year, telescope, force=force)
        if result is None:
            n_skip += 1
        else:
            filt = entry.get("filter", "")
            if filt not in master_flats:
                n_noflat += 1
            n_ok += 1
            if i % 50 == 0 or i == len(science):
                log.info("  [%d/%d] reduced", i, len(science))

    # ── Summary ───────────────────────────────────────────────────────
    qc_dir = output_dir(year, telescope, "qc")
    summary = {
        "year": year,
        "telescope": telescope,
        "bias": {
            "n_input": bias_qc.get("n_input", 0) if bias_qc else 0,
            "n_accepted": bias_qc.get("n_accepted", 0) if bias_qc else 0,
            "quality": bias_qc.get("bias_qual", "N/A") if bias_qc else "N/A",
            "warnings": bias_qc.get("warnings", []) if bias_qc else [],
        },
        "flats": {filt: {"n_input": fq["n_input"], "n_accepted": fq["n_accepted"]}
                  for filt, fq in flat_qc.items()},
        "science": {
            "n_total": len(science),
            "n_reduced": n_ok,
            "n_skipped": n_skip,
            "n_bias_only": n_noflat,
        },
    }
    write_json(summary, qc_dir / "qc_report.json")
    log.info("Done %s %s: %d reduced, %d skipped, %d bias-only (no flat)",
             year, telescope, n_ok, n_skip, n_noflat)
    return summary
