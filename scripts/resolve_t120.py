#!/usr/bin/env python3
"""Re-solve all T120 2025 frames with astrometry.net solve-field, then re-run photometry.

Usage:
    python resolve_t120.py                # solve astrometry + photometry
    python resolve_t120.py --astrom-only  # solve astrometry only
    python resolve_t120.py --phot-only    # photometry only (assumes WCS already solved)
"""

import json
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import config

YEAR = "2025"
TELESCOPE = "T120"

logging.basicConfig(level=logging.WARNING, format="%(message)s")
log = logging.getLogger("pipeline")
log.setLevel(logging.INFO)


def run_astrometry_solvefield():
    """Solve all T120 reduced frames with solve-field + Gaia validation."""
    from pipeline.astrometry import (
        solve_frame_solvefield, validate_with_gaia,
        get_initial_wcs, resolve_object_coords,
        _parse_sexagesimal, _write_wcs_to_fits,
    )
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    import numpy as np
    import shutil

    if shutil.which("solve-field") is None:
        print("ERROR: solve-field not found on PATH")
        return False

    base = config.DATA_ROOT / YEAR / TELESCOPE
    reduced_dir = base / "reduced"
    gaia_cache_dir = base / "astrometry" / "gaia_cache"
    gaia_cache_dir.mkdir(parents=True, exist_ok=True)

    # Load simbad cache for name resolution
    simbad_cache = {}
    if config.SIMBAD_CACHE.exists():
        with open(config.SIMBAD_CACHE) as f:
            simbad_cache = json.load(f)

    fits_files = sorted(reduced_dir.rglob("*.fits")) + sorted(reduced_dir.rglob("*.fit"))
    if not fits_files:
        print("No reduced FITS files found")
        return False

    print(f"\n{'='*70}")
    print(f"  T120 Astrometry: {len(fits_files)} frames via solve-field")
    print(f"{'='*70}\n")

    # Pre-resolve object coordinates (one lookup per unique target)
    object_coords = {}
    objects_per_file = {}
    for fpath in fits_files:
        hdr = fits.getheader(str(fpath))
        obj = hdr.get("OBJECT", "unknown")
        objects_per_file[str(fpath)] = obj

        if obj not in object_coords:
            # Try header coords first
            wcs_init = get_initial_wcs(hdr, TELESCOPE)
            if wcs_init is not None:
                object_coords[obj] = SkyCoord(
                    ra=wcs_init.wcs.crval[0], dec=wcs_init.wcs.crval[1], unit="deg")
            else:
                coord = resolve_object_coords(obj, simbad_cache)
                if coord is not None:
                    object_coords[obj] = coord

    resolved = len([v for v in object_coords.values() if v is not None])
    print(f"  Resolved {resolved}/{len(object_coords)} unique targets\n")

    # Solve all frames
    results = []
    t_start = time.time()
    n_ok = 0
    n_fail = 0

    for i, fpath in enumerate(fits_files, 1):
        obj = objects_per_file[str(fpath)]
        center = object_coords.get(obj)

        t0 = time.time()
        result = solve_frame_solvefield(fpath, TELESCOPE, center=center)
        dt = time.time() - t0

        if result["status"] == "OK" and result.get("wcs") is not None:
            # Validate with Gaia
            try:
                nm, rms, ng = validate_with_gaia(
                    fpath, result["wcs"], TELESCOPE, gaia_cache_dir)
                result["n_matched"] = nm
                result["n_gaia"] = ng
                if rms is not None:
                    result["rms_arcsec"] = rms
                    if rms > config.ASTROM_MAX_RMS_ARCSEC:
                        result["status"] = "FAILED"
                    elif rms > config.ASTROM_WARN_RMS_ARCSEC:
                        result["status"] = "HIGH_RMS"
            except Exception as exc:
                log.warning("Gaia validation failed for %s: %s", fpath.name, exc)

            # Write WCS to FITS header
            _write_wcs_to_fits(fpath, result)
            n_ok += 1
        else:
            n_fail += 1

        # Store result for report
        res = {
            "filename": fpath.name,
            "night": fpath.parent.name,
            "object": obj,
            "filter": fits.getheader(str(fpath)).get("FILTER", ""),
            "status": result["status"],
            "n_matched": result.get("n_matched", 0),
            "n_gaia": result.get("n_gaia", 0),
            "rms_arcsec": result.get("rms_arcsec"),
            "crval1": result.get("crval1"),
            "crval2": result.get("crval2"),
            "rotation_deg": result.get("rotation_deg"),
            "pixel_scale": result.get("pixel_scale"),
            "fail_reason": result.get("fail_reason"),
        }
        results.append(res)

        # Progress every 10 frames, or on failure
        rms_str = f'{result.get("rms_arcsec", 0):.3f}"' if result.get("rms_arcsec") else "—"
        nm = result.get("n_matched", 0)
        if i % 10 == 0 or result["status"] != "OK":
            print(f"  [{i:3d}/{len(fits_files)}] {fpath.name:35s} "
                  f"{result['status']:8s}  {nm:3d} Gaia  RMS={rms_str:>7s}  "
                  f"{dt:.1f}s  [{n_ok} OK, {n_fail} fail]")

    elapsed = time.time() - t_start
    print(f"\n  Done in {elapsed:.0f}s — {n_ok} solved, {n_fail} failed\n")

    # Write report
    report_dir = base / "astrometry"
    report_dir.mkdir(parents=True, exist_ok=True)
    report_file = report_dir / "astrom_report.json"
    with open(report_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"  Report: {report_file}")

    # Summary by object
    from collections import Counter
    by_obj = {}
    for r in results:
        by_obj.setdefault(r["object"], []).append(r)

    print(f"\n  {'Object':25s}  {'N':>4s}  {'OK':>4s}  {'Fail':>4s}  {'Med RMS\"':>9s}")
    print(f"  {'─'*25}  {'─'*4}  {'─'*4}  {'─'*4}  {'─'*9}")
    import numpy as np
    for obj in sorted(by_obj):
        frames = by_obj[obj]
        n = len(frames)
        ok = sum(1 for r in frames if r["status"] in ("OK", "HIGH_RMS"))
        fail = sum(1 for r in frames if r["status"] == "FAILED")
        rms_vals = [r["rms_arcsec"] for r in frames if r.get("rms_arcsec")]
        med_rms = f"{np.median(rms_vals):.3f}" if rms_vals else "—"
        print(f"  {obj:25s}  {n:4d}  {ok:4d}  {fail:4d}  {med_rms:>9s}")

    return True


def run_photometry():
    """Re-run photometric calibration for T120."""
    from pipeline.photometry import run_photometry
    print(f"\n{'='*70}")
    print(f"  T120 Photometry: re-running with new WCS solutions")
    print(f"{'='*70}\n")
    run_photometry(YEAR, TELESCOPE, force=True)


if __name__ == "__main__":
    astrom = True
    phot = True

    if "--astrom-only" in sys.argv:
        phot = False
    elif "--phot-only" in sys.argv:
        astrom = False

    if astrom:
        ok = run_astrometry_solvefield()
        if not ok and phot:
            print("Astrometry failed — skipping photometry")
            sys.exit(1)

    if phot:
        run_photometry()

    print("\nDone.")
