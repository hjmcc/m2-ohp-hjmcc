#!/usr/bin/env python3
"""Run SCAMP astrometric refinement on all solved T120 frames.

For each target/filter group:
  1. SExtractor → LDAC catalogues
  2. SCAMP → refined WCS .head files (against Gaia DR3)

The .head files are saved alongside the reduced FITS frames and can be
read directly by SWarp for stacking.

Usage:
    python run_scamp.py                      # all years, T120
    python run_scamp.py --years 2025         # single year
    python run_scamp.py --target M67         # single target
    python run_scamp.py --dry-run            # show inventory only
"""

import argparse
import json
import logging
import shutil
import subprocess
import tempfile
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import config
from pipeline.stacking import (
    build_target_inventory, normalise_target, _normalise_filter,
    TARGET_ALIASES, SKIP_PATTERNS, ASTEROID_PATTERNS,
)

logging.basicConfig(level=logging.WARNING, format="%(message)s")
log = logging.getLogger("pipeline")
log.setLevel(logging.INFO)

SCAMP_CMD = "scamp"
TELESCOPE = "T120"


def _ensure_single_hdu(filepath, tmpdir):
    """If a FITS file has multiple extensions (MASK, UNCERT from ccdproc),
    write a single-HDU copy so SExtractor produces a clean 1-extension LDAC.

    Returns path to the single-HDU file (may be the original if already single).
    """
    with fits.open(str(filepath)) as hdul:
        if len(hdul) == 1:
            return str(filepath)
        # Write primary HDU only
        out = Path(tmpdir) / Path(filepath).name
        fits.writeto(str(out), hdul[0].data, hdul[0].header, overwrite=True)
        return str(out)


def run_sextractor_ldac(filepath, ldac_path, telescope="T120"):
    """Run SExtractor on a frame, producing LDAC catalogue for SCAMP."""
    tel = config.TELESCOPES[telescope]
    pixel_scale = tel["pixel_scale"]
    saturation = tel["saturation"]

    tmpdir = Path(ldac_path).parent

    # Ensure single-HDU input (multi-extension frames confuse SCAMP)
    input_path = _ensure_single_hdu(filepath, tmpdir)

    param_path = tmpdir / "scamp.param"
    if not param_path.exists():
        param_path.write_text(
            "XWIN_IMAGE\nYWIN_IMAGE\n"
            "ERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\n"
            "FLUX_AUTO\nFLUXERR_AUTO\n"
            "MAG_AUTO\nMAGERR_AUTO\n"
            "FLAGS\nFLUX_RADIUS\nELONGATION\nSNR_WIN\n"
        )

    sex_path = tmpdir / f"{Path(filepath).stem}.sex"
    sex_path.write_text(
        f"CATALOG_NAME     {ldac_path}\n"
        f"CATALOG_TYPE     FITS_LDAC\n"
        f"PARAMETERS_NAME  {param_path}\n"
        f"DETECT_TYPE      CCD\n"
        f"DETECT_THRESH    5.0\n"
        f"ANALYSIS_THRESH  5.0\n"
        f"FILTER           N\n"
        f"DEBLEND_NTHRESH  32\n"
        f"DEBLEND_MINCONT  0.005\n"
        f"CLEAN            Y\n"
        f"SATUR_LEVEL      {saturation}\n"
        f"PIXEL_SCALE      {pixel_scale}\n"
        f"SEEING_FWHM      3.0\n"
        f"BACK_SIZE        64\n"
        f"BACK_FILTERSIZE  3\n"
        f"WEIGHT_TYPE      NONE\n"
        f"CHECKIMAGE_TYPE  NONE\n"
        f"VERBOSE_TYPE     QUIET\n"
    )

    proc = subprocess.run(
        [config.SEXTRACTOR_CMD, str(input_path), "-c", str(sex_path)],
        capture_output=True, text=True, timeout=120,
    )
    if proc.returncode != 0:
        return False
    return Path(ldac_path).exists()


def run_scamp_group(ldac_paths, head_mapping, workdir):
    """Run SCAMP on a group of LDAC catalogues.

    Parameters
    ----------
    ldac_paths : list of Path
        LDAC catalogue paths in the working directory.
    head_mapping : dict
        Maps LDAC stem → final .head destination path.
    workdir : Path
        SCAMP working directory.

    Returns
    -------
    dict with n_input, n_heads, internal_rms, external_rms, status
    """
    scamp_args = [
        SCAMP_CMD,
        *[str(p) for p in ldac_paths],
        "-ASTREF_CATALOG", "GAIA-DR3",
        "-DISTORT_DEGREES", "1",
        "-STABILITY_TYPE", "INSTRUMENT",
        "-CROSSID_RADIUS", "3.0",
        "-POSITION_MAXERR", "2.0",
        "-POSANGLE_MAXERR", "2.0",
        "-PIXSCALE_MAXERR", "1.05",
        "-SN_THRESHOLDS", "5.0,20.0",
        "-ASTREFMAG_LIMITS", "10.0,19.0",
        "-MATCH", "Y",
        "-MATCH_FLIPPED", "N",
        "-MOSAIC_TYPE", "UNCHANGED",
        "-SOLVE_PHOTOM", "N",
        "-WRITE_XML", "N",
        "-CHECKPLOT_TYPE", "NONE",
        "-VERBOSE_TYPE", "NORMAL",
    ]

    try:
        proc = subprocess.run(
            scamp_args, capture_output=True,
            timeout=600, cwd=str(workdir),
        )
    except subprocess.TimeoutExpired:
        return {"status": "TIMEOUT", "n_input": len(ldac_paths), "n_heads": 0}

    # Parse SCAMP output for stats (stderr may contain non-UTF-8 bytes)
    stderr_text = proc.stderr.decode("utf-8", errors="replace")
    ext_rms = None
    for line in stderr_text.splitlines():
        if "Astrometric stats (external)" in line:
            pass
        # Lines like: Group  1: 0.065"  0.107"  ...
        if "Group" in line and '"' in line:
            parts = line.split()
            for i, p in enumerate(parts):
                if p.endswith('"') and i >= 2:
                    try:
                        ext_rms = float(p.rstrip('"'))
                    except ValueError:
                        pass
                    break

    if proc.returncode != 0:
        return {
            "status": "FAILED",
            "n_input": len(ldac_paths),
            "n_heads": 0,
            "stderr": stderr_text[-300:] if stderr_text else "",
        }

    # Copy .head files to final destinations
    n_heads = 0
    for ldac in ldac_paths:
        head_src = ldac.with_suffix(".head")
        stem = ldac.stem
        if head_src.exists() and stem in head_mapping:
            dest = head_mapping[stem]
            shutil.copy2(str(head_src), str(dest))
            n_heads += 1

    return {
        "status": "OK",
        "n_input": len(ldac_paths),
        "n_heads": n_heads,
        "external_rms": ext_rms,
    }


def main():
    parser = argparse.ArgumentParser(description="Run SCAMP on solved T120 frames")
    parser.add_argument("--years", nargs="+", type=int,
                        default=[2018, 2019, 2020, 2021, 2023, 2024, 2025])
    parser.add_argument("--target", type=str, default=None,
                        help="Process only this target (canonical name)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show inventory without processing")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if .head files exist")
    args = parser.parse_args()

    if shutil.which(SCAMP_CMD) is None:
        print("ERROR: scamp not found on PATH")
        return
    if shutil.which(config.SEXTRACTOR_CMD) is None:
        print("ERROR: SExtractor not found on PATH")
        return

    # Build inventory of all T120 frames with good astrometry
    inventory = build_target_inventory(years=args.years, telescopes=[TELESCOPE])

    # Filter to frames with WCS and good astrometry
    filtered = {}
    for target, frames in inventory.items():
        good = [f for f in frames
                if f["has_wcs"] and f.get("astrstat") in ("OK", "HIGH_RMS")]
        if good:
            filtered[target] = good

    if args.target:
        matches = {k: v for k, v in filtered.items()
                   if k.lower() == args.target.lower()}
        if not matches:
            print(f"Target '{args.target}' not found. Available targets:")
            for t in sorted(filtered):
                print(f"  {t} ({len(filtered[t])} frames)")
            return
        filtered = matches

    # Summary
    total_frames = sum(len(v) for v in filtered.values())
    print(f"\n{'='*70}")
    print(f"  SCAMP astrometric refinement — T120")
    print(f"  {len(filtered)} targets, {total_frames} frames")
    print(f"{'='*70}\n")

    if args.dry_run:
        print(f"{'Target':25s}  {'Filter':8s}  {'Years':20s}  {'N':>4s}")
        print("-" * 65)
        for target in sorted(filtered):
            frames = filtered[target]
            by_filt = defaultdict(list)
            for f in frames:
                by_filt[f["filter"]].append(f)
            for filt in sorted(by_filt):
                grp = by_filt[filt]
                yrs = ",".join(str(y) for y in sorted(set(f["year"] for f in grp)))
                print(f"{target:25s}  {filt:8s}  {yrs:20s}  {len(grp):4d}")
        return

    # Process each target
    t_start = time.time()
    all_results = {}

    for target in sorted(filtered):
        frames = filtered[target]

        # Group by filter (don't mix filters in SCAMP)
        by_filt = defaultdict(list)
        for f in frames:
            by_filt[f["filter"]].append(f)

        for filt in sorted(by_filt):
            group = by_filt[filt]
            group_key = f"{target}/{filt}"

            # Check if .head files already exist
            if not args.force:
                n_existing = sum(
                    1 for f in group
                    if Path(f["path"]).with_suffix(".head").exists()
                )
                if n_existing == len(group):
                    all_results[group_key] = {
                        "status": "SKIPPED", "n_input": len(group),
                        "n_heads": n_existing,
                    }
                    continue

            years_str = ",".join(str(y) for y in sorted(set(f["year"] for f in group)))
            print(f"  {target:20s}  {filt:5s}  {len(group):3d} frames  "
                  f"({years_str})", end="", flush=True)

            if len(group) < 1:
                print("  — skipped (too few)")
                continue

            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)

                # Step 1: Run SExtractor on each frame → LDAC
                ldac_paths = []
                head_mapping = {}  # ldac_stem → final .head path

                for f in group:
                    src = Path(f["path"])
                    # Use unique stem: year_night_filename
                    stem = f"{f['year']}_{src.parent.name}_{src.stem}"
                    ldac = tmpdir / f"{stem}.ldac"

                    if run_sextractor_ldac(str(src), str(ldac), TELESCOPE):
                        ldac_paths.append(ldac)
                        # .head goes next to the original FITS file
                        head_mapping[stem] = src.with_suffix(".head")

                if len(ldac_paths) < 1:
                    print("  — SExtractor failed on all frames")
                    all_results[group_key] = {
                        "status": "SEX_FAILED", "n_input": len(group), "n_heads": 0,
                    }
                    continue

                # Step 2: Run SCAMP
                result = run_scamp_group(ldac_paths, head_mapping, tmpdir)
                all_results[group_key] = result

                rms_str = f"{result['external_rms']:.3f}\"" if result.get("external_rms") else "—"
                print(f"  → {result['n_heads']}/{len(ldac_paths)} .head  "
                      f"RMS={rms_str}  [{result['status']}]")

    elapsed = time.time() - t_start

    # Summary
    n_ok = sum(1 for r in all_results.values() if r["status"] == "OK")
    n_skip = sum(1 for r in all_results.values() if r["status"] == "SKIPPED")
    n_fail = sum(1 for r in all_results.values()
                 if r["status"] not in ("OK", "SKIPPED"))
    total_heads = sum(r.get("n_heads", 0) for r in all_results.values())

    print(f"\n{'='*70}")
    print(f"  Done in {elapsed:.0f}s")
    print(f"  {n_ok} groups processed, {n_skip} skipped, {n_fail} failed")
    print(f"  {total_heads} .head files written")
    print(f"{'='*70}")

    # Save report
    report_path = Path(__file__).parent / "scamp_report.json"
    with open(report_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"  Report: {report_path}")


if __name__ == "__main__":
    main()
