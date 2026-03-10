#!/usr/bin/env python3
"""Batch SWarp stacking with shared WCS grids and flat-field weight maps.

Stacks all T120 targets (2018-2025) using existing SCAMP .head files.
For each target group, computes a common bounding box from ALL frames
across ALL filters, then forces SWarp to use the same center/size/pixel_scale
for every filter — guaranteeing pixel-aligned multi-band stacks.

Usage:
    python scripts/stack_all.py                    # stack everything
    python scripts/stack_all.py --target Coma      # one merged group
    python scripts/stack_all.py --dry-run          # just list what would be done
    python scripts/stack_all.py --force            # re-stack even if output exists
"""

import argparse
import json
import logging
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

# ── Pipeline imports ────────────────────────────────────────────────────

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config
from pipeline.stacking import (
    TARGET_ALIASES, SKIP_PATTERNS, ASTEROID_PATTERNS,
    build_target_inventory, normalise_target,
    _normalise_filter, _reject_spatial_outliers, _reject_outliers,
    _reject_spatial_outliers, _reject_outliers,
)

log = logging.getLogger(__name__)

# ── Constants ───────────────────────────────────────────────────────────

DATA_ROOT = config.DATA_ROOT
SWARP_CMD = config.SWARP_CMD
PIXEL_SCALE = 0.770       # arcsec/px — T120 output pixel scale
TARGET_ZP = 30.0          # target zero point for all stacks
YEARS = list(range(2018, 2026))
TELESCOPE = "T120"

# Spatially overlapping targets that share one output grid
MERGED_GROUPS = {
    "Coma": ["Abell1656", "Abell1656_border", "NGC4874"],
    "M105_group": ["M105", "NGC3384"],
    "M38_group": ["M38", "NGC1912"],
    "Markarian": ["M84", "M86"],
}

# Targets to skip: too few frames, single-band only, or no .head files
SKIP_TARGETS = {
    "J0254", "Gaia18aak", "IC10", "NGC628", "Abell1228",
    "NGC4636", "NGC4278", "NGC5322",
}


# ── Utilities ───────────────────────────────────────────────────────────

def find_master_flat(frame_path, year):
    """Find the master flat used for a given frame via its FLATCORR keyword."""
    try:
        hdr = fits.getheader(str(frame_path))
    except Exception:
        return None

    flatcorr = hdr.get("FLATCORR")
    if not flatcorr or flatcorr == "NONE":
        return None

    flat_path = (DATA_ROOT / str(year) / TELESCOPE / "master_flat" /
                 f"master_flat_{year}_{TELESCOPE}_2x2_{flatcorr}.fits")
    if flat_path.exists():
        return flat_path
    return None


def read_head_wcs(head_path, naxis1, naxis2):
    """Parse a SCAMP .head file into an astropy WCS object.

    The .head file contains FITS header cards but no NAXIS keywords,
    so we inject them before constructing the WCS.
    """
    text = Path(head_path).read_text()
    # Build a Header from the text
    hdr = fits.Header()
    hdr["NAXIS"] = 2
    hdr["NAXIS1"] = naxis1
    hdr["NAXIS2"] = naxis2

    for line in text.splitlines():
        # Skip HISTORY/COMMENT/blank/END
        if (line.startswith("HISTORY") or line.startswith("COMMENT") or
                line.strip() == "" or line.startswith("END")):
            continue
        # Parse standard 80-char FITS header card
        try:
            card = fits.Card.fromstring(line)
            hdr.append(card)
        except Exception:
            continue

    return WCS(hdr)


def compute_common_grid(frame_infos, max_offset_deg=1.0):
    """Compute a shared WCS grid (center, size) from all frames.

    First extracts CRVAL from each .head file, rejects positional outliers
    (frames with bad SCAMP solutions), then computes the bounding box from
    corners of surviving frames only.

    Parameters
    ----------
    frame_infos : list of dict
        Each dict has 'path', 'head_path', 'naxis1', 'naxis2'.
    max_offset_deg : float
        Maximum offset from median CRVAL to keep a frame.

    Returns
    -------
    (grid_dict, surviving_frame_infos) or (None, []) if failed.
    grid_dict has center_ra, center_dec, nx, ny.
    """
    # Expected |det(CD)| for T120: (0.770"/3600)^2 ≈ 4.6e-8
    expected_detcd = (PIXEL_SCALE / 3600.0) ** 2

    # Pass 1: read CRVAL from each .head, validate CD matrix
    crvals = []
    frame_wcs_list = []
    for fi in frame_infos:
        try:
            wcs = read_head_wcs(fi["head_path"], fi["naxis1"], fi["naxis2"])
        except Exception as e:
            log.warning("  Cannot read WCS from %s: %s", fi["head_path"], e)
            continue

        # Validate CD matrix: determinant should be close to expected
        try:
            cd = wcs.wcs.cd
            det = abs(cd[0, 0] * cd[1, 1] - cd[0, 1] * cd[1, 0])
            if det < expected_detcd * 0.1 or det > expected_detcd * 10:
                log.info("  Rejecting bad SCAMP WCS %s: |det(CD)|=%.2e "
                         "(expected ~%.2e)",
                         Path(fi["path"]).name, det, expected_detcd)
                continue
        except Exception:
            log.warning("  No CD matrix in %s", fi["head_path"])
            continue

        # CRVAL = reference point sky coords
        ra0 = wcs.wcs.crval[0]
        dec0 = wcs.wcs.crval[1]
        crvals.append((ra0, dec0))
        frame_wcs_list.append((fi, wcs))

    if len(crvals) < 2:
        return None, []

    ras0 = np.array([c[0] for c in crvals])
    decs0 = np.array([c[1] for c in crvals])

    # Handle RA wrapping for median
    if ras0.max() - ras0.min() > 180:
        ras0_wrap = np.where(ras0 > 180, ras0 - 360, ras0)
    else:
        ras0_wrap = ras0

    med_ra = float(np.median(ras0_wrap))
    med_dec = float(np.median(decs0))
    cos_dec = np.cos(np.radians(med_dec))

    # Reject .head WCS outliers
    good_frames = []
    all_ras = []
    all_decs = []
    for (fi, wcs), ra0, dec0 in zip(frame_wcs_list, ras0_wrap, decs0):
        offset = np.sqrt(((ra0 - med_ra) * cos_dec) ** 2 +
                         (dec0 - med_dec) ** 2)
        if offset > max_offset_deg:
            log.info("  Rejecting .head outlier %s: CRVAL=(%.3f, %.3f) "
                     "offset=%.2f°",
                     Path(fi["path"]).name, ra0, dec0, offset)
            continue

        good_frames.append(fi)
        # Compute corner coords
        nx, ny = fi["naxis1"], fi["naxis2"]
        corners_x = [0, nx - 1, 0, nx - 1]
        corners_y = [0, 0, ny - 1, ny - 1]
        try:
            ra_corners, dec_corners = wcs.pixel_to_world_values(
                corners_x, corners_y)
            all_ras.extend(ra_corners)
            all_decs.extend(dec_corners)
        except Exception as e:
            log.warning("  WCS corner transform failed for %s: %s",
                        fi["head_path"], e)

    if not all_ras:
        return None, []

    all_ras = np.array(all_ras)
    all_decs = np.array(all_decs)

    # Handle RA wrapping near 0/360
    ra_range = all_ras.max() - all_ras.min()
    if ra_range > 180:
        all_ras = np.where(all_ras > 180, all_ras - 360, all_ras)

    center_ra = float((all_ras.min() + all_ras.max()) / 2)
    center_dec = float((all_decs.min() + all_decs.max()) / 2)

    # Unwrap center back to [0, 360)
    if center_ra < 0:
        center_ra += 360

    # Image size in pixels
    cos_dec = np.cos(np.radians(center_dec))
    ra_span_deg = float(all_ras.max() - all_ras.min())
    dec_span_deg = float(all_decs.max() - all_decs.min())

    pixel_scale_deg = PIXEL_SCALE / 3600.0
    nx = int(np.ceil(ra_span_deg * cos_dec / pixel_scale_deg)) + 40
    ny = int(np.ceil(dec_span_deg / pixel_scale_deg)) + 40

    grid = {
        "center_ra": center_ra,
        "center_dec": center_dec,
        "nx": nx,
        "ny": ny,
    }
    return grid, good_frames


def gather_frame_info(frames):
    """For each frame, find its .head file and read NAXIS.

    Returns list of augmented frame dicts (only those with .head files).
    """
    result = []
    for f in frames:
        fpath = Path(f["path"])
        head_path = fpath.with_suffix(".head")
        if not head_path.exists():
            continue

        try:
            hdr = fits.getheader(str(fpath))
            naxis1 = hdr["NAXIS1"]
            naxis2 = hdr["NAXIS2"]
            exptime = float(hdr.get("EXPTIME", 1.0))
        except Exception:
            continue

        info = dict(f)
        info["head_path"] = str(head_path)
        info["naxis1"] = naxis1
        info["naxis2"] = naxis2
        info["exptime"] = exptime
        result.append(info)

    return result


def stack_filter_group(frames, common_grid, output_path, weight_path, tmpdir):
    """Run SWarp on a single filter group with shared grid.

    Parameters
    ----------
    frames : list of dict
        Augmented frame dicts with head_path, naxis1/2, flxscale.
    common_grid : dict
        From compute_common_grid: center_ra, center_dec, nx, ny.
    output_path : Path
        Output coadd FITS file.
    weight_path : Path
        Output weight FITS file.
    tmpdir : Path
        Temporary working directory.

    Returns
    -------
    dict with status info.
    """
    tmpdir = Path(tmpdir)
    swarp_inputs = []
    has_weights = []

    for f in frames:
        src = Path(f["path"])
        # Unique stem: year_night_originalname
        night = src.parent.name
        stem = f"{f['year']}_{night}_{src.stem}"

        # Symlink FITS to tmpdir
        dst = tmpdir / f"{stem}.fits"
        if not dst.exists():
            os.symlink(str(src.resolve()), str(dst))

        # Copy .head file (WCS override from SCAMP) and append FLXSCALE
        head_src = Path(f["head_path"])
        head_dst = tmpdir / f"{stem}.head"
        if not head_dst.exists():
            head_text = head_src.read_text()
            # Strip trailing END card and whitespace, then append FLXSCALE + END
            head_text = head_text.replace("END", "").rstrip()
            head_text += (
                f"\nFLXSCALE=   {f['flxscale']:.15E}"
                f" / Flux scale to ZP={TARGET_ZP}\n"
                f"END\n"
            )
            head_dst.write_text(head_text)

        # Find master flat for weight map
        flat_path = find_master_flat(src, f["year"])
        if flat_path:
            weight_dst = tmpdir / f"{stem}.weight.fits"
            if not weight_dst.exists():
                os.symlink(str(flat_path.resolve()), str(weight_dst))
            has_weights.append(True)
        else:
            has_weights.append(False)

        swarp_inputs.append(str(dst))

    if not swarp_inputs:
        return {"status": "NO_INPUTS"}

    # Build SWarp command
    grid = common_grid
    swarp_args = [
        SWARP_CMD,
        *swarp_inputs,
        "-IMAGEOUT_NAME", str(output_path),
        "-WEIGHTOUT_NAME", str(weight_path),
        "-COMBINE", "Y",
        "-COMBINE_TYPE", "MEDIAN",
        "-RESAMPLE", "Y",
        "-RESAMPLE_DIR", str(tmpdir),
        "-SUBTRACT_BACK", "Y",
        "-BACK_SIZE", "256",
        "-FSCALE_KEYWORD", "FLXSCALE",
        "-FSCALASTRO_TYPE", "NONE",
        "-CELESTIAL_TYPE", "NATIVE",
        "-PROJECTION_TYPE", "TAN",
        "-CENTER_TYPE", "MANUAL",
        "-CENTER", f"{grid['center_ra']:.8f},{grid['center_dec']:.8f}",
        "-IMAGE_SIZE", f"{grid['nx']},{grid['ny']}",
        "-PIXEL_SCALE", str(PIXEL_SCALE),
        "-PIXELSCALE_TYPE", "MANUAL",
        "-RESAMPLING_TYPE", "LANCZOS3",
        "-VMEM_DIR", str(tmpdir),
        "-DELETE_TMPFILES", "Y",
        "-VERBOSE_TYPE", "QUIET",
        "-WRITE_XML", "N",
    ]

    # Weight map handling
    if all(has_weights):
        swarp_args.extend(["-WEIGHT_TYPE", "MAP_WEIGHT",
                           "-WEIGHT_SUFFIX", ".weight.fits"])
    elif any(has_weights):
        # Mixed: build explicit weight image list
        weight_list = []
        for inp, has_w in zip(swarp_inputs, has_weights):
            if has_w:
                weight_list.append(inp.replace(".fits", ".weight.fits"))
            else:
                # Create uniform weight map
                uniform = Path(inp).with_name(
                    Path(inp).stem + ".weight.fits")
                if not uniform.exists():
                    # Read frame shape and make all-ones weight
                    hdr = fits.getheader(inp)
                    data = np.ones((hdr["NAXIS2"], hdr["NAXIS1"]),
                                   dtype=np.float32)
                    fits.writeto(str(uniform), data, overwrite=True)
                weight_list.append(str(uniform))
        swarp_args.extend(["-WEIGHT_TYPE", "MAP_WEIGHT",
                           "-WEIGHT_IMAGE", ",".join(weight_list)])
    else:
        swarp_args.extend(["-WEIGHT_TYPE", "NONE"])

    log.info("    SWarp: %d frames → %s", len(swarp_inputs), output_path.name)
    try:
        proc = subprocess.run(swarp_args, capture_output=True, text=True,
                              timeout=600)
        if proc.returncode != 0:
            log.error("    SWarp failed: %s",
                      proc.stderr[-500:] if proc.stderr else "")
            return {"status": "SWARP_FAILED", "n_input": len(swarp_inputs)}
    except subprocess.TimeoutExpired:
        log.error("    SWarp timed out")
        return {"status": "TIMEOUT", "n_input": len(swarp_inputs)}

    return {"status": "OK", "n_input": len(swarp_inputs)}


def stack_target_group(group_name, targets, inventory, output_base, force=False,
                       dry_run=False):
    """Stack all filters for a target group with a shared WCS grid.

    Parameters
    ----------
    group_name : str
        Output name for the group (e.g. "Coma", "M67").
    targets : list of str
        Canonical target names to merge.
    inventory : dict
        Full inventory from build_target_inventory.
    output_base : Path
        Base output directory.
    force : bool
        Re-stack even if output exists.
    dry_run : bool
        Just report what would be done.

    Returns
    -------
    dict of filter -> result dict.
    """
    # Collect all frames for constituent targets
    all_frames = []
    for t in targets:
        if t in inventory:
            # T120 only
            all_frames.extend(
                f for f in inventory[t] if f["telescope"] == TELESCOPE
            )

    if not all_frames:
        log.info("  %s: no T120 frames found", group_name)
        return {}

    # Filter to frames with .head files and read NAXIS
    frame_infos = gather_frame_info(all_frames)
    if not frame_infos:
        log.info("  %s: no frames with .head files", group_name)
        return {}

    # Reject spatial outliers (0.5 deg for merged groups, tighter for single)
    max_offset = 0.5 if len(targets) > 1 else 0.5
    frame_infos = _reject_spatial_outliers(frame_infos, max_offset_deg=max_offset)
    if len(frame_infos) < 2:
        log.info("  %s: too few frames after spatial filtering (%d)",
                 group_name, len(frame_infos))
        return {}

    # Compute common grid from ALL frames (all filters), rejecting .head outliers
    grid, frame_infos = compute_common_grid(frame_infos)
    if grid is None:
        log.info("  %s: could not compute grid", group_name)
        return {}
    if len(frame_infos) < 2:
        log.info("  %s: too few frames after .head outlier rejection", group_name)
        return {}

    # Group by normalised filter
    by_filter = defaultdict(list)
    for f in frame_infos:
        by_filter[f["filter"]].append(f)

    filters_str = ", ".join(f"{k}({len(v)})" for k, v in
                            sorted(by_filter.items()))
    log.info("  %s: %d frames, grid %dx%d @ (%.4f, %.4f), filters: %s",
             group_name, len(frame_infos), grid["nx"], grid["ny"],
             grid["center_ra"], grid["center_dec"], filters_str)

    if dry_run:
        return {filt: {"status": "DRY_RUN", "n_input": len(frs)}
                for filt, frs in by_filter.items()}

    out_dir = output_base / group_name / TELESCOPE
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    for filt in sorted(by_filter):
        filt_frames = by_filter[filt]

        if len(filt_frames) < 2:
            log.info("    %s: skipping (only %d frame)", filt, len(filt_frames))
            results[filt] = {"status": "TOO_FEW", "n_input": len(filt_frames)}
            continue

        output_file = out_dir / f"{group_name}_{filt}.fits"
        weight_file = out_dir / f"{group_name}_{filt}.weight.fits"

        if not force and output_file.exists():
            log.info("    %s: exists, skipping (use --force to redo)", filt)
            results[filt] = {"status": "EXISTS", "n_input": len(filt_frames)}
            continue

        # Reject ZP outliers
        filt_frames = _reject_outliers(filt_frames, sigma=2.0)
        if len(filt_frames) < 2:
            log.info("    %s: too few frames after ZP rejection", filt)
            results[filt] = {"status": "ZP_REJECTED", "n_input": 0}
            continue

        # Compute flux scaling: scale each frame to TARGET_ZP=30
        # Pipeline ZP is for counts/sec: mag = -2.5*log10(flux/exptime) + ZP
        # So FLXSCALE = 10^(0.4*(TARGET_ZP - ZP)) / exptime
        has_zps = any(f["photzp"] is not None for f in filt_frames)
        if has_zps:
            for f in filt_frames:
                if f["photzp"] is not None:
                    f["flxscale"] = (10 ** (0.4 * (TARGET_ZP - f["photzp"]))
                                     / f["exptime"])
                else:
                    f["flxscale"] = 1.0
        else:
            for f in filt_frames:
                f["flxscale"] = 1.0
            log.warning("    %s: no ZPs — stacking without flux scaling", filt)

        with tempfile.TemporaryDirectory() as td:
            res = stack_filter_group(filt_frames, grid,
                                     output_file, weight_file, td)

        # Update output header with provenance
        if res["status"] == "OK" and output_file.exists():
            with fits.open(str(output_file), mode="update") as hdul:
                h = hdul[0].header
                h["NCOMBINE"] = (res["n_input"], "Number of input frames")
                h["FILTER"] = (filt, "Filter band")
                h["OBJECT"] = (group_name, "Target group name")
                if has_zps:
                    h["PHOTZP"] = (TARGET_ZP,
                                   "Zero point (all frames scaled to this)")
                h["COMBTYPE"] = ("MEDIAN", "Combine method")
                h["BACKSUB"] = (True, "Background subtracted before combine")
                h["BACKSIZE"] = (256, "Background mesh size (pixels)")
                h["PIXSCALE"] = (PIXEL_SCALE, "Output pixel scale (arcsec/px)")
                # Record input frames
                for i, f in enumerate(filt_frames[:999]):
                    h[f"INP{i+1:04d}"] = (Path(f["path"]).name,
                                           f"Input frame {i+1}")
                hdul.flush()

            res["zp_ref"] = TARGET_ZP if has_zps else None
            res["output"] = str(output_file)
            res["weight"] = str(weight_file)
            log.info("    %s: OK (%d frames, ZP=%s)",
                     filt, res["n_input"],
                     f"{TARGET_ZP:.1f}" if has_zps else "none")

        results[filt] = res

    return results


# ── Main ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Batch SWarp stacking with shared WCS grids")
    parser.add_argument("--target", "-t",
                        help="Stack only this target/group (canonical name)")
    parser.add_argument("--dry-run", "-n", action="store_true",
                        help="List what would be done without running SWarp")
    parser.add_argument("--force", "-f", action="store_true",
                        help="Re-stack even if output exists")
    parser.add_argument("--years", type=int, nargs="+", default=YEARS,
                        help="Years to include (default: all)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%H:%M:%S",
    )

    log.info("Building target inventory for T120, years %s...", args.years)
    inventory = build_target_inventory(years=args.years, telescopes=[TELESCOPE])
    log.info("Found %d canonical targets", len(inventory))

    # Build group mapping: group_name -> list of canonical targets
    # Merged groups first, then individual targets not in any group
    merged_members = set()
    for members in MERGED_GROUPS.values():
        merged_members.update(members)

    groups = {}
    for gname, members in MERGED_GROUPS.items():
        # Only include if at least one member has frames
        if any(m in inventory for m in members):
            groups[gname] = members

    for target in sorted(inventory):
        if target in merged_members:
            continue
        if target in SKIP_TARGETS:
            continue
        groups[target] = [target]

    # Filter to requested target
    if args.target:
        t = args.target
        if t in groups:
            groups = {t: groups[t]}
        else:
            # Case-insensitive match
            matches = {k: v for k, v in groups.items()
                       if k.lower() == t.lower()}
            if matches:
                groups = matches
            else:
                # Check if it's a member of a merged group
                for gname, members in MERGED_GROUPS.items():
                    if t in members or t.lower() in [m.lower() for m in members]:
                        groups = {gname: members}
                        break
                else:
                    log.error("Target '%s' not found. Available groups:", t)
                    for g in sorted(groups):
                        log.info("  %s → %s", g, groups[g])
                    return

    stack_dir = DATA_ROOT / "stacks"
    all_results = {}

    log.info("Processing %d target groups%s",
             len(groups), " (dry run)" if args.dry_run else "")

    for group_name in sorted(groups):
        targets = groups[group_name]
        log.info("─── %s %s", group_name,
                 f"({', '.join(targets)})" if len(targets) > 1 else "")

        results = stack_target_group(
            group_name, targets, inventory, stack_dir,
            force=args.force, dry_run=args.dry_run,
        )
        if results:
            all_results[group_name] = results

    # Write JSON report
    if all_results and not args.dry_run:
        report_path = stack_dir / "stack_all_report.json"
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path, "w") as f:
            json.dump(all_results, f, indent=2, default=str)
        log.info("Wrote %s", report_path)

    # Print summary
    print(f"\n{'Group':25s}  {'Filter':8s}  {'N':>4s}  {'ZP':>7s}  {'Status':10s}")
    print("─" * 62)
    for group in sorted(all_results):
        for filt, res in sorted(all_results[group].items()):
            zp_str = f"{res['zp_ref']:.2f}" if res.get("zp_ref") else "—"
            n = res.get("n_input", 0)
            print(f"{group:25s}  {filt:8s}  {n:4d}  {zp_str:>7s}  "
                  f"{res['status']:10s}")


if __name__ == "__main__":
    main()
