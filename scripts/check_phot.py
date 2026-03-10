#!/usr/bin/env python3
"""Validate and correct photometric calibration of stacks against Pan-STARRS1.

For each stack with a PS1-compatible filter:
1. Measure FWHM the same way pipeline QC does (sep.flux_radius)
2. Extract aperture photometry the same way pipeline photometry does
   (sep.sum_circle with r = PHOT_APERTURE_MULT * fwhm)
3. Cross-match against PS1 DR2, compute zero point
4. Report the offset from the target ZP (30.0)

With --fix: apply a multiplicative correction so the stack's ZP becomes
exactly 30.0. Pixel data *= 10^(0.4 * offset), weight /= factor^2.

Usage:
    python scripts/check_phot.py                    # report only
    python scripts/check_phot.py --fix              # report + apply corrections
    python scripts/check_phot.py --target M67       # single target
    python scripts/check_phot.py --plot             # save diagnostic plots
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import sep
sep.set_sub_object_limit(4096)
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pipeline import config
from pipeline.photometry import query_ps1, compute_zeropoint

log = logging.getLogger(__name__)

TARGET_ZP = 30.0
STACK_DIR = config.DATA_ROOT / "stacks"
PS1_CACHE = config.DATA_ROOT / "stacks" / "ps1_cache"

# Filter → PS1 band mapping (normalised filter names used in stacks)
FILTER_TO_PS1 = {
    "g": "gmag",
    "r": "rmag",
    "R": "rmag",
    "i": "imag",
    "V": "gmag",
    "B": "gmag",
}
# Filters with no PS1 counterpart — skip
SKIP_FILTERS = {"Ha", "OIII", "SII", "U"}

# Minimum offset (mag) to bother correcting
MIN_OFFSET_TO_FIX = 0.01


def measure_fwhm(data, mask):
    """Measure FWHM exactly as pipeline/quality.py does: sep.flux_radius.

    Returns median FWHM in pixels, or 4.0 as fallback.
    """
    work = np.ascontiguousarray(data, dtype=np.float64)
    bkg = sep.Background(work, mask=mask)
    data_sub = work - bkg

    objects = sep.extract(data_sub, config.SEP_THRESH, err=bkg.globalrms,
                          mask=mask, minarea=config.SEP_MINAREA)
    if len(objects) == 0:
        return 4.0

    x, y = objects["x"], objects["y"]
    a, b = objects["a"], objects["b"]

    # Half-light radius (same as quality.py)
    r_half, flag_r = sep.flux_radius(data_sub, x, y, 6.0 * a, 0.5, subpix=5)
    fwhm_all = 2.0 * r_half

    # Elongation
    elong_all = np.where(b > 0, a / b, np.nan)

    # Same quality cuts as quality.py (except saturation — not applicable)
    good_mask = (
        (objects["flag"] == 0)
        & (flag_r == 0)
        & (fwhm_all >= config.FWHM_MIN_PX)
        & (fwhm_all <= config.FWHM_MAX_PX)
        & np.isfinite(fwhm_all)
        & np.isfinite(elong_all)
    )

    if good_mask.sum() == 0:
        return 4.0

    return float(np.median(fwhm_all[good_mask]))


def extract_photometry_stack(data, header, fwhm_px):
    """Extract aperture photometry exactly as pipeline/photometry.py does.

    Same as extract_photometry() but adapted for stacks:
    - Masks zero-coverage regions
    - No saturation cut (pixel values are flux-scaled)
    """
    # Mask zero regions (outside coverage)
    mask = (data == 0).astype(np.uint8)

    work = np.ascontiguousarray(data, dtype=np.float64)
    bkg = sep.Background(work, mask=mask)
    data_sub = work - bkg

    # Same detection parameters as pipeline photometry
    objects = sep.extract(data_sub, config.PHOT_SEP_THRESH,
                          err=bkg.globalrms, mask=mask,
                          minarea=config.SEP_MINAREA)
    if len(objects) == 0:
        return None

    # Same quality cuts as pipeline photometry (minus saturation)
    obj_mask = (
        (objects["flag"] == 0)
        & (objects["a"] > 0) & (objects["b"] > 0)
        & (objects["a"] / objects["b"] < 3.0)
    )
    obj = objects[obj_mask]
    if len(obj) < config.PHOT_ZP_MIN_STARS:
        return None

    x = obj["x"]
    y = obj["y"]

    # Same aperture as pipeline photometry: PHOT_APERTURE_MULT * fwhm_px
    r_aper = config.PHOT_APERTURE_MULT * fwhm_px

    flux, fluxerr, flag = sep.sum_circle(data_sub, x, y, r_aper,
                                          mask=mask,
                                          err=bkg.globalrms, subpix=5)

    # Keep positive flux, low flags
    good = (flux > 0) & (flag == 0)
    if good.sum() < config.PHOT_ZP_MIN_STARS:
        return None

    x, y = x[good], y[good]
    flux, fluxerr = flux[good], fluxerr[good]

    # Convert pixel → sky coords
    try:
        wcs = WCS(header)
        sky = wcs.pixel_to_world(x, y)
        ra = np.asarray(sky.ra.deg)
        dec = np.asarray(sky.dec.deg)
    except Exception:
        return None

    inst_mag = -2.5 * np.log10(flux)
    inst_mag_err = 2.5 / np.log(10) * fluxerr / flux

    return {
        "x": x, "y": y,
        "ra": ra, "dec": dec,
        "flux": flux, "fluxerr": fluxerr,
        "inst_mag": inst_mag, "inst_mag_err": inst_mag_err,
        "fwhm": fwhm_px,
        "r_aper": r_aper,
        "n_sources": len(x),
    }


def cross_match(phot_ra, phot_dec, ps1_stars, match_radius_arcsec=3.0):
    """Cross-match extracted sources against PS1 catalogue.

    Returns (idx_phot, idx_ps1) arrays of matched pairs.
    """
    if not ps1_stars:
        return np.array([], dtype=int), np.array([], dtype=int)

    cat_ra = np.array([s["ra"] for s in ps1_stars])
    cat_dec = np.array([s["dec"] for s in ps1_stars])

    cos_dec = np.cos(np.radians(np.median(phot_dec)))

    cat_xy = np.column_stack([cat_ra * cos_dec, cat_dec])
    phot_xy = np.column_stack([phot_ra * cos_dec, phot_dec])

    tree = cKDTree(cat_xy)
    radius_deg = match_radius_arcsec / 3600.0

    dist, idx = tree.query(phot_xy, k=1)

    matched = dist < radius_deg
    idx_phot = np.where(matched)[0]
    idx_ps1 = idx[matched]

    return idx_phot, idx_ps1


def check_one_stack(stack_path, fix=False, make_plot=False):
    """Check photometric calibration of one stack file.

    Returns dict with results, or None if not applicable.
    """
    stack_path = Path(stack_path)
    hdr = fits.getheader(str(stack_path))

    filt = hdr.get("FILTER", "")
    obj = hdr.get("OBJECT", "")
    photzp = hdr.get("PHOTZP")
    ncombine = hdr.get("NCOMBINE", 0)

    if filt in SKIP_FILTERS:
        return None

    ps1_band = FILTER_TO_PS1.get(filt)
    if ps1_band is None:
        log.info("  %s: filter '%s' has no PS1 mapping, skipping",
                 stack_path.name, filt)
        return None

    # Read data
    data = fits.getdata(str(stack_path)).astype(np.float64)
    mask = (data == 0).astype(np.uint8)

    # Measure FWHM the same way QC does
    fwhm_px = measure_fwhm(data, mask)

    # Extract photometry the same way the pipeline does
    phot = extract_photometry_stack(data, hdr, fwhm_px)
    if phot is None:
        log.warning("  %s: source extraction failed", stack_path.name)
        return {"file": stack_path.name, "filter": filt, "object": obj,
                "status": "EXTRACTION_FAILED"}

    # Query PS1
    wcs = WCS(hdr)
    ny, nx = data.shape
    center_sky = wcs.pixel_to_world(nx / 2, ny / 2)
    center_ra = float(center_sky.ra.deg)
    center_dec = float(center_sky.dec.deg)

    # Search radius: half the diagonal in arcmin + margin
    pixel_scale = abs(hdr.get("CD1_1", 2.14e-4)) * 3600  # arcsec/px
    diag_arcmin = np.sqrt(nx**2 + ny**2) * pixel_scale / 60.0 / 2 + 1.0

    ps1_stars = query_ps1(center_ra, center_dec, diag_arcmin,
                          ps1_band, PS1_CACHE)
    if not ps1_stars:
        log.warning("  %s: no PS1 stars returned", stack_path.name)
        return {"file": stack_path.name, "filter": filt, "object": obj,
                "status": "NO_PS1"}

    # Cross-match
    idx_phot, idx_ps1 = cross_match(phot["ra"], phot["dec"], ps1_stars)
    if len(idx_phot) < 5:
        log.warning("  %s: only %d matches", stack_path.name, len(idx_phot))
        return {"file": stack_path.name, "filter": filt, "object": obj,
                "status": "FEW_MATCHES", "n_matches": len(idx_phot)}

    inst_mags = phot["inst_mag"][idx_phot]
    cat_mags = np.array([ps1_stars[i]["mag"] for i in idx_ps1])

    # Compute measured ZP (same sigma-clipped median as pipeline)
    zp_meas, zp_err, n_used = compute_zeropoint(inst_mags, cat_mags)
    if zp_meas is None:
        log.warning("  %s: ZP fit failed", stack_path.name)
        return {"file": stack_path.name, "filter": filt, "object": obj,
                "status": "ZP_FIT_FAILED", "n_matches": len(idx_phot)}

    offset = zp_meas - TARGET_ZP
    residuals = cat_mags - (inst_mags + zp_meas)

    result = {
        "file": stack_path.name,
        "filter": filt,
        "ps1_band": ps1_band,
        "object": obj,
        "ncombine": ncombine,
        "header_zp": photzp,
        "measured_zp": round(zp_meas, 4),
        "zp_err": round(zp_err, 4),
        "offset": round(offset, 4),
        "n_matched": len(idx_phot),
        "n_used": n_used,
        "rms": round(float(np.std(residuals)), 4),
        "fwhm_px": round(phot["fwhm"], 2),
        "r_aper_px": round(phot["r_aper"], 1),
        "n_extracted": phot["n_sources"],
        "status": "OK",
    }

    # Save diagnostic plot
    if make_plot:
        _save_plot(stack_path, inst_mags, cat_mags, zp_meas, n_used,
                   idx_phot, phot, ps1_stars, idx_ps1, result)

    # Apply correction
    if fix and abs(offset) >= MIN_OFFSET_TO_FIX:
        factor = 10 ** (0.4 * offset)
        log.info("  %s: applying correction factor %.6f (offset=%.4f mag)",
                 stack_path.name, factor, offset)

        with fits.open(str(stack_path), mode="update") as hdul:
            hdul[0].data = (hdul[0].data * factor).astype(np.float32)
            hdul[0].header["PHOTZP"] = (TARGET_ZP,
                                         "Zero point (corrected to PS1)")
            hdul[0].header["ZPOFFSET"] = (round(offset, 4),
                                           "ZP correction applied (mag)")
            hdul[0].header["ZPFACTOR"] = (round(factor, 6),
                                           "Flux correction factor applied")
            hdul[0].header["ZPMEAS"] = (round(zp_meas, 4),
                                         "Measured ZP before correction")
            hdul[0].header["ZPNSTAR"] = (n_used,
                                          "N stars used for ZP measurement")
            hdul[0].header["ZPRMS"] = (round(float(np.std(residuals)), 4),
                                        "ZP residual RMS (mag)")
            hdul.flush()

        # Also correct weight map
        weight_path = stack_path.with_name(
            stack_path.name.replace(".fits", ".weight.fits"))
        if weight_path.exists():
            with fits.open(str(weight_path), mode="update") as whdul:
                whdul[0].data = (whdul[0].data / (factor ** 2)).astype(
                    np.float32)
                whdul.flush()

        result["fixed"] = True
        result["factor"] = round(factor, 6)
    elif fix and abs(offset) < MIN_OFFSET_TO_FIX:
        result["fixed"] = False
        result["reason"] = f"offset {offset:.4f} < {MIN_OFFSET_TO_FIX}"

    return result


def _save_plot(stack_path, inst_mags, cat_mags, zp_meas, n_used,
               idx_phot, phot, ps1_stars, idx_ps1, result):
    """Save a diagnostic plot: PS1 mag vs calibrated mag + residuals."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    cal_mags = inst_mags + zp_meas
    residuals = cat_mags - cal_mags

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left: 1:1 plot
    ax1.scatter(cat_mags, cal_mags, s=8, alpha=0.5, c="steelblue")
    xlim = ax1.get_xlim()
    ax1.plot(xlim, xlim, "k--", lw=0.8, alpha=0.5)
    ax1.set_xlabel(f"PS1 {result['ps1_band']}")
    ax1.set_ylabel(f"Calibrated mag (ZP={zp_meas:.2f})")
    ax1.set_title(f"{result['object']} {result['filter']}: "
                  f"ZP={zp_meas:.3f}±{result['zp_err']:.3f}, "
                  f"offset={result['offset']:+.3f}")
    ax1.set_xlim(xlim)
    ax1.invert_xaxis()
    ax1.invert_yaxis()

    # Right: residuals vs magnitude
    ax2.scatter(cat_mags, residuals, s=8, alpha=0.5, c="steelblue")
    ax2.axhline(0, color="k", ls="--", lw=0.8, alpha=0.5)
    ax2.axhline(result["rms"], color="r", ls=":", lw=0.8, alpha=0.5)
    ax2.axhline(-result["rms"], color="r", ls=":", lw=0.8, alpha=0.5)
    ax2.set_xlabel(f"PS1 {result['ps1_band']}")
    ax2.set_ylabel("Residual (PS1 - calibrated)")
    ax2.set_title(f"N={result['n_used']}, RMS={result['rms']:.3f} mag")
    ax2.set_ylim(-0.5, 0.5)
    ax2.invert_xaxis()

    fig.tight_layout()
    plot_dir = stack_path.parent / "plots"
    plot_dir.mkdir(exist_ok=True)
    plot_path = plot_dir / f"{stack_path.stem}_phot.png"
    fig.savefig(str(plot_path), dpi=120)
    plt.close(fig)
    result["plot"] = str(plot_path)


def main():
    parser = argparse.ArgumentParser(
        description="Check and correct stack photometric calibration vs PS1")
    parser.add_argument("--fix", action="store_true",
                        help="Apply flux corrections to match ZP=30")
    parser.add_argument("--plot", action="store_true",
                        help="Save diagnostic plots")
    parser.add_argument("--target", "-t",
                        help="Check only this target group")
    parser.add_argument("--filter", "-f", dest="filt",
                        help="Check only this filter")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%H:%M:%S",
    )

    # Find all stack files
    stack_files = sorted(STACK_DIR.rglob("*.fits"))
    stack_files = [f for f in stack_files
                   if not f.name.endswith(".weight.fits")
                   and f.parent.name == "T120"]

    if args.target:
        stack_files = [f for f in stack_files
                       if f.parent.parent.name.lower() == args.target.lower()]
    if args.filt:
        stack_files = [f for f in stack_files
                       if f"_{args.filt}." in f.name]

    if not stack_files:
        log.error("No stack files found in %s", STACK_DIR)
        return

    log.info("Checking %d stacks%s", len(stack_files),
             " (will fix)" if args.fix else "")

    results = []
    for sf in stack_files:
        log.info("─── %s", sf.name)
        res = check_one_stack(sf, fix=args.fix, make_plot=args.plot)
        if res is not None:
            results.append(res)

    # Write JSON report
    report_path = STACK_DIR / "phot_check_report.json"
    with open(report_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    log.info("Wrote %s", report_path)

    # Print summary table
    ok_results = [r for r in results if r.get("status") == "OK"]
    fail_results = [r for r in results if r.get("status") != "OK"]

    print(f"\n{'File':35s}  {'Filt':5s}  {'PS1':5s}  {'N':>4s}  "
          f"{'ZP_meas':>8s}  {'Offset':>8s}  {'RMS':>6s}  {'FWHM':>5s}  {'Fixed':6s}")
    print("─" * 93)

    for r in sorted(ok_results, key=lambda x: abs(x.get("offset", 0)),
                    reverse=True):
        fixed_str = ""
        if r.get("fixed"):
            fixed_str = f"×{r['factor']:.4f}"
        elif r.get("fixed") is False:
            fixed_str = "skip"

        print(f"{r['file']:35s}  {r['filter']:5s}  {r['ps1_band']:5s}  "
              f"{r['n_used']:4d}  {r['measured_zp']:8.4f}  "
              f"{r['offset']:+8.4f}  {r['rms']:6.4f}  {r['fwhm_px']:5.2f}  {fixed_str:6s}")

    if fail_results:
        print(f"\nFailed ({len(fail_results)}):")
        for r in fail_results:
            print(f"  {r['file']:35s}  {r.get('status', '?')}")

    # Summary stats
    if ok_results:
        offsets = [r["offset"] for r in ok_results]
        print(f"\nSummary: {len(ok_results)} stacks checked, "
              f"{len(fail_results)} failed")
        print(f"  Offset: median={np.median(offsets):+.4f}, "
              f"mean={np.mean(offsets):+.4f}, "
              f"std={np.std(offsets):.4f}, "
              f"range=[{min(offsets):+.4f}, {max(offsets):+.4f}]")
        rms_vals = [r["rms"] for r in ok_results]
        print(f"  RMS:    median={np.median(rms_vals):.4f}, "
              f"range=[{min(rms_vals):.4f}, {max(rms_vals):.4f}]")
        if args.fix:
            n_fixed = sum(1 for r in ok_results if r.get("fixed"))
            n_skipped = sum(1 for r in ok_results if r.get("fixed") is False)
            print(f"  Fixed: {n_fixed}, skipped (small offset): {n_skipped}")


if __name__ == "__main__":
    main()
