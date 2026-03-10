#!/usr/bin/env python3
"""Generate per-frame photometric calibration diagnostic plots for 2025 OHP data.

For each frame with PHOTSTAT=OK, produces a scatter plot of
PS1 catalog mag vs instrumental mag with the fitted ZP line.

Plots saved to plots/photometry/{telescope}/{night}/ as PNGs.
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sep
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pipeline import config

YEAR = "2025"
PLOT_DIR = Path(__file__).resolve().parent / "plots" / "photometry"


def extract_and_match(filepath, telescope, ps1_cache_dir, qc_stats=None):
    """Re-run source extraction + PS1 cross-match for a single frame.

    Returns dict with per-star arrays + frame metadata, or None on failure.
    """
    hdr = fits.getheader(str(filepath))
    exptime = hdr.get("EXPTIME", 0.0)
    filt = hdr.get("FILTER", "")
    zp = hdr.get("PHOTZP")
    zp_err = hdr.get("PHOTZPER", 0.0)
    n_used = hdr.get("PHOTNSTR", 0)

    if exptime <= 0 or zp is None:
        return None

    ps1_band = config.PHOT_FILTER_TO_PS1.get(filt)
    if ps1_band is None:
        return None

    fwhm_px = 4.0
    if qc_stats and qc_stats.get("fwhm_px") is not None:
        fwhm_px = qc_stats["fwhm_px"]

    data = fits.getdata(str(filepath)).astype(np.float64)
    if not data.flags["C_CONTIGUOUS"]:
        data = np.ascontiguousarray(data)

    tel = config.TELESCOPES[telescope]
    bkg = sep.Background(data)
    data_sub = data - bkg
    try:
        objects = sep.extract(data_sub, config.PHOT_SEP_THRESH,
                              err=bkg.globalrms, minarea=config.SEP_MINAREA)
    except Exception:
        return None

    if len(objects) == 0:
        return None

    sat = 0.9 * tel["saturation"]
    mask = (
        (objects["flag"] == 0)
        & (objects["peak"] < sat)
        & (objects["a"] > 0) & (objects["b"] > 0)
        & (objects["a"] / objects["b"] < 3.0)
    )
    obj = objects[mask]
    if len(obj) < config.PHOT_ZP_MIN_STARS:
        return None

    x, y = obj["x"], obj["y"]
    r_aper = config.PHOT_APERTURE_MULT * fwhm_px
    flux, fluxerr, flag = sep.sum_circle(data_sub, x, y, r_aper,
                                          err=bkg.globalrms, subpix=5)

    try:
        wcs = WCS(hdr)
        sky = wcs.pixel_to_world(x, y)
        ra_src = sky.ra.deg
        dec_src = sky.dec.deg
    except Exception:
        return None

    good = (flux > 0) & (flag == 0)
    if np.sum(good) < config.PHOT_ZP_MIN_STARS:
        return None

    ra_src = np.asarray(ra_src)[good]
    dec_src = np.asarray(dec_src)[good]
    flux_src = flux[good]
    inst_mag = -2.5 * np.log10(flux_src / exptime)

    # Load cached PS1
    ny, nx = data.shape
    center = wcs.pixel_to_world(nx / 2.0, ny / 2.0)
    fov_arcmin = max(ny, nx) * tel["pixel_scale"] / 60.0
    query_radius = fov_arcmin * 0.75

    ra_r = round(center.ra.deg, 2)
    dec_r = round(center.dec.deg, 2)
    rad_r = round(query_radius, 1)
    cache_key = f"ps1_{ra_r}_{dec_r}_{rad_r}_{ps1_band}"
    cache_file = ps1_cache_dir / f"{cache_key}.json"

    if not cache_file.exists():
        return None

    with open(cache_file) as f:
        ps1_stars = json.load(f)

    if len(ps1_stars) < config.PHOT_ZP_MIN_STARS:
        return None

    cat_ra = np.array([s["ra"] for s in ps1_stars])
    cat_dec = np.array([s["dec"] for s in ps1_stars])
    cat_mag = np.array([s["mag"] for s in ps1_stars])

    match_radius_deg = config.PHOT_MATCH_RADIUS_ARCSEC / 3600.0
    cos_dec = np.cos(np.radians(np.mean(dec_src)))
    tree_cat = cKDTree(np.column_stack([cat_ra * cos_dec, cat_dec]))
    tree_src = np.column_stack([ra_src * cos_dec, dec_src])
    dists, idxs = tree_cat.query(tree_src, distance_upper_bound=match_radius_deg)
    matched = np.isfinite(dists)

    if np.sum(matched) < config.PHOT_ZP_MIN_STARS:
        return None

    m_inst = inst_mag[matched]
    m_cat = cat_mag[idxs[matched]]

    return {
        "inst_mag": m_inst,
        "cat_mag": m_cat,
        "zp": zp,
        "zp_err": zp_err,
        "n_used": n_used,
        "ps1_band": ps1_band,
        "filter": filt,
        "object": hdr.get("OBJECT", ""),
        "exptime": exptime,
        "fwhm_px": fwhm_px,
        "aper_px": r_aper,
    }


def make_frame_plot(tel, night, filename, match_data, outdir):
    """One scatter plot per frame: catalog mag vs instrumental mag."""
    m_inst = match_data["inst_mag"]
    m_cat = match_data["cat_mag"]
    zp = match_data["zp"]
    zp_err = match_data["zp_err"]
    n_used = match_data["n_used"]
    ps1_band = match_data["ps1_band"]
    filt = match_data["filter"]
    obj = match_data["object"]
    exptime = match_data["exptime"]
    fwhm = match_data["fwhm_px"]
    aper = match_data["aper_px"]

    # Sigma-clip to identify inliers (reproduce pipeline logic)
    diff = m_cat - m_inst
    inlier = np.ones(len(diff), dtype=bool)
    for _ in range(5):
        med = np.median(diff[inlier])
        mad = np.median(np.abs(diff[inlier] - med))
        sigma = mad * 1.4826
        if sigma < 1e-6:
            break
        new_inlier = np.abs(diff - med) < config.PHOT_ZP_SIGMA_CLIP * sigma
        if np.sum(new_inlier) < config.PHOT_ZP_MIN_STARS:
            break
        inlier = new_inlier

    n_inlier = int(np.sum(inlier))
    resid_in = diff[inlier] - zp
    rms = np.std(resid_in) if n_inlier > 1 else 0.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    # --- Left panel: catalog mag vs instrumental mag ---
    ax1.scatter(m_inst[~inlier], m_cat[~inlier], s=80, alpha=0.5,
                color="C3", edgecolors="k", linewidths=0.3, zorder=2,
                label=f"clipped ({np.sum(~inlier)})")
    ax1.scatter(m_inst[inlier], m_cat[inlier], s=80, alpha=0.7,
                color="C0", edgecolors="k", linewidths=0.3, zorder=3,
                label=f"used for ZP ({n_inlier})")

    # ZP line: cat = inst + ZP
    x_lo, x_hi = np.min(m_inst) - 0.5, np.max(m_inst) + 0.5
    ax1.plot([x_lo, x_hi], [x_lo + zp, x_hi + zp], "r-", lw=2, alpha=0.8,
             zorder=4, label=f"ZP = {zp:.3f} ± {zp_err:.3f}")

    ax1.set_xlabel("Instrumental mag  $[-2.5\\,\\log_{10}(\\mathrm{flux}/t_\\mathrm{exp})]$",
                    fontsize=13)
    ax1.set_ylabel(f"PS1 {ps1_band}  (catalog)", fontsize=13)
    ax1.legend(fontsize=11, loc="upper left")
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(labelsize=11)

    # --- Right panel: residuals vs catalog mag ---
    ax2.scatter(m_cat[~inlier], diff[~inlier] - zp, s=80, alpha=0.5,
                color="C3", edgecolors="k", linewidths=0.3, zorder=2)
    ax2.scatter(m_cat[inlier], diff[inlier] - zp, s=80, alpha=0.7,
                color="C0", edgecolors="k", linewidths=0.3, zorder=3)
    ax2.axhline(0, color="r", lw=2, alpha=0.7)
    if n_inlier > 1:
        ax2.axhline(+rms, color="r", lw=1, ls="--", alpha=0.5)
        ax2.axhline(-rms, color="r", lw=1, ls="--", alpha=0.5,
                     label=f"σ = {rms:.3f} mag")
    ylim = max(3.0, 1.5 * np.max(np.abs(diff - zp)))
    ax2.set_ylim(-ylim, ylim)

    ax2.set_xlabel(f"PS1 {ps1_band}  (catalog)", fontsize=13)
    ax2.set_ylabel("Residual  (cat − inst − ZP)  [mag]", fontsize=13)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(labelsize=11)

    fig.suptitle(
        f"{tel}  {night}  {filename}\n"
        f"{obj}  |  filter {filt}  |  {exptime:.0f}s  |  "
        f"FWHM {fwhm:.1f} px  |  aper {aper:.1f} px  |  "
        f"{len(m_inst)} matched → {n_inlier} used  |  "
        f"ZP = {zp:.3f} ± {zp_err:.3f}  |  σ = {rms:.3f}",
        fontsize=13, fontweight="bold", y=1.02)

    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    stem = Path(filename).stem
    outfile = outdir / f"{stem}_phot.png"
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return outfile.name


def main():
    print("Per-frame photometric calibration plots — 2025 data\n")

    for telescope in ["T080", "T120"]:
        base = config.DATA_ROOT / YEAR / telescope
        reduced_dir = base / "reduced"
        ps1_cache_dir = base / "photometry" / "ps1_cache"
        report_file = base / "photometry" / "phot_report.json"

        if not report_file.exists():
            print(f"  {telescope}: no photometry report, skipping")
            continue

        with open(report_file) as f:
            report = json.load(f)
        ok_frames = {r["filename"]: r for r in report if r.get("status") == "OK"}

        qc_file = base / "qc" / "frame_stats.json"
        qc_lookup = {}
        if qc_file.exists():
            with open(qc_file) as f:
                for entry in json.load(f):
                    qc_lookup[entry["filename"]] = entry

        print(f"  {telescope}: {len(ok_frames)} calibrated frames")

        n_plotted = 0
        for fname, rep in ok_frames.items():
            night = rep["night"]
            fpath = reduced_dir / night / fname
            if not fpath.exists():
                continue

            qc = qc_lookup.get(fname)
            match_data = extract_and_match(fpath, telescope, ps1_cache_dir, qc_stats=qc)
            if match_data is None:
                continue

            outdir = PLOT_DIR / telescope / night
            pname = make_frame_plot(telescope, night, fname, match_data, outdir)
            n_plotted += 1

            if n_plotted % 25 == 0:
                print(f"    [{n_plotted}] ...")

        print(f"    {n_plotted} plots saved under {PLOT_DIR / telescope}/\n")

    print(f"Done. All plots in {PLOT_DIR}/")


if __name__ == "__main__":
    main()
