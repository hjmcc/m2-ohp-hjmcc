"""Quality control checks for bias, flat, and science frames."""

import csv
import json
import logging
from pathlib import Path

import numpy as np
import sep
from astropy.io import fits

from . import config

log = logging.getLogger("pipeline")


def check_bias_frames(bias_entries: list, telescope: str) -> dict:
    """QC-check bias frames. Returns {accepted, rejected, warnings}."""
    tel = config.TELESCOPES[telescope]
    accepted = []
    rejected = []
    warnings = []

    # First pass: reject on EXPTIME, flag temperature
    candidates = []
    for entry in bias_entries:
        path = entry["path"]
        hdr = fits.getheader(path)
        exptime = hdr.get("EXPTIME", 0.0)

        if exptime > 0:
            rejected.append({
                "filename": entry["filename"], "path": path,
                "reason": f"EXPTIME={exptime}s (test exposure, not bias)"
            })
            log.warning("REJECT bias %s: EXPTIME=%.1f", entry["filename"], exptime)
            continue

        ccd_temp = hdr.get("CCD-TEMP")
        temp_flag = None
        if ccd_temp is not None:
            delta = abs(ccd_temp - tel["nominal_temp"])
            if delta > tel["temp_tolerance"]:
                temp_flag = "TEMP_MISMATCH"
                warnings.append(
                    f"{entry['filename']}: CCD-TEMP={ccd_temp:.1f}°C, "
                    f"nominal={tel['nominal_temp']:.0f}°C (Δ={delta:.1f}°C)"
                )

        # Read data for statistics
        data = fits.getdata(path).astype(np.float32)
        med = float(np.median(data))
        std = float(np.std(data))

        candidates.append({
            "filename": entry["filename"], "path": path,
            "median": med, "stddev": std, "temp_flag": temp_flag,
            "ccd_temp": ccd_temp,
        })

    # Second pass: sigma-clip outliers by median level
    if len(candidates) >= 3:
        medians = np.array([c["median"] for c in candidates])
        group_med = np.median(medians)
        group_std = np.std(medians)
        clip = config.BIAS_SIGMA_CLIP

        for c in candidates:
            if group_std > 0 and abs(c["median"] - group_med) > clip * group_std:
                rejected.append({
                    "filename": c["filename"], "path": c["path"],
                    "reason": f"median={c['median']:.1f} deviates >{clip}σ from "
                              f"group median={group_med:.1f} (σ={group_std:.1f})"
                })
                log.warning("REJECT bias %s: outlier median=%.1f",
                            c["filename"], c["median"])
            else:
                accepted.append(c)
    else:
        accepted = candidates

    # Summary
    temp_mismatch = any(c.get("temp_flag") == "TEMP_MISMATCH" for c in accepted)
    if temp_mismatch:
        warnings.insert(0,
            f"ALL {telescope} biases have temperature mismatch — "
            f"readout pattern still usable but bias level may differ from science"
        )

    result = {
        "telescope": telescope,
        "n_input": len(bias_entries),
        "n_accepted": len(accepted),
        "n_rejected": len(rejected),
        "accepted": [{"filename": c["filename"], "median": c["median"],
                       "stddev": c["stddev"], "ccd_temp": c.get("ccd_temp"),
                       "temp_flag": c.get("temp_flag")}
                      for c in accepted],
        "rejected": rejected,
        "warnings": warnings,
        "bias_qual": "TEMP_MISMATCH" if temp_mismatch else "OK",
    }

    if accepted:
        meds = [c["median"] for c in accepted]
        result["group_median"] = float(np.median(meds))
        result["group_stddev"] = float(np.std(meds))

    log.info("Bias QC %s: %d/%d accepted, %d rejected, qual=%s",
             telescope, len(accepted), len(bias_entries),
             len(rejected), result["bias_qual"])
    return result


def check_flat_frames(flat_entries_by_filter: dict, telescope: str,
                      master_bias_data=None) -> dict:
    """QC-check flat frames per filter. Returns {filter: {accepted, rejected, ...}}."""
    results = {}

    for filt, entries in flat_entries_by_filter.items():
        accepted = []
        rejected = []

        # First pass: read data, compute stats after bias subtraction
        candidates = []
        for entry in entries:
            data = fits.getdata(entry["path"]).astype(np.float32)
            if master_bias_data is not None:
                data = data - master_bias_data

            med = float(np.median(data))

            # Signal range check
            if med < config.FLAT_MIN_ADU:
                rejected.append({
                    "filename": entry["filename"],
                    "reason": f"median={med:.0f} < {config.FLAT_MIN_ADU} (too low)"
                })
                log.warning("REJECT flat %s [%s]: low signal %.0f",
                            entry["filename"], filt, med)
                continue
            if med > config.FLAT_MAX_ADU:
                rejected.append({
                    "filename": entry["filename"],
                    "reason": f"median={med:.0f} > {config.FLAT_MAX_ADU} (too high)"
                })
                log.warning("REJECT flat %s [%s]: high signal %.0f",
                            entry["filename"], filt, med)
                continue

            # Gradient check: quadrant median ratio
            ny, nx = data.shape
            quads = [
                data[:ny//2, :nx//2], data[:ny//2, nx//2:],
                data[ny//2:, :nx//2], data[ny//2:, nx//2:],
            ]
            quad_meds = [float(np.median(q)) for q in quads]
            if min(quad_meds) > 0:
                ratio = max(quad_meds) / min(quad_meds)
            else:
                ratio = 999.0
            if ratio > config.FLAT_GRADIENT_MAX:
                rejected.append({
                    "filename": entry["filename"],
                    "reason": f"gradient ratio={ratio:.3f} > {config.FLAT_GRADIENT_MAX}"
                })
                log.warning("REJECT flat %s [%s]: gradient %.3f",
                            entry["filename"], filt, ratio)
                continue

            candidates.append({
                "filename": entry["filename"], "path": entry["path"],
                "median": med, "gradient_ratio": ratio,
            })

        # Second pass: sigma-clip within filter group
        if len(candidates) >= 3:
            medians = np.array([c["median"] for c in candidates])
            group_med = np.median(medians)
            group_std = np.std(medians)
            clip = config.FLAT_SIGMA_CLIP

            for c in candidates:
                if group_std > 0 and abs(c["median"] - group_med) > clip * group_std:
                    rejected.append({
                        "filename": c["filename"],
                        "reason": f"median={c['median']:.0f} deviates >{clip}σ from "
                                  f"group median={group_med:.0f}"
                    })
                    log.warning("REJECT flat %s [%s]: outlier", c["filename"], filt)
                else:
                    accepted.append(c)
        else:
            accepted = candidates

        results[filt] = {
            "n_input": len(entries),
            "n_accepted": len(accepted),
            "n_rejected": len(rejected),
            "accepted": [{"filename": c["filename"], "median": c["median"],
                           "gradient_ratio": c["gradient_ratio"]}
                          for c in accepted],
            "rejected": rejected,
        }
        log.info("Flat QC %s [%s]: %d/%d accepted",
                 telescope, filt, len(accepted), len(entries))

    return results


# ── Science frame QC ─────────────────────────────────────────────────────

def measure_frame_stats(path: str | Path, telescope: str) -> dict:
    """Measure QC statistics for a single reduced science frame.

    Uses sep for background estimation and source extraction.
    Returns a dict of per-frame measurements and quality flags.
    """
    path = Path(path)
    tel = config.TELESCOPES[telescope]
    pixel_scale = tel["pixel_scale"]
    saturation = tel["saturation"]

    hdr = fits.getheader(path)
    data = fits.getdata(path).astype(np.float64)

    # Ensure C-contiguous for sep
    if not data.flags["C_CONTIGUOUS"]:
        data = np.ascontiguousarray(data)

    # ── Image stats ──────────────────────────────────────────────────
    img_median = float(np.median(data))
    # Robust sigma via MAD
    mad = float(np.median(np.abs(data - img_median)))
    img_sigma = mad * 1.4826
    img_min = float(np.min(data))
    img_max = float(np.max(data))
    sat_threshold = 0.9 * saturation
    saturated_frac = float(np.sum(data > sat_threshold)) / data.size

    # ── sep background + source extraction ───────────────────────────
    bkg = sep.Background(data)
    data_sub = data - bkg
    objects = sep.extract(data_sub, config.SEP_THRESH, err=bkg.globalrms,
                          minarea=config.SEP_MINAREA)

    n_sources = len(objects)

    # ── FWHM and elongation from extracted sources ───────────────────
    fwhm_px = np.nan
    fwhm_arcsec = np.nan
    elongation = np.nan
    n_good = 0

    if n_sources > 0:
        # Compute half-light radii (kronrad not needed — use flux_radius)
        x = objects["x"]
        y = objects["y"]
        a = objects["a"]
        b = objects["b"]
        theta = objects["theta"]

        # Half-light radius for each source
        r_half, flag_r = sep.flux_radius(
            data_sub, x, y, 6.0 * a, 0.5,
            subpix=5,
        )
        fwhm_all = 2.0 * r_half

        # Elongation per source
        elong_all = np.where(b > 0, a / b, np.nan)

        # Filter "good" sources: unflagged, unsaturated peak, FWHM in range
        peak = objects["peak"]
        src_flag = objects["flag"]
        good_mask = (
            (src_flag == 0)
            & (flag_r == 0)
            & (peak < sat_threshold)
            & (fwhm_all >= config.FWHM_MIN_PX)
            & (fwhm_all <= config.FWHM_MAX_PX)
            & np.isfinite(fwhm_all)
            & np.isfinite(elong_all)
        )
        n_good = int(np.sum(good_mask))

        if n_good > 0:
            fwhm_px = float(np.median(fwhm_all[good_mask]))
            fwhm_arcsec = fwhm_px * pixel_scale
            elongation = float(np.median(elong_all[good_mask]))

    # ── Extract night from path ──────────────────────────────────────
    # reduced/{night}/filename.fits
    night = path.parent.name

    # ── Flagging ─────────────────────────────────────────────────────
    flags = []
    if n_good < config.QC_MIN_GOOD_SOURCES:
        flags.append("LOW_SOURCES")
    if np.isfinite(fwhm_arcsec) and fwhm_arcsec > config.QC_MAX_FWHM_ARCSEC:
        flags.append("BAD_SEEING")
    if np.isfinite(elongation) and elongation > config.QC_MAX_ELONGATION:
        flags.append("ELONGATED")
    if saturated_frac > config.QC_MAX_SATURATED_FRAC:
        flags.append("SATURATED")

    flag_str = ",".join(flags) if flags else "OK"

    return {
        "filename": path.name,
        "night": night,
        "object": hdr.get("OBJECT", ""),
        "filter": hdr.get("FILTER", ""),
        "exptime": hdr.get("EXPTIME", 0.0),
        "airmass": hdr.get("AIRMASS"),
        "sky_median": round(img_median, 1),
        "sky_sigma": round(img_sigma, 1),
        "n_sources": n_sources,
        "n_good": n_good,
        "fwhm_px": round(fwhm_px, 2) if np.isfinite(fwhm_px) else None,
        "fwhm_arcsec": round(fwhm_arcsec, 2) if np.isfinite(fwhm_arcsec) else None,
        "elongation": round(elongation, 3) if np.isfinite(elongation) else None,
        "peak_adu": round(img_max, 1),
        "saturated_frac": round(saturated_frac, 6),
        "flag": flag_str,
    }


def run_science_qc(year: str, telescope: str, force: bool = False):
    """Run QC stats on all reduced science frames for a telescope/year.

    Writes qc/frame_stats.json and qc/frame_stats.csv.
    Prints a summary table to console.
    """
    base = config.DATA_ROOT / year / telescope
    reduced_dir = base / "reduced"
    qc_dir = base / "qc"

    if not reduced_dir.exists():
        log.warning("No reduced directory: %s", reduced_dir)
        return

    # Check if output exists and skip unless forced
    json_out = qc_dir / "frame_stats.json"
    csv_out = qc_dir / "frame_stats.csv"
    if json_out.exists() and not force:
        log.info("QC stats already exist: %s (use --force to rerun)", json_out)
        print(f"  QC stats already exist for {telescope} — use --force to rerun")
        return

    # Collect all reduced FITS files
    fits_files = sorted(reduced_dir.rglob("*.fits")) + sorted(reduced_dir.rglob("*.fit"))
    if not fits_files:
        log.warning("No reduced FITS files in %s", reduced_dir)
        return

    print(f"\n  {telescope}: measuring QC stats for {len(fits_files)} reduced frames...")

    # Process each frame
    all_stats = []
    for i, fpath in enumerate(fits_files, 1):
        log.debug("QC [%d/%d] %s", i, len(fits_files), fpath.name)
        try:
            stats = measure_frame_stats(fpath, telescope)
            all_stats.append(stats)
        except Exception as exc:
            log.error("QC failed for %s: %s", fpath.name, exc)
            all_stats.append({
                "filename": fpath.name,
                "night": fpath.parent.name,
                "flag": f"ERROR: {exc}",
            })
        if i % 50 == 0:
            print(f"    {i}/{len(fits_files)} done...")

    # Write outputs
    qc_dir.mkdir(parents=True, exist_ok=True)

    with open(json_out, "w") as f:
        json.dump(all_stats, f, indent=2, default=str)
    log.info("Wrote %s (%d frames)", json_out, len(all_stats))

    # CSV
    if all_stats:
        fieldnames = list(all_stats[0].keys())
        with open(csv_out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in all_stats:
                writer.writerow({k: row.get(k, "") for k in fieldnames})
        log.info("Wrote %s", csv_out)

    # ── Console summary ──────────────────────────────────────────────
    _print_qc_summary(all_stats, telescope)


def _print_qc_summary(all_stats: list, telescope: str):
    """Print per-filter summary table to console."""
    # Group by filter
    by_filter = {}
    for s in all_stats:
        filt = s.get("filter", "?")
        by_filter.setdefault(filt, []).append(s)

    n_flagged = sum(1 for s in all_stats if s.get("flag", "OK") != "OK")
    n_total = len(all_stats)

    print(f"\n  ── {telescope} QC Summary: {n_total} frames, {n_flagged} flagged ──")
    print(f"  {'Filter':>10s}  {'N':>4s}  {'FWHM″':>7s}  {'Sky':>7s}  {'Sources':>7s}  {'Flagged':>7s}")
    print(f"  {'─'*10}  {'─'*4}  {'─'*7}  {'─'*7}  {'─'*7}  {'─'*7}")

    for filt in sorted(by_filter):
        frames = by_filter[filt]
        n = len(frames)
        fwhms = [s["fwhm_arcsec"] for s in frames
                 if s.get("fwhm_arcsec") is not None]
        skys = [s["sky_median"] for s in frames
                if s.get("sky_median") is not None]
        srcs = [s["n_good"] for s in frames
                if s.get("n_good") is not None]
        nflag = sum(1 for s in frames if s.get("flag", "OK") != "OK")

        med_fwhm = f"{np.median(fwhms):.2f}" if fwhms else "—"
        med_sky = f"{np.median(skys):.0f}" if skys else "—"
        med_src = f"{np.median(srcs):.0f}" if srcs else "—"

        print(f"  {filt:>10s}  {n:4d}  {med_fwhm:>7s}  {med_sky:>7s}  {med_src:>7s}  {nflag:7d}")

    # Show flag breakdown
    flag_counts = {}
    for s in all_stats:
        flag = s.get("flag", "OK")
        if flag != "OK":
            for f in flag.split(","):
                flag_counts[f] = flag_counts.get(f, 0) + 1

    if flag_counts:
        print(f"\n  Flag breakdown:")
        for f, cnt in sorted(flag_counts.items(), key=lambda x: -x[1]):
            print(f"    {f}: {cnt}")
    print()
