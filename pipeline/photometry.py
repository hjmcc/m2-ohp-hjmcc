"""Photometric calibration: aperture photometry + PS1 cross-match → per-frame zero points."""

import csv
import json
import logging
from pathlib import Path

import numpy as np
import sep
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree

from . import config

log = logging.getLogger("pipeline")


# ── PS1 catalogue query ─────────────────────────────────────────────────

def query_ps1(center_ra, center_dec, radius_arcmin, ps1_band, cache_dir):
    """Query PS1 DR2 via VizieR TAP with local JSON file caching.

    Parameters
    ----------
    center_ra, center_dec : float
        Field centre in degrees (ICRS).
    radius_arcmin : float
        Cone-search radius in arcminutes.
    ps1_band : str
        PS1 column name, e.g. 'gmag', 'rmag', 'imag'.
    cache_dir : Path
        Directory for cached query results.

    Returns
    -------
    list of dict with keys: ra, dec, mag
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    ra_r = round(center_ra, 2)
    dec_r = round(center_dec, 2)
    rad_r = round(radius_arcmin, 1)
    cache_key = f"ps1_{ra_r}_{dec_r}_{rad_r}_{ps1_band}"
    cache_file = cache_dir / f"{cache_key}.json"

    if cache_file.exists():
        with open(cache_file) as f:
            stars = json.load(f)
        log.debug("PS1 cache hit: %d stars (%s)", len(stars), cache_file.name)
        return stars

    # VizieR TAP query for Pan-STARRS1 DR2
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    coord = SkyCoord(ra=center_ra, dec=center_dec, unit="deg")
    v = Vizier(
        columns=["RAJ2000", "DEJ2000", ps1_band, f"e_{ps1_band}"],
        column_filters={
            ps1_band: f"{config.PHOT_PS1_MAG_BRIGHT}..{config.PHOT_PS1_MAG_FAINT}",
        },
        row_limit=5000,
    )

    try:
        result = v.query_region(coord, radius=radius_arcmin * u.arcmin,
                                catalog="II/349/ps1")
    except Exception as exc:
        log.warning("PS1 VizieR query failed: %s", exc)
        return []

    if result is None or len(result) == 0:
        log.warning("PS1 query returned no results at RA=%.3f Dec=%.3f",
                    center_ra, center_dec)
        return []

    tab = result[0]
    stars = []
    for row in tab:
        mag = float(row[ps1_band])
        if np.isnan(mag):
            continue
        e_mag = float(row[f"e_{ps1_band}"]) if f"e_{ps1_band}" in tab.colnames else 0.0
        if np.isnan(e_mag):
            e_mag = 0.0
        stars.append({
            "ra": float(row["RAJ2000"]),
            "dec": float(row["DEJ2000"]),
            "mag": mag,
            "e_mag": e_mag,
        })

    # Cache result
    with open(cache_file, "w") as f:
        json.dump(stars, f)
    log.debug("PS1 query: %d stars at %s<%.1f, cached to %s",
              len(stars), ps1_band, config.PHOT_PS1_MAG_FAINT, cache_file.name)
    return stars


# ── Source extraction + aperture photometry ─────────────────────────────

def extract_photometry(data, header, telescope, fwhm_px=None):
    """Extract sources with sep and perform aperture photometry.

    Parameters
    ----------
    data : 2D array
        Reduced science frame.
    header : FITS header
    telescope : str
    fwhm_px : float or None
        Median FWHM in pixels from QC. Falls back to 4.0 if None.

    Returns
    -------
    dict with keys: x, y, ra, dec, flux, fluxerr, flag (all arrays)
        Returns None if WCS is not available or too few sources.
    """
    tel = config.TELESCOPES[telescope]
    if fwhm_px is None or not np.isfinite(fwhm_px):
        fwhm_px = 4.0

    work = np.ascontiguousarray(data, dtype=np.float64)
    bkg = sep.Background(work)
    data_sub = work - bkg

    objects = sep.extract(data_sub, config.PHOT_SEP_THRESH,
                          err=bkg.globalrms, minarea=config.SEP_MINAREA)
    if len(objects) == 0:
        return None

    # Filter: unflagged, not saturated, reasonable shape
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

    x = obj["x"]
    y = obj["y"]

    # Aperture photometry on background-subtracted data (no local annulus —
    # sep.Background already provides a smooth sky model)
    r_aper = config.PHOT_APERTURE_MULT * fwhm_px

    flux, fluxerr, flag = sep.sum_circle(data_sub, x, y, r_aper,
                                          err=bkg.globalrms, subpix=5)

    # Convert pixel → sky coords
    try:
        wcs = WCS(header)
        sky = wcs.pixel_to_world(x, y)
        ra = sky.ra.deg
        dec = sky.dec.deg
    except Exception:
        return None

    return {
        "x": x, "y": y,
        "ra": np.asarray(ra), "dec": np.asarray(dec),
        "flux": flux, "fluxerr": fluxerr,
        "flag": flag,
    }


# ── Zero-point fitting ──────────────────────────────────────────────────

def compute_zeropoint(inst_mags, cat_mags):
    """Compute zero point as sigma-clipped median of (cat_mag - inst_mag).

    Returns
    -------
    (zp, zp_err, n_used) or (None, None, 0) if insufficient stars.
    """
    diff = cat_mags - inst_mags
    n = len(diff)
    if n < config.PHOT_ZP_MIN_STARS:
        return None, None, 0

    # Iterative sigma clipping
    mask = np.ones(n, dtype=bool)
    for _ in range(5):
        med = np.median(diff[mask])
        mad = np.median(np.abs(diff[mask] - med))
        sigma = mad * 1.4826
        if sigma < 1e-6:
            break
        new_mask = mask & (np.abs(diff - med) < config.PHOT_ZP_SIGMA_CLIP * sigma)
        if np.sum(new_mask) < config.PHOT_ZP_MIN_STARS:
            break
        mask = new_mask

    n_used = int(np.sum(mask))
    if n_used < config.PHOT_ZP_MIN_STARS:
        return None, None, 0

    zp = float(np.median(diff[mask]))
    zp_err = float(np.std(diff[mask]) / np.sqrt(n_used))
    return zp, zp_err, n_used


# ── Per-frame orchestrator ──────────────────────────────────────────────

def calibrate_frame(filepath, telescope, ps1_cache_dir, force=False,
                    qc_stats=None):
    """Photometrically calibrate a single reduced frame.

    Parameters
    ----------
    filepath : Path
        Path to reduced FITS file.
    telescope : str
    ps1_cache_dir : Path
    force : bool
        Re-calibrate even if PHOTZP already present.
    qc_stats : dict or None
        QC entry for this frame (from frame_stats.json) for FWHM.

    Returns
    -------
    dict with calibration results for the report.
    """
    filepath = Path(filepath)
    hdr = fits.getheader(str(filepath))
    filt = hdr.get("FILTER", "")
    obj = hdr.get("OBJECT", "")
    exptime = hdr.get("EXPTIME", 0.0)
    night = filepath.parent.name

    result = {
        "filename": filepath.name,
        "night": night,
        "object": obj,
        "filter": filt,
        "exptime": exptime,
    }

    # Skip if already calibrated (unless forced)
    if not force and hdr.get("PHOTZP") is not None:
        result["status"] = "SKIPPED"
        result["zp"] = hdr["PHOTZP"]
        return result

    # Check astrometry status
    astr_stat = hdr.get("ASTRSTAT")
    if astr_stat not in ("OK", "HIGH_RMS"):
        result["status"] = "NO_WCS"
        return result

    # Filter → PS1 band mapping
    ps1_band = config.PHOT_FILTER_TO_PS1.get(filt)
    if ps1_band is None:
        result["status"] = "NO_PS1_BAND"
        result["ps1_band"] = None
        return result
    result["ps1_band"] = ps1_band

    # Get FWHM from QC stats
    fwhm_px = None
    if qc_stats and qc_stats.get("fwhm_px") is not None:
        fwhm_px = qc_stats["fwhm_px"]

    # Read data
    data = fits.getdata(str(filepath)).astype(np.float64)
    if not data.flags["C_CONTIGUOUS"]:
        data = np.ascontiguousarray(data)

    # Extract sources + aperture photometry
    phot = extract_photometry(data, hdr, telescope, fwhm_px=fwhm_px)
    if phot is None:
        result["status"] = "FEW_SOURCES"
        return result

    # Filter: positive flux, unflagged
    good = (phot["flux"] > 0) & (phot["flag"] == 0)
    if np.sum(good) < config.PHOT_ZP_MIN_STARS:
        result["status"] = "FEW_SOURCES"
        return result

    ra_src = phot["ra"][good]
    dec_src = phot["dec"][good]
    flux_src = phot["flux"][good]

    # Instrumental magnitudes: -2.5 * log10(flux / exptime)
    if exptime <= 0:
        result["status"] = "BAD_EXPTIME"
        return result
    inst_mag = -2.5 * np.log10(flux_src / exptime)

    # Query PS1
    wcs = WCS(hdr)
    ny, nx = data.shape
    center = wcs.pixel_to_world(nx / 2.0, ny / 2.0)
    tel = config.TELESCOPES[telescope]
    fov_arcmin = max(ny, nx) * tel["pixel_scale"] / 60.0
    query_radius = fov_arcmin * 0.75  # slightly larger than half-diagonal

    ps1_stars = query_ps1(center.ra.deg, center.dec.deg, query_radius,
                          ps1_band, ps1_cache_dir)
    if len(ps1_stars) < config.PHOT_ZP_MIN_STARS:
        result["status"] = "FEW_PS1"
        result["n_ps1"] = len(ps1_stars)
        return result

    # Cross-match with cKDTree
    cat_ra = np.array([s["ra"] for s in ps1_stars])
    cat_dec = np.array([s["dec"] for s in ps1_stars])
    cat_mag = np.array([s["mag"] for s in ps1_stars])

    match_radius_deg = config.PHOT_MATCH_RADIUS_ARCSEC / 3600.0

    # Build trees in (ra*cos(dec), dec) space for proper angular distances
    cos_dec = np.cos(np.radians(np.mean(dec_src)))
    tree_cat = cKDTree(np.column_stack([cat_ra * cos_dec, cat_dec]))
    tree_src = np.column_stack([ra_src * cos_dec, dec_src])

    dists, idxs = tree_cat.query(tree_src, distance_upper_bound=match_radius_deg)
    matched = np.isfinite(dists)
    n_matched = int(np.sum(matched))

    if n_matched < config.PHOT_ZP_MIN_STARS:
        result["status"] = "FEW_MATCHES"
        result["n_matched"] = n_matched
        return result

    matched_inst = inst_mag[matched]
    matched_cat = cat_mag[idxs[matched]]

    # Fit zero point
    zp, zp_err, n_used = compute_zeropoint(matched_inst, matched_cat)
    if zp is None:
        result["status"] = "ZP_FIT_FAILED"
        result["n_matched"] = n_matched
        return result

    # Write to FITS header
    with fits.open(str(filepath), mode="update") as hdul:
        h = hdul[0].header
        h["PHOTZP"] = (round(zp, 4), "Photometric zero point (mag)")
        h["PHOTZPER"] = (round(zp_err, 4), "Zero-point uncertainty (mag)")
        h["PHOTNSTR"] = (n_used, "N stars used for ZP fit")
        h["PHOTNMAT"] = (n_matched, "N cross-matched stars")
        h["PHOTCAT"] = (f"PS1_{ps1_band}", "Photometric reference catalogue")
        h["PHOTSTAT"] = ("OK", "Photometric calibration status")
        h["PHOTAPER"] = (round(config.PHOT_APERTURE_MULT * (fwhm_px or 4.0), 1),
                         "Aperture radius (px)")
        hdul.flush()

    result["status"] = "OK"
    result["zp"] = round(zp, 4)
    result["zp_err"] = round(zp_err, 4)
    result["n_matched"] = n_matched
    result["n_used"] = n_used
    return result


# ── Batch runner ────────────────────────────────────────────────────────

def run_photometry(year, telescope, force=False):
    """Run photometric calibration on all reduced frames for a telescope/year.

    Writes photometry/phot_report.json and phot_report.csv.
    """
    base = config.DATA_ROOT / year / telescope
    reduced_dir = base / "reduced"
    phot_dir = base / "photometry"
    ps1_cache_dir = phot_dir / "ps1_cache"

    if not reduced_dir.exists():
        log.warning("No reduced directory: %s", reduced_dir)
        return

    # Check if output already exists
    json_out = phot_dir / "phot_report.json"
    if json_out.exists() and not force:
        log.info("Photometry report already exists: %s (use --force to rerun)",
                 json_out)
        print(f"  Photometry report already exists for {telescope} — use --force to rerun")
        return

    # Load QC frame stats for FWHM lookup
    frame_stats_file = base / "qc" / "frame_stats.json"
    qc_lookup = {}
    if frame_stats_file.exists():
        with open(frame_stats_file) as f:
            for entry in json.load(f):
                qc_lookup[entry["filename"]] = entry

    # Collect all reduced FITS files
    fits_files = sorted(reduced_dir.rglob("*.fits")) + sorted(reduced_dir.rglob("*.fit"))
    if not fits_files:
        log.warning("No reduced FITS files in %s", reduced_dir)
        return

    print(f"\n  {telescope}: photometric calibration for {len(fits_files)} frames")

    results = []
    for i, fpath in enumerate(fits_files, 1):
        qc = qc_lookup.get(fpath.name)
        try:
            result = calibrate_frame(fpath, telescope, ps1_cache_dir,
                                     force=force, qc_stats=qc)
        except Exception as exc:
            log.error("Photometry failed for %s: %s", fpath.name, exc)
            result = {
                "filename": fpath.name,
                "night": fpath.parent.name,
                "status": f"ERROR: {exc}",
            }
        results.append(result)

        if i % 20 == 0:
            n_ok = sum(1 for r in results if r.get("status") == "OK")
            print(f"    [{i}/{len(fits_files)}] {n_ok} calibrated so far...")

    # Write reports
    phot_dir.mkdir(parents=True, exist_ok=True)

    with open(json_out, "w") as f:
        json.dump(results, f, indent=2, default=str)
    log.info("Wrote %s (%d frames)", json_out, len(results))

    csv_out = phot_dir / "phot_report.csv"
    if results:
        fieldnames = list(results[0].keys())
        # Collect all keys across results for uniform CSV
        for r in results:
            for k in r:
                if k not in fieldnames:
                    fieldnames.append(k)
        with open(csv_out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow({k: row.get(k, "") for k in fieldnames})
        log.info("Wrote %s", csv_out)

    # Console summary
    _print_summary(results, telescope)


def _print_summary(results, telescope):
    """Print per-filter photometry summary."""
    by_filter = {}
    for r in results:
        filt = r.get("filter", "?")
        by_filter.setdefault(filt, []).append(r)

    n_ok = sum(1 for r in results if r.get("status") == "OK")
    n_total = len(results)
    status_counts = {}
    for r in results:
        s = r.get("status", "?")
        status_counts[s] = status_counts.get(s, 0) + 1

    print(f"\n  ── {telescope} Photometry Summary: {n_ok}/{n_total} calibrated ──")
    print(f"  {'Filter':>10s}  {'N':>4s}  {'OK':>4s}  {'ZP':>7s}  {'σ_ZP':>6s}  {'N★':>5s}")
    print(f"  {'─'*10}  {'─'*4}  {'─'*4}  {'─'*7}  {'─'*6}  {'─'*5}")

    for filt in sorted(by_filter):
        frames = by_filter[filt]
        n = len(frames)
        ok_frames = [r for r in frames if r.get("status") == "OK"]
        n_ok_f = len(ok_frames)
        if ok_frames:
            zps = [r["zp"] for r in ok_frames]
            med_zp = f"{np.median(zps):.2f}"
            std_zp = f"{np.std(zps):.3f}" if len(zps) > 1 else "—"
            n_stars = [r.get("n_used", 0) for r in ok_frames]
            med_n = f"{int(np.median(n_stars))}"
        else:
            med_zp = "—"
            std_zp = "—"
            med_n = "—"
        print(f"  {filt:>10s}  {n:4d}  {n_ok_f:4d}  {med_zp:>7s}  {std_zp:>6s}  {med_n:>5s}")

    # Status breakdown
    non_ok = {k: v for k, v in status_counts.items() if k != "OK"}
    if non_ok:
        print(f"\n  Status breakdown:")
        for s, cnt in sorted(non_ok.items(), key=lambda x: -x[1]):
            print(f"    {s}: {cnt}")
    print()
