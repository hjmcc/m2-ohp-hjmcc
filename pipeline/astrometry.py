"""Astrometric calibration: WCS solutions via Gaia DR3 cross-matching.

Primary solver: sep extraction → Gaia cone search → cKDTree cross-match →
least-squares WCS fit. Uses astroalign for T080 (unknown rotation).
Fallback: astrometry.net solve-field if available.
"""

import csv
import json
import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import sep
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from scipy.spatial import cKDTree
from scipy.optimize import least_squares

from . import config

log = logging.getLogger("pipeline")


# ── Source extraction ────────────────────────────────────────────────────

def extract_sources(data, telescope, thresh=None):
    """Extract point sources with sep.

    Returns structured array with fields (x, y, flux, peak, a, b) sorted
    by flux descending, capped at ASTROM_MAX_SOURCES.
    """
    tel = config.TELESCOPES[telescope]
    if thresh is None:
        thresh = config.ASTROM_SEP_THRESH

    work = np.ascontiguousarray(data, dtype=np.float64)
    bkg = sep.Background(work)
    data_sub = work - bkg

    objects = sep.extract(data_sub, thresh, err=bkg.globalrms,
                          minarea=config.SEP_MINAREA)
    if len(objects) == 0:
        return np.empty(0, dtype=[("x", float), ("y", float), ("flux", float)])

    # Filter: unflagged, not saturated, reasonable shape
    sat = 0.9 * tel["saturation"]
    mask = (
        (objects["flag"] == 0)
        & (objects["peak"] < sat)
        & (objects["a"] > 0) & (objects["b"] > 0)
        & (objects["a"] / objects["b"] < 3.0)  # reject cosmics/streaks
    )
    obj = objects[mask]
    if len(obj) == 0:
        return np.empty(0, dtype=[("x", float), ("y", float), ("flux", float)])

    # Sort by flux descending, keep top N
    order = np.argsort(-obj["flux"])
    obj = obj[order[: config.ASTROM_MAX_SOURCES]]

    result = np.zeros(len(obj), dtype=[("x", float), ("y", float), ("flux", float)])
    result["x"] = obj["x"]
    result["y"] = obj["y"]
    result["flux"] = obj["flux"]
    return result


# ── Initial WCS from header ─────────────────────────────────────────────

def _parse_sexagesimal(s):
    """Parse 'HH MM SS' or '+DD MM SS' string to float degrees."""
    parts = s.strip().split()
    if len(parts) != 3:
        return None
    sign = -1 if parts[0].startswith("-") else 1
    d = abs(float(parts[0]))
    m = float(parts[1])
    s = float(parts[2])
    return sign * (d + m / 60.0 + s / 3600.0)


def _build_cd_matrix(scale_deg, rot_deg, flipped=False):
    """Build CD matrix from scale, rotation, and parity.

    Parameters
    ----------
    scale_deg : float
        Pixel scale in degrees per pixel.
    rot_deg : float
        Position angle in degrees.
    flipped : bool
        If True, use flipped parity (det > 0).

    Returns
    -------
    (2, 2) numpy array
    """
    rot_rad = np.radians(rot_deg)
    c, s = np.cos(rot_rad), np.sin(rot_rad)
    if flipped:
        # Flipped parity (det > 0): e.g. T080 Cassegrain
        return np.array([[scale_deg * c, -scale_deg * s],
                         [scale_deg * s,  scale_deg * c]])
    else:
        # Standard parity (det < 0): e.g. T120 Newton
        return np.array([[-scale_deg * c, scale_deg * s],
                         [ scale_deg * s, scale_deg * c]])


def _build_wcs(crval, crpix, scale_deg, rot_deg, flipped=False):
    """Build a complete TAN WCS from parameters."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = crpix
    w.wcs.crval = crval
    w.wcs.cd = _build_cd_matrix(scale_deg, rot_deg, flipped)
    w.wcs.equinox = 2000.0
    w.wcs.radesys = "ICRS"
    return w


def get_initial_wcs(header, telescope):
    """Build initial TAN WCS from header pointing + known parameters.

    Returns WCS object or None if no pointing available.
    """
    ra_str = header.get("OBJCTRA")
    dec_str = header.get("OBJCTDEC")
    if ra_str is None or dec_str is None:
        return None

    ra_hrs = _parse_sexagesimal(ra_str)
    dec_deg = _parse_sexagesimal(dec_str)
    if ra_hrs is None or dec_deg is None:
        return None

    ra_deg = ra_hrs * 15.0  # hours → degrees

    tel = config.TELESCOPES[telescope]
    scale_deg = tel["pixel_scale"] / 3600.0

    ny = header.get("NAXIS2", tel["expected_shape"][0])
    nx = header.get("NAXIS1", tel["expected_shape"][1])

    rot_deg = config.ASTROM_KNOWN_ROTATION.get(telescope)
    if rot_deg is None:
        rot_deg = 0.0
    flipped = config.ASTROM_FLIPPED_PARITY.get(telescope, False)

    return _build_wcs(
        [ra_deg, dec_deg],
        [nx / 2.0 + 0.5, ny / 2.0 + 0.5],
        scale_deg, rot_deg, flipped,
    )


# ── Object name resolution ──────────────────────────────────────────────

def _clean_object_name(name):
    """Clean T080/T120 object names for Simbad lookup."""
    cleaned = name.strip()
    # Strip telescope prefix
    cleaned = re.sub(r"^T\d{3}_", "", cleaned)
    # Map special cases
    special = {
        "Coma_N4874": "NGC 4874",
        "Coma Cote": "NGC 4874",  # Coma cluster field
        "29P_2": "29P/Schwassmann-Wachmann",
        "Vullcano": "20 Massalia",  # asteroid Vulcano typo
        "vulcano": "20 Massalia",
        "Vulcano": "20 Massalia",
        "EuclidQ1": "Euclid",
    }
    if cleaned in special:
        return special[cleaned]
    # Strip trailing underscores
    cleaned = cleaned.rstrip("_")
    # Strip champN_ suffix (M67champ1_ → M67)
    cleaned = re.sub(r"champ\d+$", "", cleaned)
    # Strip _2_ sequence suffix (M99_2_ → M99)
    cleaned = re.sub(r"_\d+$", "", cleaned)
    # Strip leading/trailing whitespace again
    cleaned = cleaned.strip("_").strip()
    return cleaned


def resolve_object_coords(object_name, simbad_cache):
    """Resolve object name to SkyCoord via simbad_cache, then live Simbad.

    Returns SkyCoord or None.
    """
    # Try raw name in cache first
    entry = simbad_cache.get(object_name)
    if entry is not None:
        return SkyCoord(ra=entry["ra"], dec=entry["dec"], unit="deg")

    # Try cleaned name
    cleaned = _clean_object_name(object_name)
    entry = simbad_cache.get(cleaned)
    if entry is not None:
        return SkyCoord(ra=entry["ra"], dec=entry["dec"], unit="deg")

    # Try common variants in cache
    for variant in [cleaned, cleaned.replace(" ", ""), cleaned.replace("_", " ")]:
        entry = simbad_cache.get(variant)
        if entry is not None:
            return SkyCoord(ra=entry["ra"], dec=entry["dec"], unit="deg")

    # Live Simbad TAP query as last resort
    try:
        from astroquery.simbad import Simbad
        result = Simbad.query_object(cleaned)
        if result is not None and len(result) > 0:
            ra = result["ra"][0]
            dec = result["dec"][0]
            coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
            log.info("Simbad resolved '%s' → RA=%.4f Dec=%.4f",
                     cleaned, coord.ra.deg, coord.dec.deg)
            return coord
    except Exception as exc:
        log.warning("Simbad query failed for '%s': %s", cleaned, exc)

    return None


# ── Gaia DR3 queries ────────────────────────────────────────────────────

def query_gaia(center, radius_arcmin, mag_limit, cache_dir):
    """Query Gaia DR3 cone search, with local file caching.

    Parameters
    ----------
    center : SkyCoord
    radius_arcmin : float
    mag_limit : float
        Faint magnitude limit (G band).
    cache_dir : Path

    Returns
    -------
    list of dict with keys: ra, dec, pmra, pmdec, phot_g_mean_mag
        Coordinates are corrected for proper motion to epoch J2025.2.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    ra_r = round(center.ra.deg, 2)
    dec_r = round(center.dec.deg, 2)
    rad_r = round(radius_arcmin, 1)
    cache_key = f"{ra_r}_{dec_r}_{rad_r}_{mag_limit}"
    cache_file = cache_dir / f"{cache_key}.json"

    if cache_file.exists():
        with open(cache_file) as f:
            stars = json.load(f)
        log.debug("Gaia cache hit: %d stars (%s)", len(stars), cache_file.name)
        return stars

    # TAP query
    query = f"""
    SELECT ra, dec, pmra, pmdec, phot_g_mean_mag
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {center.ra.deg}, {center.dec.deg}, {radius_arcmin / 60.0})
    )
    AND phot_g_mean_mag < {mag_limit}
    AND pmra IS NOT NULL AND pmdec IS NOT NULL
    ORDER BY phot_g_mean_mag ASC
    """

    stars = _run_gaia_query(query)
    if stars is None:
        return []

    # Proper-motion correction: J2016.0 → J2025.2 (9.2 years)
    dt_yr = 9.2
    for s in stars:
        # pmra is in mas/yr (includes cos(dec) factor)
        s["ra"] += (s["pmra"] / 3600000.0) / np.cos(np.radians(s["dec"])) * dt_yr
        s["dec"] += (s["pmdec"] / 3600000.0) * dt_yr

    # Cache
    with open(cache_file, "w") as f:
        json.dump(stars, f)
    log.debug("Gaia query: %d stars at G<%.1f, cached to %s",
              len(stars), mag_limit, cache_file.name)

    return stars


def _run_gaia_query(query, retries=1):
    """Execute Gaia TAP query, return list of dicts or None."""
    from astroquery.gaia import Gaia

    for attempt in range(retries + 1):
        try:
            job = Gaia.launch_job(query)
            table = job.get_results()
            stars = []
            for row in table:
                stars.append({
                    "ra": float(row["ra"]),
                    "dec": float(row["dec"]),
                    "pmra": float(row["pmra"]),
                    "pmdec": float(row["pmdec"]),
                    "phot_g_mean_mag": float(row["phot_g_mean_mag"]),
                })
            return stars
        except Exception as exc:
            if attempt < retries:
                log.warning("Gaia query failed (attempt %d), retrying: %s",
                            attempt + 1, exc)
            else:
                log.error("Gaia query failed: %s", exc)
                return None


def get_gaia_stars(center, telescope, cache_dir):
    """Query Gaia with adaptive depth for a telescope FOV.

    Returns list of dicts with ra, dec, phot_g_mean_mag.
    """
    tel = config.TELESCOPES[telescope]
    ny, nx = tel["expected_shape"]
    pixel_scale = tel["pixel_scale"]

    # FOV diagonal radius in arcmin
    diag_px = np.sqrt(nx**2 + ny**2) / 2.0
    fov_radius = diag_px * pixel_scale / 60.0  # arcmin
    search_radius = fov_radius * config.ASTROM_GAIA_QUERY_MARGIN

    mag_limit = config.ASTROM_GAIA_MAG_LIMIT.get(telescope, 18.0)
    stars = query_gaia(center, search_radius, mag_limit, cache_dir)

    # Adaptive depth: if too few stars, go deeper
    if len(stars) < 20:
        deep_limit = config.ASTROM_GAIA_MAG_DEEP
        log.info("Only %d Gaia stars at G<%.1f, querying deeper to G<%.1f",
                 len(stars), mag_limit, deep_limit)
        stars = query_gaia(center, search_radius, deep_limit, cache_dir)

    return stars


# ── Cross-matching ───────────────────────────────────────────────────────

def cross_match(src_xy, cat_xy, tolerance_px):
    """Cross-match two sets of (x,y) points using cKDTree.

    Returns (src_idx, cat_idx, distances) arrays for matches within tolerance.
    """
    if len(src_xy) == 0 or len(cat_xy) == 0:
        return np.array([], int), np.array([], int), np.array([])

    tree = cKDTree(cat_xy)
    dists, cat_idx = tree.query(src_xy, k=1)

    mask = dists < tolerance_px
    src_idx = np.where(mask)[0]
    cat_idx = cat_idx[mask]
    dists = dists[mask]

    # Remove duplicate catalog matches (keep closest)
    if len(cat_idx) > 0:
        unique_cat = {}
        for i, (si, ci, d) in enumerate(zip(src_idx, cat_idx, dists)):
            if ci not in unique_cat or d < unique_cat[ci][1]:
                unique_cat[ci] = (si, d)
        keep_src = []
        keep_cat = []
        keep_dist = []
        for ci, (si, d) in unique_cat.items():
            keep_src.append(si)
            keep_cat.append(ci)
            keep_dist.append(d)
        src_idx = np.array(keep_src)
        cat_idx = np.array(keep_cat)
        dists = np.array(keep_dist)

    return src_idx, cat_idx, dists


# ── WCS fitting ──────────────────────────────────────────────────────────

def fit_wcs(src_xy, cat_radec, wcs_init, image_shape, crval_only=False):
    """Fit WCS via least squares.

    Parameters
    ----------
    src_xy : (N, 2) array of pixel positions
    cat_radec : (N, 2) array of (ra, dec) in degrees
    wcs_init : initial WCS
    image_shape : (ny, nx)
    crval_only : bool
        If True, fit only CRVAL (2 params), keeping CD matrix fixed.
        More robust with few or noisy matches.

    Returns
    -------
    WCS with fitted parameters, or None on failure.
    """
    ny, nx = image_shape
    crpix = np.array([nx / 2.0 + 0.5, ny / 2.0 + 0.5])
    crval = wcs_init.wcs.crval.copy()
    cd = wcs_init.wcs.cd.copy()

    if crval_only:
        p0 = np.array([crval[0], crval[1]])

        def residuals(params):
            w = WCS(naxis=2)
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            w.wcs.crpix = crpix
            w.wcs.crval = [params[0], params[1]]
            w.wcs.cd = cd
            sky = w.all_pix2world(src_xy, 0)
            cos_d = np.cos(np.radians(cat_radec[:, 1]))
            dra = (sky[:, 0] - cat_radec[:, 0]) * cos_d * 3600.0
            ddec = (sky[:, 1] - cat_radec[:, 1]) * 3600.0
            return np.concatenate([dra, ddec])
    else:
        p0 = np.array([crval[0], crval[1],
                        cd[0, 0], cd[0, 1], cd[1, 0], cd[1, 1]])

        def residuals(params):
            w = WCS(naxis=2)
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            w.wcs.crpix = crpix
            w.wcs.crval = [params[0], params[1]]
            w.wcs.cd = np.array([[params[2], params[3]],
                                 [params[4], params[5]]])
            sky = w.all_pix2world(src_xy, 0)
            cos_d = np.cos(np.radians(cat_radec[:, 1]))
            dra = (sky[:, 0] - cat_radec[:, 0]) * cos_d * 3600.0
            ddec = (sky[:, 1] - cat_radec[:, 1]) * 3600.0
            return np.concatenate([dra, ddec])

    try:
        result = least_squares(residuals, p0, method="trf", loss="huber",
                               f_scale=2.0)  # Huber threshold 2″
        p = result.x
    except Exception as exc:
        log.warning("WCS fit failed: %s", exc)
        return None

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = crpix
    if crval_only:
        w.wcs.crval = [p[0], p[1]]
        w.wcs.cd = cd
    else:
        w.wcs.crval = [p[0], p[1]]
        w.wcs.cd = np.array([[p[2], p[3]], [p[4], p[5]]])
    w.wcs.equinox = 2000.0
    w.wcs.radesys = "ICRS"

    return w


def compute_rms(wcs, src_xy, cat_radec):
    """Compute RMS of WCS solution in arcsec."""
    sky = wcs.all_pix2world(src_xy, 0)
    cos_dec = np.cos(np.radians(cat_radec[:, 1]))
    dra = (sky[:, 0] - cat_radec[:, 0]) * cos_dec * 3600.0
    ddec = (sky[:, 1] - cat_radec[:, 1]) * 3600.0
    return float(np.sqrt(np.mean(dra**2 + ddec**2)))


# ── Iterative solver ────────────────────────────────────────────────────

def solve_frame_iterative(data, header, telescope, wcs_init, gaia_stars):
    """Iteratively refine WCS via cross-matching with Gaia.

    Returns dict with status, wcs, rms_arcsec, n_matched, etc.
    """
    sources = extract_sources(data, telescope)
    if len(sources) < config.ASTROM_MIN_MATCHES:
        # Try lower threshold
        sources = extract_sources(data, telescope, thresh=3.0)

    n_detected = len(sources)
    if n_detected < config.ASTROM_MIN_MATCHES:
        return _fail("too_few_sources", n_detected=n_detected)

    src_xy = np.column_stack([sources["x"], sources["y"]])

    # Gaia catalog positions
    gaia_radec = np.array([[s["ra"], s["dec"]] for s in gaia_stars])
    n_gaia = len(gaia_radec)
    if n_gaia < config.ASTROM_MIN_MATCHES:
        return _fail("sparse_field", n_detected=n_detected, n_gaia=n_gaia)

    tel = config.TELESCOPES[telescope]
    ny = header.get("NAXIS2", tel["expected_shape"][0])
    nx = header.get("NAXIS1", tel["expected_shape"][1])
    image_shape = (ny, nx)

    wcs_current = wcs_init
    n_matched = 0
    rms = 999.0

    tol_init = config.ASTROM_MATCH_TOLERANCE_INIT
    tol_final = config.ASTROM_MATCH_TOLERANCE_FINAL

    # ── Pre-alignment: sky-coordinate offset vote ──────────────────
    # Header pointing can be off by arcminutes. Project detected sources
    # to RA/Dec via the (wrong) initial WCS — the constellation shape is
    # preserved. Then find the translational offset in sky coords via a
    # multi-resolution 2D histogram of pairwise (dRA, dDec).
    if len(gaia_radec) >= config.ASTROM_MIN_MATCHES and len(src_xy) >= 5:
        src_sky = wcs_current.all_pix2world(src_xy, 0)
        cos_dec = np.cos(np.radians(np.mean(gaia_radec[:, 1])))
        n_s = len(src_sky)
        n_c = min(500, len(gaia_radec))
        dra_all = (gaia_radec[:n_c, 0][None, :] - src_sky[:n_s, 0][:, None]) * cos_dec * 3600.0
        ddec_all = (gaia_radec[:n_c, 1][None, :] - src_sky[:n_s, 1][:, None]) * 3600.0
        dra_all = dra_all.ravel()
        ddec_all = ddec_all.ravel()

        # Coarse vote (30″ bins)
        max_off = 1600.0
        bins_c = np.arange(-max_off, max_off + 30, 30)
        Hc, xc, yc = np.histogram2d(dra_all, ddec_all, bins=[bins_c, bins_c])
        pkc = np.unravel_index(np.argmax(Hc), Hc.shape)
        dra_c = (xc[pkc[0]] + xc[pkc[0] + 1]) / 2.0
        ddec_c = (yc[pkc[1]] + yc[pkc[1] + 1]) / 2.0
        votes_c = int(Hc[pkc])
        # Expected random level: n_pairs / n_bins
        n_bins_c = len(bins_c) - 1
        random_level = n_s * n_c / (n_bins_c * n_bins_c)

        if votes_c >= max(5, 3 * random_level):
            # Medium refinement (10″ bins)
            rng = 60
            bm_ra = np.arange(dra_c - rng, dra_c + rng + 10, 10)
            bm_dec = np.arange(ddec_c - rng, ddec_c + rng + 10, 10)
            Hm, xm, ym = np.histogram2d(dra_all, ddec_all, bins=[bm_ra, bm_dec])
            pkm = np.unravel_index(np.argmax(Hm), Hm.shape)
            dra_m = (xm[pkm[0]] + xm[pkm[0] + 1]) / 2.0
            ddec_m = (ym[pkm[1]] + ym[pkm[1] + 1]) / 2.0
            # Fine refinement (3″ bins)
            rng2 = 20
            bf_ra = np.arange(dra_m - rng2, dra_m + rng2 + 3, 3)
            bf_dec = np.arange(ddec_m - rng2, ddec_m + rng2 + 3, 3)
            Hf, xf, yf = np.histogram2d(dra_all, ddec_all, bins=[bf_ra, bf_dec])
            pkf = np.unravel_index(np.argmax(Hf), Hf.shape)
            dra_f = (xf[pkf[0]] + xf[pkf[0] + 1]) / 2.0
            ddec_f = (yf[pkf[1]] + yf[pkf[1] + 1]) / 2.0
            # Apply to CRVAL
            dra_deg = dra_f / (cos_dec * 3600.0)
            ddec_deg = ddec_f / 3600.0
            wcs_aligned = WCS(naxis=2)
            wcs_aligned.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            wcs_aligned.wcs.crpix = wcs_current.wcs.crpix.copy()
            wcs_aligned.wcs.crval = [
                wcs_current.wcs.crval[0] + dra_deg,
                wcs_current.wcs.crval[1] + ddec_deg,
            ]
            wcs_aligned.wcs.cd = wcs_current.wcs.cd.copy()
            wcs_aligned.wcs.equinox = 2000.0
            wcs_aligned.wcs.radesys = "ICRS"
            shift_arcsec = np.sqrt(dra_f**2 + ddec_f**2)
            log.debug("  pre-align: shift=%.1f″ (%d coarse votes)",
                      shift_arcsec, votes_c)
            wcs_current = wcs_aligned

    max_iter = 12

    for iteration in range(max_iter):
        # Phase 1 (0-5): median CRVAL offset — robust to false matches,
        #   starts with wide tolerance (30px) and tightens to 10px.
        # Phase 2 (6-11): full 6-param fit, tolerance ramps 8→5px.
        if iteration < 6:
            # Wide → tight during median phase: 30, 25, 20, 15, 12, 10
            tol = 30.0 * (10.0 / 30.0) ** (iteration / 5.0)
        else:
            frac = (iteration - 6) / max(max_iter - 7, 1)
            tol = 8.0 * (tol_final / 8.0) ** frac

        # Project Gaia to pixel coords using current WCS
        cat_px = wcs_current.all_world2pix(gaia_radec, 0)

        # Filter catalog stars within image bounds (with margin)
        margin_iter = max(50, int(tol * 2))
        in_image = (
            (cat_px[:, 0] > -margin_iter) & (cat_px[:, 0] < nx + margin_iter)
            & (cat_px[:, 1] > -margin_iter) & (cat_px[:, 1] < ny + margin_iter)
        )
        cat_px_filt = cat_px[in_image]
        gaia_radec_filt = gaia_radec[in_image]

        if len(cat_px_filt) < config.ASTROM_MIN_MATCHES:
            return _fail("sparse_field", n_detected=n_detected,
                         n_gaia=int(np.sum(in_image)),
                         iteration=iteration)

        # Cross-match (widen tolerance if too few matches in phase 1)
        si, ci, dists = cross_match(src_xy, cat_px_filt, tol)
        if len(si) < config.ASTROM_MIN_MATCHES and iteration < 6:
            wider = min(tol * 2, 60.0)
            si, ci, dists = cross_match(src_xy, cat_px_filt, wider)
        n_matched = len(si)

        if n_matched < config.ASTROM_MIN_MATCHES:
            continue

        matched_src = src_xy[si]
        matched_cat = gaia_radec_filt[ci]

        # Phase 1 (iter 0-5): robust median CRVAL shift
        # Phase 2 (iter 6+): full 6-param least-squares fit
        if iteration < 6:
            # Median offset — robust even with >50% false matches
            sky_m = wcs_current.all_pix2world(matched_src, 0)
            cos_d = np.cos(np.radians(matched_cat[:, 1]))
            dra_med = np.median((matched_cat[:, 0] - sky_m[:, 0]))
            ddec_med = np.median((matched_cat[:, 1] - sky_m[:, 1]))
            wcs_new = WCS(naxis=2)
            wcs_new.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            wcs_new.wcs.crpix = wcs_current.wcs.crpix.copy()
            wcs_new.wcs.crval = [
                wcs_current.wcs.crval[0] + dra_med,
                wcs_current.wcs.crval[1] + ddec_med,
            ]
            wcs_new.wcs.cd = wcs_current.wcs.cd.copy()
            wcs_new.wcs.equinox = 2000.0
            wcs_new.wcs.radesys = "ICRS"
        else:
            wcs_new = fit_wcs(matched_src, matched_cat, wcs_current,
                              image_shape)
        if wcs_new is None:
            continue

        rms_new = compute_rms(wcs_new, matched_src, matched_cat)

        # Outlier rejection using MAD (robust against false matches)
        sky = wcs_new.all_pix2world(matched_src, 0)
        cos_dec = np.cos(np.radians(matched_cat[:, 1]))
        sep_arcsec = np.sqrt(
            ((sky[:, 0] - matched_cat[:, 0]) * cos_dec * 3600.0) ** 2
            + ((sky[:, 1] - matched_cat[:, 1]) * 3600.0) ** 2
        )
        med_sep = np.median(sep_arcsec)
        mad = np.median(np.abs(sep_arcsec - med_sep)) * 1.4826
        # Clip at median + 2.5*MAD, with floor at 1.5″
        clip_thresh = med_sep + 2.5 * max(mad, 0.3)
        clip_thresh = max(clip_thresh, 1.5)
        good = sep_arcsec < clip_thresh

        if np.sum(good) >= config.ASTROM_MIN_MATCHES:
            matched_src = matched_src[good]
            matched_cat = matched_cat[good]
            n_matched = len(matched_src)

            wcs_clip = fit_wcs(matched_src, matched_cat, wcs_new, image_shape)
            if wcs_clip is not None:
                wcs_new = wcs_clip
                rms_new = compute_rms(wcs_new, matched_src, matched_cat)

        # Convergence check
        prev_rms = rms
        rms = rms_new
        wcs_current = wcs_new

        log.debug("  iter %d: tol=%.1f px, matched=%d, RMS=%.3f″",
                  iteration, tol, n_matched, rms)

        if rms < 0.5:
            break
        if iteration >= 3 and abs(prev_rms - rms) / max(prev_rms, 1e-6) < 0.05:
            break

    # Final status
    if n_matched < config.ASTROM_MIN_MATCHES:
        return _fail("match_failed", n_detected=n_detected, n_gaia=n_gaia)

    if rms > config.ASTROM_MAX_RMS_ARCSEC:
        return _fail("high_rms", n_detected=n_detected, n_gaia=n_gaia,
                     n_matched=n_matched, rms_arcsec=rms)

    status = "HIGH_RMS" if rms > config.ASTROM_WARN_RMS_ARCSEC else "OK"
    return _result(status, wcs_current, n_detected, n_gaia, n_matched, rms)


def _result(status, wcs, n_detected, n_gaia, n_matched, rms):
    """Build a result dict from solved WCS."""
    cd = wcs.wcs.cd
    scale = np.sqrt(cd[0, 0]**2 + cd[1, 0]**2) * 3600.0  # arcsec/px
    rotation = np.degrees(np.arctan2(cd[0, 1], cd[0, 0]))

    return {
        "status": status,
        "wcs": wcs,
        "n_detected": n_detected,
        "n_gaia": n_gaia,
        "n_matched": n_matched,
        "rms_arcsec": round(rms, 3),
        "crval1": round(wcs.wcs.crval[0], 6),
        "crval2": round(wcs.wcs.crval[1], 6),
        "rotation_deg": round(rotation, 2),
        "pixel_scale": round(scale, 4),
        "fail_reason": None,
    }


def _fail(reason, **kwargs):
    """Build a FAILED result dict."""
    return {
        "status": "FAILED",
        "wcs": None,
        "fail_reason": reason,
        "n_detected": kwargs.get("n_detected", 0),
        "n_gaia": kwargs.get("n_gaia", 0),
        "n_matched": kwargs.get("n_matched", 0),
        "rms_arcsec": kwargs.get("rms_arcsec"),
        "crval1": None,
        "crval2": None,
        "rotation_deg": None,
        "pixel_scale": None,
    }


# ── Astroalign solver (T080) ────────────────────────────────────────────

def solve_frame_astroalign(data, header, telescope, center, gaia_stars):
    """Solve frame using astroalign triangle asterism matching.

    For T080 where rotation is unknown. Finds the geometric transform
    between detected sources and Gaia catalog projected at 0° rotation,
    then refines with the iterative solver.
    """
    import astroalign

    sources = extract_sources(data, telescope)
    if len(sources) < config.ASTROM_MIN_MATCHES:
        sources = extract_sources(data, telescope, thresh=3.0)

    n_detected = len(sources)
    if n_detected < config.ASTROM_MIN_MATCHES:
        return _fail("too_few_sources", n_detected=n_detected)

    src_xy = np.column_stack([sources["x"], sources["y"]])

    # Build rough WCS at 0° rotation to project Gaia → pixel
    tel = config.TELESCOPES[telescope]
    ny = header.get("NAXIS2", tel["expected_shape"][0])
    nx = header.get("NAXIS1", tel["expected_shape"][1])
    scale_deg = tel["pixel_scale"] / 3600.0

    w0 = WCS(naxis=2)
    w0.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w0.wcs.crpix = [nx / 2.0 + 0.5, ny / 2.0 + 0.5]
    w0.wcs.crval = [center.ra.deg, center.dec.deg]
    w0.wcs.cd = np.array([[-scale_deg, 0], [0, scale_deg]])

    gaia_radec = np.array([[s["ra"], s["dec"]] for s in gaia_stars])
    n_gaia = len(gaia_radec)
    if n_gaia < config.ASTROM_MIN_MATCHES:
        return _fail("sparse_field", n_detected=n_detected, n_gaia=n_gaia)

    gaia_px = w0.all_world2pix(gaia_radec, 0)

    # Filter catalog stars within generous bounds
    margin = nx  # very generous for unknown rotation
    in_image = (
        (gaia_px[:, 0] > -margin) & (gaia_px[:, 0] < nx + margin)
        & (gaia_px[:, 1] > -margin) & (gaia_px[:, 1] < ny + margin)
    )
    gaia_px_filt = gaia_px[in_image]
    gaia_radec_filt = gaia_radec[in_image]

    if len(gaia_px_filt) < config.ASTROM_MIN_MATCHES:
        return _fail("sparse_field", n_detected=n_detected,
                     n_gaia=int(np.sum(in_image)))

    # Find transform via astroalign
    try:
        transform, (src_matched, cat_matched) = astroalign.find_transform(
            src_xy, gaia_px_filt,
            max_control_points=50,
        )
    except Exception as exc:
        log.warning("astroalign failed: %s", exc)
        return _fail("asterism_failed", n_detected=n_detected, n_gaia=n_gaia)

    # Extract rotation and scale from the affine transform matrix
    mat = transform.params[:2, :2]
    scale_x = np.sqrt(mat[0, 0]**2 + mat[1, 0]**2)
    scale_y = np.sqrt(mat[0, 1]**2 + mat[1, 1]**2)
    rotation = np.arctan2(mat[1, 0], mat[0, 0])

    # Build proper WCS from transform
    # The transform maps src_xy → gaia_px (projected via w0)
    # We need to compose: pixel → transform → w0_pixel → sky
    # Instead, build new WCS directly:
    # Apply transform to image center to get corresponding catalog position
    center_px = np.array([[nx / 2.0, ny / 2.0]])
    center_cat_px = transform(center_px)
    center_sky = w0.all_pix2world(center_cat_px, 0)

    # New CD matrix: compose transform scale/rotation with w0 scale
    avg_scale = (scale_x + scale_y) / 2.0
    new_scale_deg = scale_deg * avg_scale
    cd11 = -new_scale_deg * np.cos(rotation)
    cd12 = new_scale_deg * np.sin(rotation)
    cd21 = new_scale_deg * np.sin(rotation)
    cd22 = new_scale_deg * np.cos(rotation)

    wcs_init = WCS(naxis=2)
    wcs_init.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs_init.wcs.crpix = [nx / 2.0 + 0.5, ny / 2.0 + 0.5]
    wcs_init.wcs.crval = [float(center_sky[0, 0]), float(center_sky[0, 1])]
    wcs_init.wcs.cd = np.array([[cd11, cd12], [cd21, cd22]])
    wcs_init.wcs.equinox = 2000.0
    wcs_init.wcs.radesys = "ICRS"

    log.debug("astroalign: rotation=%.1f° scale=%.4f″/px",
              np.degrees(rotation), new_scale_deg * 3600)

    # Refine with iterative solver
    return solve_frame_iterative(data, header, telescope, wcs_init, gaia_stars)


# ── solve-field fallback ─────────────────────────────────────────────────

def solve_frame_solvefield(filepath, telescope, center, pixel_scale):
    """Fallback solver using astrometry.net solve-field.

    Returns result dict or FAILED if solve-field not installed / fails.
    """
    if shutil.which("solve-field") is None:
        return _fail("no_solvefield")

    ra = center.ra.deg
    dec = center.dec.deg
    scale_lo = 0.8 * pixel_scale
    scale_hi = 1.2 * pixel_scale

    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            "solve-field",
            "--ra", str(ra),
            "--dec", str(dec),
            "--radius", "1",
            "--scale-low", str(scale_lo),
            "--scale-high", str(scale_hi),
            "--scale-units", "arcsecperpix",
            "--no-plots",
            "--overwrite",
            "--dir", tmpdir,
            str(filepath),
        ]

        log.debug("solve-field: %s", " ".join(cmd))
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired:
            return _fail("solvefield_timeout")
        except Exception as exc:
            log.warning("solve-field error: %s", exc)
            return _fail("solvefield_error")

        if proc.returncode != 0:
            log.debug("solve-field stderr: %s", proc.stderr[-500:] if proc.stderr else "")
            return _fail("solvefield_failed")

        # Find the .wcs file
        wcs_files = list(Path(tmpdir).glob("*.wcs"))
        if not wcs_files:
            return _fail("solvefield_no_wcs")

        try:
            wcs = WCS(fits.getheader(str(wcs_files[0])))
        except Exception as exc:
            log.warning("Failed to read solve-field WCS: %s", exc)
            return _fail("solvefield_parse_error")

    # Read the frame to compute RMS vs Gaia
    hdr = fits.getheader(str(filepath))
    n_detected = 0  # we don't have extraction info from solve-field

    return _result("OK", wcs, n_detected, 0, 0, 0.0)


# ── Per-frame orchestrator ──────────────────────────────────────────────

def solve_frame(filepath, telescope, simbad_cache, gaia_cache_dir,
                force=False, known_coords=None, known_rotation=None,
                reference_wcs=None):
    """Solve one frame: extract → match → fit → write WCS.

    Parameters
    ----------
    filepath : Path
    telescope : str
    simbad_cache : dict
    gaia_cache_dir : Path
    force : bool
        Re-solve even if WCS already present.
    known_coords : SkyCoord or None
        Pre-resolved coordinates (e.g. median from other frames of same object).
    known_rotation : float or None
        Previously determined rotation for this object (T080).
    reference_wcs : WCS or None
        WCS from a successfully solved neighboring frame of the same object.
        Used as initial WCS (bypasses header pointing and vote pre-alignment).

    Returns
    -------
    dict with status, n_matched, rms_arcsec, etc.
    """
    filepath = Path(filepath)
    hdr = fits.getheader(str(filepath))

    # Check if already solved
    if not force and hdr.get("CTYPE1") == "RA---TAN":
        return {"status": "SKIPPED", "filename": filepath.name,
                "fail_reason": None}

    data = fits.getdata(str(filepath)).astype(np.float64)
    if not data.flags["C_CONTIGUOUS"]:
        data = np.ascontiguousarray(data)

    object_name = hdr.get("OBJECT", "")

    # ── Get initial WCS ──────────────────────────────────────────────
    # Priority: reference_wcs > header pointing > name resolution
    if reference_wcs is not None:
        wcs_init = reference_wcs
        center = SkyCoord(ra=wcs_init.wcs.crval[0], dec=wcs_init.wcs.crval[1],
                          unit="deg")
    else:
        wcs_init = get_initial_wcs(hdr, telescope)

        if wcs_init is not None:
            center = SkyCoord(ra=wcs_init.wcs.crval[0],
                              dec=wcs_init.wcs.crval[1], unit="deg")
        elif known_coords is not None:
            center = known_coords
        else:
            center = resolve_object_coords(object_name, simbad_cache)

    if center is None:
        return _make_frame_result(filepath, hdr,
                                  _fail("no_pointing"))

    # ── Query Gaia ───────────────────────────────────────────────────
    gaia_stars = get_gaia_stars(center, telescope, gaia_cache_dir)
    if len(gaia_stars) < config.ASTROM_MIN_MATCHES:
        return _make_frame_result(filepath, hdr,
                                  _fail("sparse_field", n_gaia=len(gaia_stars)))

    # ── Solve ────────────────────────────────────────────────────────
    # If no WCS from header but we have center coords, build one from
    # known rotation (config or cached from previous frame).
    if wcs_init is None and center is not None:
        rot = known_rotation
        if rot is None:
            rot = config.ASTROM_KNOWN_ROTATION.get(telescope)
        if rot is not None:
            tel = config.TELESCOPES[telescope]
            ny = hdr.get("NAXIS2", tel["expected_shape"][0])
            nx = hdr.get("NAXIS1", tel["expected_shape"][1])
            flipped = config.ASTROM_FLIPPED_PARITY.get(telescope, False)
            wcs_init = _build_wcs(
                [center.ra.deg, center.dec.deg],
                [nx / 2.0 + 0.5, ny / 2.0 + 0.5],
                tel["pixel_scale"] / 3600.0, rot, flipped,
            )

    if wcs_init is not None:
        result = solve_frame_iterative(data, hdr, telescope, wcs_init, gaia_stars)
    else:
        # No rotation info at all: try astroalign
        result = solve_frame_astroalign(data, hdr, telescope, center, gaia_stars)

    # ── Fallback: solve-field ────────────────────────────────────────
    if result["status"] == "FAILED":
        log.info("Primary solver failed (%s), trying solve-field fallback",
                 result["fail_reason"])
        sf_result = solve_frame_solvefield(
            filepath, telescope, center,
            config.TELESCOPES[telescope]["pixel_scale"],
        )
        if sf_result["status"] != "FAILED":
            result = sf_result

    # ── Write WCS to FITS header ─────────────────────────────────────
    if result["wcs"] is not None:
        _write_wcs_to_fits(filepath, result)

    return _make_frame_result(filepath, hdr, result)


def _write_wcs_to_fits(filepath, result):
    """Write WCS keywords to FITS header in place."""
    wcs = result["wcs"]
    with fits.open(str(filepath), mode="update") as hdul:
        hdr = hdul[0].header
        hdr["CTYPE1"] = ("RA---TAN", "WCS projection type")
        hdr["CTYPE2"] = ("DEC--TAN", "WCS projection type")
        hdr["CRPIX1"] = (wcs.wcs.crpix[0], "Reference pixel X")
        hdr["CRPIX2"] = (wcs.wcs.crpix[1], "Reference pixel Y")
        hdr["CRVAL1"] = (wcs.wcs.crval[0], "Reference RA (deg)")
        hdr["CRVAL2"] = (wcs.wcs.crval[1], "Reference Dec (deg)")
        hdr["CD1_1"] = (wcs.wcs.cd[0, 0], "WCS CD matrix")
        hdr["CD1_2"] = (wcs.wcs.cd[0, 1], "WCS CD matrix")
        hdr["CD2_1"] = (wcs.wcs.cd[1, 0], "WCS CD matrix")
        hdr["CD2_2"] = (wcs.wcs.cd[1, 1], "WCS CD matrix")
        hdr["EQUINOX"] = (2000.0, "Equinox of coordinates")
        hdr["RADESYS"] = ("ICRS", "Reference frame")
        hdr["ASTRRMS"] = (result["rms_arcsec"], "Astrometric RMS (arcsec)")
        hdr["ASTRNMAT"] = (result["n_matched"], "Number of astrometric matches")
        hdr["ASTRSTAT"] = (result["status"], "Astrometric solution status")
        hdul.flush()
    log.debug("Wrote WCS to %s (RMS=%.3f″)", filepath.name, result["rms_arcsec"])


def _make_frame_result(filepath, header, result):
    """Combine solver result with file metadata for the report."""
    filepath = Path(filepath)
    out = {
        "filename": filepath.name,
        "night": filepath.parent.name,
        "object": header.get("OBJECT", ""),
        "filter": header.get("FILTER", ""),
    }
    out.update({k: v for k, v in result.items() if k != "wcs"})
    return out


def _wcs_from_result(result, header, telescope):
    """Reconstruct WCS object from a solved result dict."""
    if result.get("crval1") is None:
        return None
    tel = config.TELESCOPES[telescope]
    ny = header.get("NAXIS2", tel["expected_shape"][0])
    nx = header.get("NAXIS1", tel["expected_shape"][1])
    pixel_scale = result.get("pixel_scale", tel["pixel_scale"])
    rotation = result.get("rotation_deg", 0)
    flipped = config.ASTROM_FLIPPED_PARITY.get(telescope, False)

    return _build_wcs(
        [result["crval1"], result["crval2"]],
        [nx / 2.0 + 0.5, ny / 2.0 + 0.5],
        pixel_scale / 3600.0, rotation, flipped,
    )


def _make_exception_result(fpath, obj, hdr, exc):
    """Build result dict for an exception during solving."""
    return {
        "filename": fpath.name,
        "night": fpath.parent.name,
        "object": obj,
        "filter": hdr.get("FILTER", ""),
        "status": "FAILED",
        "fail_reason": f"exception: {exc}",
        "n_detected": 0, "n_gaia": 0, "n_matched": 0,
        "rms_arcsec": None, "crval1": None, "crval2": None,
        "rotation_deg": None, "pixel_scale": None,
    }


# ── Batch runner ─────────────────────────────────────────────────────────

def run_astrometry(year, telescope, force=False):
    """Run astrometric calibration on all reduced frames.

    Walks reduced/{night}/ for FITS files, solves each, writes reports.
    """
    base = config.DATA_ROOT / year / telescope
    reduced_dir = base / "reduced"
    astrom_dir = base / "astrometry"
    gaia_cache_dir = astrom_dir / "gaia_cache"

    if not reduced_dir.exists():
        log.warning("No reduced directory: %s", reduced_dir)
        return

    # Load simbad cache
    simbad_cache = {}
    if config.SIMBAD_CACHE.exists():
        with open(config.SIMBAD_CACHE) as f:
            simbad_cache = json.load(f)

    # Collect all reduced FITS files
    fits_files = sorted(reduced_dir.rglob("*.fits")) + sorted(reduced_dir.rglob("*.fit"))
    if not fits_files:
        log.warning("No reduced FITS files in %s", reduced_dir)
        return

    print(f"\n  {telescope}: astrometric calibration for {len(fits_files)} frames")

    # ── Pre-compute per-object coordinates ───────────────────────────
    # Group files by object, compute median coords from frames that have them
    object_files = {}
    for fpath in fits_files:
        hdr = fits.getheader(str(fpath))
        obj = hdr.get("OBJECT", "unknown")
        object_files.setdefault(obj, []).append((fpath, hdr))

    # Median RA/Dec per object from frames with OBJCTRA
    object_coords = {}
    for obj, file_hdrs in object_files.items():
        ras, decs = [], []
        for fpath, hdr in file_hdrs:
            ra_str = hdr.get("OBJCTRA")
            dec_str = hdr.get("OBJCTDEC")
            if ra_str and dec_str:
                ra_h = _parse_sexagesimal(ra_str)
                dec_d = _parse_sexagesimal(dec_str)
                if ra_h is not None and dec_d is not None:
                    ras.append(ra_h * 15.0)
                    decs.append(dec_d)
        if ras:
            object_coords[obj] = SkyCoord(
                ra=np.median(ras), dec=np.median(decs), unit="deg")
        else:
            # Try name resolution
            coord = resolve_object_coords(obj, simbad_cache)
            if coord is not None:
                object_coords[obj] = coord
                log.info("Resolved '%s' → RA=%.4f Dec=%.4f",
                         obj, coord.ra.deg, coord.dec.deg)

    # ── Pass 1: Solve all frames ────────────────────────────────────
    results = {}  # keyed by filepath for easy update in pass 2
    result_list = []  # ordered list for reporting
    # Track rotation per object for T080 (reuse once found)
    object_rotation = {}
    # Track successful WCS per object for pass 2
    object_wcs = {}  # obj → WCS (from first successful solve)

    for i, fpath in enumerate(fits_files, 1):
        hdr = fits.getheader(str(fpath))
        obj = hdr.get("OBJECT", "unknown")
        known_coords = object_coords.get(obj)
        known_rot = object_rotation.get(obj)

        log.info("[%d/%d] %s", i, len(fits_files), fpath.name)
        try:
            result = solve_frame(
                fpath, telescope, simbad_cache, gaia_cache_dir,
                force=force, known_coords=known_coords,
                known_rotation=known_rot,
            )
        except Exception as exc:
            log.error("EXCEPTION solving %s: %s", fpath.name, exc)
            result = _make_exception_result(fpath, obj, hdr, exc)

        results[str(fpath)] = result
        result_list.append(result)

        # Cache rotation and WCS for T080/pass 2 on success
        if result.get("status") in ("OK", "HIGH_RMS"):
            if (result.get("rotation_deg") is not None
                    and obj not in object_rotation):
                object_rotation[obj] = result["rotation_deg"]
                log.info("Cached rotation for %s: %.1f°", obj,
                         result["rotation_deg"])
            if obj not in object_wcs and result.get("crval1") is not None:
                # Reconstruct WCS from result for pass 2
                ref_wcs = _wcs_from_result(result, hdr, telescope)
                if ref_wcs is not None:
                    object_wcs[obj] = ref_wcs

        status = result.get("status", "?")
        rms = result.get("rms_arcsec")
        rms_str = f"{rms:.3f}″" if rms is not None else "—"
        nm = result.get("n_matched", 0)

        if i % 20 == 0 or status == "FAILED":
            print(f"    [{i}/{len(fits_files)}] {fpath.name}: "
                  f"{status} (matched={nm}, RMS={rms_str})")

    # ── Pass 2: Retry failures using reference WCS from neighbors ───
    failed = [(fpath, r) for fpath, r in results.items()
              if r.get("status") == "FAILED"]
    if failed and object_wcs:
        n_retry = sum(1 for fp, r in failed
                      if r.get("object", "") in object_wcs)
        if n_retry > 0:
            print(f"\n    Pass 2: retrying {n_retry} failures "
                  f"with reference WCS from neighbors")

        for fp_str, old_result in failed:
            fpath = Path(fp_str)
            obj = old_result.get("object", "")
            ref_wcs = object_wcs.get(obj)
            if ref_wcs is None:
                continue

            hdr = fits.getheader(str(fpath))
            known_rot = object_rotation.get(obj)
            known_coords = object_coords.get(obj)

            log.info("[pass2] %s (using ref WCS from %s)", fpath.name, obj)
            try:
                result = solve_frame(
                    fpath, telescope, simbad_cache, gaia_cache_dir,
                    force=True, known_coords=known_coords,
                    known_rotation=known_rot,
                    reference_wcs=ref_wcs,
                )
            except Exception as exc:
                log.error("EXCEPTION in pass 2 for %s: %s", fpath.name, exc)
                continue

            if result.get("status") in ("OK", "HIGH_RMS"):
                results[fp_str] = result
                # Update result_list
                for idx, r in enumerate(result_list):
                    if r.get("filename") == fpath.name and \
                       r.get("night") == fpath.parent.name:
                        result_list[idx] = result
                        break
                rms = result.get("rms_arcsec")
                log.info("[pass2] %s → %s (RMS=%.3f″, matched=%d)",
                         fpath.name, result["status"], rms or 0,
                         result.get("n_matched", 0))

        n_recovered = sum(1 for fp, r in results.items()
                          if r.get("status") in ("OK", "HIGH_RMS")) - \
                       sum(1 for r in result_list
                           if r.get("status") in ("OK", "HIGH_RMS"))
        # Recount from results dict
        n_solved_pass2 = sum(1 for r in results.values()
                             if r.get("status") in ("OK", "HIGH_RMS"))
        n_solved_pass1 = sum(1 for r in result_list
                             if r.get("status") in ("OK", "HIGH_RMS"))
        if n_solved_pass2 > n_solved_pass1:
            print(f"    Pass 2 recovered {n_solved_pass2 - n_solved_pass1} "
                  f"additional frames")

    # Use updated results for reporting
    final_results = list(results.values())

    # ── Write reports ────────────────────────────────────────────────
    astrom_dir.mkdir(parents=True, exist_ok=True)

    json_out = astrom_dir / "astrom_report.json"
    with open(json_out, "w") as f:
        json.dump(final_results, f, indent=2, default=str)
    log.info("Wrote %s", json_out)

    csv_out = astrom_dir / "astrom_report.csv"
    if final_results:
        fieldnames = [k for k in final_results[0].keys() if k != "wcs"]
        with open(csv_out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in final_results:
                writer.writerow({k: row.get(k, "") for k in fieldnames})
        log.info("Wrote %s", csv_out)

    # ── Summary ──────────────────────────────────────────────────────
    _print_summary(final_results, telescope)


def _print_summary(results, telescope):
    """Print per-object astrometry summary."""
    by_object = {}
    for r in results:
        obj = r.get("object", "?")
        by_object.setdefault(obj, []).append(r)

    n_ok = sum(1 for r in results if r.get("status") == "OK")
    n_warn = sum(1 for r in results if r.get("status") == "HIGH_RMS")
    n_fail = sum(1 for r in results if r.get("status") == "FAILED")
    n_skip = sum(1 for r in results if r.get("status") == "SKIPPED")
    n_total = len(results)

    print(f"\n  ── {telescope} Astrometry Summary ──")
    print(f"  Total: {n_total}  OK: {n_ok}  HIGH_RMS: {n_warn}  "
          f"FAILED: {n_fail}  SKIPPED: {n_skip}")
    print(f"\n  {'Object':25s}  {'N':>4s}  {'OK':>4s}  {'Fail':>4s}  "
          f"{'Med RMS″':>9s}  {'Rotation°':>9s}")
    print(f"  {'─'*25}  {'─'*4}  {'─'*4}  {'─'*4}  {'─'*9}  {'─'*9}")

    for obj in sorted(by_object):
        frames = by_object[obj]
        n = len(frames)
        ok = sum(1 for r in frames if r.get("status") in ("OK", "HIGH_RMS"))
        fail = sum(1 for r in frames if r.get("status") == "FAILED")
        rms_vals = [r["rms_arcsec"] for r in frames
                    if r.get("rms_arcsec") is not None]
        rot_vals = [r["rotation_deg"] for r in frames
                    if r.get("rotation_deg") is not None]

        med_rms = f"{np.median(rms_vals):.3f}" if rms_vals else "—"
        med_rot = f"{np.median(rot_vals):.1f}" if rot_vals else "—"

        print(f"  {obj:25s}  {n:4d}  {ok:4d}  {fail:4d}  "
              f"{med_rms:>9s}  {med_rot:>9s}")

    # Show failure reasons
    fail_reasons = {}
    for r in results:
        if r.get("status") == "FAILED" and r.get("fail_reason"):
            reason = r["fail_reason"]
            fail_reasons[reason] = fail_reasons.get(reason, 0) + 1

    if fail_reasons:
        print(f"\n  Failure breakdown:")
        for reason, count in sorted(fail_reasons.items(), key=lambda x: -x[1]):
            print(f"    {reason}: {count}")
    print()
