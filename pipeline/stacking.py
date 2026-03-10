"""Image stacking: SExtractor → SCAMP → SWarp pipeline.

SCAMP handles astrometric alignment against Gaia, SWarp does resampling
and sigma-clipped coaddition.
"""

import json
import logging
import os
import shutil
import statistics
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

from . import config

log = logging.getLogger(__name__)

SWARP_CMD = config.SWARP_CMD
SCAMP_CMD = "scamp"

# ── Target name normalisation ───────────────────────────────────────────
# Map variant OBJECT names to a canonical target identifier.
# Add entries here as new aliases are discovered.
TARGET_ALIASES = {
    # Coma cluster fields
    "ACO 1656": "Abell1656",
    "ACO1656": "Abell1656",
    "ACO1656V2": "Abell1656",
    "COMA-cluster": "Abell1656",
    "COMA-cluster-border": "Abell1656_border",
    "Coma Cote": "Abell1656_border",
    "coma": "Abell1656",
    "coma cote": "Abell1656_border",
    # Pandora's cluster (Abell 2744)
    "Pandora": "Abell2744",
    "pandora": "Abell2744",
    "Pandora 8": "Abell2744",
    "pandora 05": "Abell2744",
    "pandora 4": "Abell2744",
    "pandora 6": "Abell2744",
    "pandora 7": "Abell2744",
    "Pandora r 2": "Abell2744",
    "Pandora r 3": "Abell2744",
    # M67 = NGC 2682
    "M67": "M67",
    "M67 40s": "M67",
    "NGC 2682": "M67",
    "NGC2682": "M67",
    "NGC 2682 30sec": "M67",
    "NGC 2682 5sec": "M67",
    # M48 = NGC 2548
    "M48": "M48",
    "NGC 2548": "M48",
    "NGC2548": "M48",
    "NGC 2548 30sec": "M48",
    "NGC 2548 5sec": "M48",
    # NGC/Messier with extra text
    "N628 R": "NGC628",
    # Vulcano asteroid variants
    "Vulcano": "Vulcano",
    "Vullcano": "Vulcano",
    "vulcano": "Vulcano",
    # 29P variants
    "29P": "29P",
    "29P 2": "29P",
    # 62P variants
    "62P": "62P",
    "62P 2": "62P",
    "62P joli": "62P",
    # Euphrosyne variants
    "Euphrosyne": "Euphrosyne",
    "Euphrosyne 1": "Euphrosyne",
    "Euphrosyne 2": "Euphrosyne",
    "Euphrosyne 3": "Euphrosyne",
    "Euphrosyne 4": "Euphrosyne",
    "Euphrosyne 5": "Euphrosyne",
    # J0254
    "J0254": "J0254",
    "J0254 B": "J0254",
    "J0254 R": "J0254",
    # Abell 1703
    "A1703": "Abell1703",
    # aco1228
    "aco1228": "Abell1228",
    # NGC variants with prefix/suffix
    "T080 NGC2683": "NGC2683",
    "T080_NGC2683": "NGC2683",
    "N628 R": "NGC628",
    # NGC 2548 = M48 variants
    "NGC 2548 30sec": "M48",
    "NGC 2548 5sec": "M48",
    # NGC 2682 = M67 variants
    "NGC 2682 30sec": "M67",
    "NGC 2682 5sec": "M67",
    # NGC with underscored spaces
    "NGC 2301": "NGC2301",
    "NGC 2420": "NGC2420",
    "NGC 2539": "NGC2539",
    # M104 with trailing underscore
    "M104": "M104",
    "Coma N4874": "NGC4874",
    "Coma_N4874": "NGC4874",
    # M99 variants
    "M99 2": "M99",
    "M99_2": "M99",
    # NGC with spaces
    "NGC 3718": "NGC3718",
    # Tourbillon = M51
    "Tourbillon": "M51",
    # M97 = Owl nebula
    "Nebuleuse de la chouette": "M97",
    # M64 = Black Eye
    "Black Eye": "M64",
}

# Targets to skip (asteroids, comets, transits, calibration, planets, unresolvable)
SKIP_PATTERNS = [
    # Transits (single-band time series)
    "HAT-P", "Hat-P", "Hat P", "Hat_P", "WASP", "Wasp", "TrES-", "TOI-",
    # Calibration / test
    "Champ", "champ", "CHAMP", "Calibration", "reference", "bias", "flat",
    "A1 80s", "A2 200s",
    # Solar system (moving targets — can't stack)
    "moon", "Moon", "lune", "Lune", "jupiter", "Jupiter",
    # Unresolvable student names
    "McCracken", "McCracker", "saMar30", "SCAX", "X81876", "Kirik",
    "2023_EO", "2023 EO", "2023 E0", "C2023", "C2019",
    "ceres", "Cheval", "a10tqqp", "EuclidQ1",
    "C0QQC", "C11ZDU", "C2AJT", "A10k", "ZTs",
    "PGZ1", "PSZ1",
    "Autosave Image", "Autosave_Image",
    "A1 80s", "A2 200s", "A1_80s", "A2_200s",
]

# Asteroids / comets — skip for deep stacking (moving targets)
ASTEROID_PATTERNS = [
    "62P", "29P", "144P", "12P", "118p",
    "CK10", "2726", "2003 CK",
    "Euphrosyne", "Frigga", "Kozai", "Comete", "comete",
    "Vulcano", "Vullcano", "vulcano",
]


def normalise_target(obj_name):
    """Return canonical target name, or None if target should be skipped."""
    # Strip trailing underscores/spaces
    name = obj_name.strip().rstrip("_").strip()
    # Also check with spaces restored (scanner replaces spaces with _)
    name_sp = name.replace("_", " ")

    for check in (name, name_sp):
        if any(p.lower() in check.lower() for p in SKIP_PATTERNS):
            return None
        if any(p.lower() in check.lower() for p in ASTEROID_PATTERNS):
            return None

    # Try alias lookup with original, stripped, and space-restored forms
    for key in (obj_name, name, name_sp):
        if key in TARGET_ALIASES:
            return TARGET_ALIASES[key]
    return name


def _normalise_filter(filt):
    """Normalise filter names so R/r and similar variants group together."""
    mapping = {
        "R": "R", "r": "R", "r'": "R",
        "R Cousins": "R", "R_Cousins": "R", "R_cousins": "R",
        "r_Gunn": "R", "r Gunn": "R",
        "G": "g", "g": "g", "g'": "g",
        "g_Gunn": "g", "g Gunn": "g",
        "v_Gunn": "V", "v Gunn": "V",
        "i_Gunn": "i", "i Gunn": "i",
        "B": "B", "B Cousins": "B", "B_Cousins": "B", "B_cousins": "B",
        "V": "V", "V Cousins": "V", "V_Cousins": "V", "V_cousins": "V",
        "I": "I", "i": "i", "i'": "i",
        "H alpha OHP": "Ha", "H_alpha": "Ha", "Halpha": "Ha", "Ha": "Ha",
        "Halpha_OHP": "Ha",
    }
    return mapping.get(filt, filt)


# ── Inventory ───────────────────────────────────────────────────────────

def build_target_inventory(years=None, telescopes=None):
    """Scan all reduced frames and group by (canonical_target, telescope, filter).

    Returns dict: canonical_target -> list of frame dicts, each with:
        path, year, telescope, filter, object, exptime, has_wcs, photzp
    """
    if years is None:
        years = range(2018, 2026)
    if telescopes is None:
        telescopes = list(config.TELESCOPES)

    inventory = defaultdict(list)

    for year in years:
        for tel in telescopes:
            reduced_dir = config.DATA_ROOT / str(year) / tel / "reduced"
            if not reduced_dir.is_dir():
                continue
            fits_files = sorted(reduced_dir.rglob("*.fits")) + \
                         sorted(reduced_dir.rglob("*.fit"))
            for fpath in fits_files:
                try:
                    hdr = fits.getheader(str(fpath))
                except Exception:
                    continue
                obj = hdr.get("OBJECT", "")
                canonical = normalise_target(obj)
                if canonical is None:
                    continue

                filt = hdr.get("FILTER", "?")
                has_wcs = "CRVAL1" in hdr and "CD1_1" in hdr
                photzp = hdr.get("PHOTZP")
                exptime = hdr.get("EXPTIME", 0)
                fwhm = hdr.get("QCFWHM")
                astrstat = hdr.get("ASTRSTAT", "")
                astrrms = hdr.get("ASTRRMS")

                inventory[canonical].append({
                    "path": str(fpath),
                    "year": year,
                    "telescope": tel,
                    "filter": _normalise_filter(filt),
                    "filter_raw": filt,
                    "object": obj,
                    "exptime": exptime,
                    "has_wcs": has_wcs,
                    "astrstat": astrstat,
                    "astrrms": astrrms,
                    "photzp": photzp,
                    "fwhm": fwhm,
                })

    return dict(inventory)


def print_inventory(inventory):
    """Print a summary table of the target inventory."""
    print(f"\n{'Target':25s}  {'Tel':4s}  {'Filter':8s}  {'Years':15s}  "
          f"{'N':>4s}  {'WCS':>4s}  {'ZP':>4s}")
    print("-" * 75)
    for target in sorted(inventory):
        frames = inventory[target]
        # Group by (telescope, filter)
        groups = defaultdict(list)
        for fr in frames:
            groups[(fr["telescope"], fr["filter"])].append(fr)
        for (tel, filt), group in sorted(groups.items()):
            years = sorted(set(f["year"] for f in group))
            n_wcs = sum(1 for f in group if f["has_wcs"])
            n_zp = sum(1 for f in group if f["photzp"] is not None)
            year_str = ",".join(str(y) for y in years)
            print(f"{target:25s}  {tel}   {filt:8s}  {year_str:15s}  "
                  f"{len(group):4d}  {n_wcs:4d}  {n_zp:4d}")


# ── Stacking ────────────────────────────────────────────────────────────

def _select_reference(frames):
    """Pick the best frame as astrometric reference.

    Prefers: has WCS, lowest FWHM, highest ZP (most transparent).
    """
    wcs_frames = [f for f in frames if f["has_wcs"]]
    if not wcs_frames:
        return None

    def score(f):
        # Lower is better: prioritise good seeing, then transparency
        fwhm = f["fwhm"] if f["fwhm"] is not None else 99.0
        # Negate ZP so higher ZP = lower score
        zp = -(f["photzp"] or 0)
        return (fwhm, zp)

    return min(wcs_frames, key=score)


def _compute_flxscale(frames, zp_ref):
    """Compute FLXSCALE for each frame relative to reference ZP.

    Frames without ZP get FLXSCALE=1.0 (no correction).
    """
    for f in frames:
        if f["photzp"] is not None and zp_ref is not None:
            f["flxscale"] = 10 ** (0.4 * (zp_ref - f["photzp"]))
        else:
            f["flxscale"] = 1.0


def _reject_outliers(frames, sigma=2.0):
    """Reject frames with ZP far from the median (likely clouded)."""
    zps = [f["photzp"] for f in frames if f["photzp"] is not None]
    if len(zps) < 3:
        return frames  # not enough to clip

    med = statistics.median(zps)
    mad = statistics.median(abs(z - med) for z in zps)
    if mad == 0:
        return frames

    # 1.4826 converts MAD to sigma-equivalent
    sigma_est = mad * 1.4826
    keep = []
    for f in frames:
        if f["photzp"] is None:
            keep.append(f)  # keep frames without ZP (can't judge)
        elif abs(f["photzp"] - med) <= sigma * sigma_est:
            keep.append(f)
        else:
            log.info("  Rejecting %s: ZP=%.2f (median=%.2f, %.1fσ)",
                     Path(f["path"]).name, f["photzp"], med,
                     abs(f["photzp"] - med) / sigma_est)
    return keep


def _run_sextractor_for_scamp(filepath, telescope, ldac_path):
    """Run SExtractor on a frame, producing LDAC catalogue for SCAMP."""
    tel = config.TELESCOPES[telescope]
    pixel_scale = tel["pixel_scale"]
    saturation = tel["saturation"]

    tmpdir = Path(ldac_path).parent

    param_path = tmpdir / "scamp.param"
    param_path.write_text(
        "XWIN_IMAGE\nYWIN_IMAGE\n"
        "ERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\n"
        "FLUX_AUTO\nFLUXERR_AUTO\n"
        "MAG_AUTO\nMAGERR_AUTO\n"
        "FLAGS\n"
        "FLUX_RADIUS\nELONGATION\nSNR_WIN\n"
    )

    sex_path = tmpdir / "scamp.sex"
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

    cmd = [config.SEXTRACTOR_CMD, str(filepath), "-c", str(sex_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        log.warning("    SExtractor failed on %s: %s",
                    Path(filepath).name, proc.stderr[-200:] if proc.stderr else "")
        return False
    return Path(ldac_path).exists()


def _run_scamp(ldac_paths, telescope, workdir):
    """Run SCAMP on a set of LDAC catalogues.

    Writes .head files alongside each LDAC catalogue.
    Returns True on success.
    """
    tel = config.TELESCOPES[telescope]
    pixel_scale = tel["pixel_scale"]

    scamp_args = [
        SCAMP_CMD,
        *[str(p) for p in ldac_paths],
        "-ASTREF_CATALOG", "GAIA-DR3",
        "-DISTORT_DEGREES", "1",
        "-STABILITY_TYPE", "INSTRUMENT",
        "-CROSSID_RADIUS", "3.0",
        "-POSITION_MAXERR", "2.0",        # arcmin — pipeline WCS is close
        "-POSANGLE_MAXERR", "2.0",        # degrees
        "-PIXSCALE_MAXERR", "1.01",
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

    log.info("    Running SCAMP on %d catalogues", len(ldac_paths))
    try:
        proc = subprocess.run(scamp_args, capture_output=True, text=True,
                              timeout=600, cwd=str(workdir))
        # SCAMP writes progress to stderr
        if proc.stderr:
            # Log contrast info
            for line in proc.stderr.splitlines():
                if "contrast" in line.lower() or "warning" in line.lower():
                    log.info("    SCAMP: %s", line.strip())
        if proc.returncode != 0:
            log.error("    SCAMP failed: %s",
                      proc.stderr[-500:] if proc.stderr else "no stderr")
            return False
    except subprocess.TimeoutExpired:
        log.error("    SCAMP timed out")
        return False

    # Check that .head files were created
    n_heads = sum(1 for p in ldac_paths
                  if Path(str(p).replace(".ldac", ".head")).exists())
    log.info("    SCAMP produced %d/%d .head files", n_heads, len(ldac_paths))
    return n_heads > 0


def _reject_spatial_outliers(frames, max_offset_deg=0.5):
    """Reject frames whose CRVAL is far from the median position.

    Catches mis-identified targets (e.g. Simbad resolving to wrong object).
    """
    from astropy.io import fits as afits

    ras, decs = [], []
    for f in frames:
        h = afits.getheader(f["path"])
        ras.append(h.get("CRVAL1", 0))
        decs.append(h.get("CRVAL2", 0))

    if not ras:
        return frames

    med_ra = statistics.median(ras)
    med_dec = statistics.median(decs)
    cos_dec = np.cos(np.radians(med_dec))

    keep = []
    for f, ra, dec in zip(frames, ras, decs):
        offset = np.sqrt(((ra - med_ra) * cos_dec) ** 2 + (dec - med_dec) ** 2)
        if offset <= max_offset_deg:
            keep.append(f)
        else:
            log.info("  Rejecting spatial outlier %s: (%.3f, %.3f) "
                     "offset=%.1f° from median",
                     Path(f["path"]).name, ra, dec, offset)
    return keep


def stack_group(frames, output_path, reference=None, reject_sigma=2.0,
                combine_type="CLIPPED"):
    """Stack frames using SExtractor → SCAMP → SWarp pipeline.

    SCAMP handles astrometric alignment against Gaia DR3.
    SWarp does resampling and sigma-clipped coaddition.

    Parameters
    ----------
    frames : list of dict
        Frame dicts from build_target_inventory.
    output_path : Path
        Output directory for coadded images.
    reference : dict or None
        Reference frame dict. If None, auto-selected.
    reject_sigma : float
        Sigma-clip threshold for ZP outlier rejection.
    combine_type : str
        SWarp COMBINE_TYPE (CLIPPED, WEIGHTED, MEDIAN, AVERAGE).
    """
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # SCAMP needs a rough WCS in the headers to know where to query Gaia.
    # Filter to frames that have at least a CRVAL (even if astrometry was poor).
    all_frames = [f for f in frames if f["has_wcs"]]

    # Reject spatial outliers (frames with bogus positions)
    all_frames = _reject_spatial_outliers(all_frames)

    # Reject frames with poor astrometric RMS (> 0.5")
    MAX_ASTRRMS = 0.5
    n_before = len(all_frames)
    all_frames = [f for f in all_frames
                  if f.get("astrrms") is not None and f["astrrms"] <= MAX_ASTRRMS]
    n_rejected = n_before - len(all_frames)
    if n_rejected:
        log.info("  Rejected %d frames with ASTRRMS > %.1f\" (or missing RMS)",
                 n_rejected, MAX_ASTRRMS)
    if len(all_frames) < 2:
        log.warning("  Need at least 2 frames with WCS to stack")
        return {}

    telescope = all_frames[0]["telescope"]

    # Pick reference for flux scaling
    if reference is None:
        reference = _select_reference(all_frames)
    ref_name = Path(reference["path"]).name if reference else "none"
    log.info("  Reference: %s (FWHM=%.1f\", ZP=%s)",
             ref_name,
             (reference["fwhm"] or 0) if reference else 0,
             f"{reference['photzp']:.2f}" if reference and reference["photzp"] else "—")

    # Group by filter
    by_filter = defaultdict(list)
    for f in all_frames:
        by_filter[f["filter"]].append(f)

    results = {}
    for filt in sorted(by_filter):
        filt_frames = by_filter[filt]
        log.info("  Filter %s: %d frames", filt, len(filt_frames))

        if len(filt_frames) < 2:
            log.info("    Skipping (need >= 2 frames)")
            continue

        # Reject ZP outliers
        filt_frames = _reject_outliers(filt_frames, sigma=reject_sigma)
        if not filt_frames:
            log.warning("    All frames rejected for %s", filt)
            continue

        # Compute flux scaling
        zps = [f["photzp"] for f in filt_frames if f["photzp"] is not None]
        if zps:
            zp_ref = max(zps)
            _compute_flxscale(filt_frames, zp_ref)
            has_fluxscale = True
        else:
            zp_ref = None
            has_fluxscale = False
            for f in filt_frames:
                f["flxscale"] = 1.0
            log.warning("    No ZPs — stacking without flux scaling")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Step 1: Copy frames (primary HDU only) and run SExtractor
            image_paths = []
            ldac_paths = []
            frame_map = {}  # image_path -> frame dict

            for f in filt_frames:
                src = Path(f["path"])
                stem = f"{f['year']}_{src.stem}"
                dst = tmpdir / f"{stem}.fits"
                ldac = tmpdir / f"{stem}.ldac"

                with fits.open(str(src)) as hdul:
                    h = hdul[0].header.copy()
                    h["FLXSCALE"] = (f["flxscale"],
                                     "Flux scale factor for coaddition")
                    if "GAIN" not in h:
                        h["GAIN"] = (1.0, "Effective gain (e-/ADU)")
                    # Standardise keywords so SCAMP sees one instrument
                    h["FILTER"] = filt
                    h["INSTRUME"] = "Andor"
                    fits.writeto(str(dst), hdul[0].data, h, overwrite=True)

                if _run_sextractor_for_scamp(str(dst), telescope, str(ldac)):
                    image_paths.append(dst)
                    ldac_paths.append(ldac)
                    frame_map[str(dst)] = f
                else:
                    log.warning("    Skipping %s (SExtractor failed)", src.name)

            if len(ldac_paths) < 2:
                log.warning("    Too few catalogues for SCAMP")
                results[filt] = {"status": "SEX_FAILED",
                                 "n_input": len(filt_frames)}
                continue

            # Step 2: SCAMP astrometric alignment
            scamp_ok = _run_scamp(ldac_paths, telescope, tmpdir)
            if not scamp_ok:
                results[filt] = {"status": "SCAMP_FAILED",
                                 "n_input": len(ldac_paths)}
                continue

            # Step 3: SWarp coaddition
            # SCAMP writes .head files named after the LDAC catalogues.
            # SWarp reads .head files matching the image filename.
            # Since images are {stem}.fits and LDACs are {stem}.ldac,
            # SCAMP's {stem}.head already matches the image.
            swarp_inputs = []
            for img, ldac in zip(image_paths, ldac_paths):
                head_file = img.with_suffix(".head")
                if head_file.exists():
                    swarp_inputs.append(str(img))
                else:
                    log.warning("    No .head for %s — skipping", img.name)

            if not swarp_inputs:
                results[filt] = {"status": "NO_HEADS",
                                 "n_input": len(ldac_paths)}
                continue

            output_file = output_path / f"coadd_{filt}.fits"
            weight_file = output_path / f"coadd_{filt}.weight.fits"

            swarp_args = [
                SWARP_CMD,
                *swarp_inputs,
                "-IMAGEOUT_NAME", str(output_file),
                "-WEIGHTOUT_NAME", str(weight_file),
                "-COMBINE", "Y",
                "-COMBINE_TYPE", combine_type,
                "-CLIP_SIGMA", "3.0",
                "-RESAMPLE", "Y",
                "-RESAMPLE_DIR", str(tmpdir),
                "-SUBTRACT_BACK", "N",
                "-BACK_SIZE", "128",
                "-FSCALE_KEYWORD", "FLXSCALE",
                "-FSCALASTRO_TYPE", "NONE",
                "-CELESTIAL_TYPE", "NATIVE",
                "-PROJECTION_TYPE", "TAN",
                "-CENTER_TYPE", "ALL",
                "-PIXEL_SCALE", "0",
                "-PIXELSCALE_TYPE", "MEDIAN",
                "-RESAMPLING_TYPE", "LANCZOS3",
                "-VMEM_DIR", str(tmpdir),
                "-DELETE_TMPFILES", "Y",
                "-VERBOSE_TYPE", "QUIET",
                "-WRITE_XML", "N",
            ]

            log.info("    Running SWarp on %d frames → %s",
                     len(swarp_inputs), output_file.name)
            try:
                proc = subprocess.run(swarp_args, capture_output=True,
                                      text=True, timeout=300)
                if proc.returncode != 0:
                    log.error("    SWarp failed: %s",
                              proc.stderr[-500:] if proc.stderr else "")
                    results[filt] = {"status": "SWARP_FAILED",
                                     "n_input": len(swarp_inputs)}
                    continue
            except subprocess.TimeoutExpired:
                log.error("    SWarp timed out")
                results[filt] = {"status": "TIMEOUT",
                                 "n_input": len(swarp_inputs)}
                continue

        # Update coadd header with provenance
        if output_file.exists():
            with fits.open(str(output_file), mode="update") as hdul:
                h = hdul[0].header
                h["NCOMBINE"] = (len(swarp_inputs), "Number of input frames")
                h["FILTER"] = (filt, "Filter band")
                if zp_ref is not None:
                    h["PHOTZP"] = (round(zp_ref, 4),
                                   "Zero point of coadd (brightest input)")
                h["COMBTYPE"] = (combine_type, "Combine method")
                h["FLXSCALD"] = (has_fluxscale,
                                 "Flux scaling applied before coaddition")
                for i, f in enumerate(filt_frames[:999]):
                    h[f"INP{i+1:04d}"] = (Path(f["path"]).name,
                                           f"Input frame {i+1}")
                hdul.flush()

            if zp_ref:
                log.info("    OK: %s (%d frames, ZP=%.2f)",
                         output_file.name, len(swarp_inputs), zp_ref)
            else:
                log.info("    OK: %s (%d frames, no ZP)",
                         output_file.name, len(swarp_inputs))

            results[filt] = {
                "status": "OK",
                "output": str(output_file),
                "weight": str(weight_file),
                "n_input": len(swarp_inputs),
                "zp_ref": zp_ref,
                "filters": filt,
            }
        else:
            results[filt] = {"status": "NO_OUTPUT",
                             "n_input": len(swarp_inputs)}

    return results


# ── Batch entry point ───────────────────────────────────────────────────

def run_stacking(target=None, years=None, telescope=None, force=False):
    """Stack all (or selected) targets.

    Parameters
    ----------
    target : str or None
        Stack only this target (canonical name). If None, stack all.
    years : list of int or None
        Years to include. Default: all.
    telescope : str or None
        Restrict to one telescope.
    force : bool
        Re-stack even if output exists.
    """
    telescopes = [telescope] if telescope else list(config.TELESCOPES)
    inventory = build_target_inventory(years=years, telescopes=telescopes)

    if target:
        # Filter to requested target
        if target in inventory:
            inventory = {target: inventory[target]}
        else:
            # Try case-insensitive match
            matches = {k: v for k, v in inventory.items()
                       if k.lower() == target.lower()}
            if matches:
                inventory = matches
            else:
                log.error("Target '%s' not found in inventory", target)
                print_inventory(build_target_inventory(years=years,
                                                       telescopes=telescopes))
                return

    stack_dir = config.DATA_ROOT / "stacks"
    all_results = {}

    for canonical in sorted(inventory):
        frames = inventory[canonical]
        # Group by telescope — don't mix pixel scales
        by_tel = defaultdict(list)
        for f in frames:
            by_tel[f["telescope"]].append(f)

        for tel in sorted(by_tel):
            tel_frames = by_tel[tel]
            n_wcs = sum(1 for f in tel_frames if f["has_wcs"])
            if n_wcs < 2:
                continue

            out_dir = stack_dir / canonical / tel
            if not force and out_dir.exists() and list(out_dir.glob("coadd_*.fits")):
                log.info("Skipping %s/%s (output exists)", canonical, tel)
                continue

            years_present = sorted(set(f["year"] for f in tel_frames))
            filters_present = sorted(set(f["filter"] for f in tel_frames))
            log.info("Stacking %s / %s: %d frames (%d with WCS), "
                     "filters=%s, years=%s",
                     canonical, tel, len(tel_frames), n_wcs,
                     ",".join(filters_present),
                     ",".join(str(y) for y in years_present))

            results = stack_group(tel_frames, out_dir)
            all_results[f"{canonical}/{tel}"] = results

    # Write summary report
    if all_results:
        report_path = stack_dir / "stack_report.json"
        with open(report_path, "w") as f:
            json.dump(all_results, f, indent=2)
        log.info("Wrote %s", report_path)

    # Print summary
    print(f"\n{'Target/Tel':35s}  {'Filter':8s}  {'N':>4s}  {'ZP':>7s}  {'Status':10s}")
    print("-" * 72)
    for key in sorted(all_results):
        for filt, res in sorted(all_results[key].items()):
            zp_str = f"{res['zp_ref']:.2f}" if res.get("zp_ref") else "—"
            print(f"{key:35s}  {filt:8s}  {res['n_input']:4d}  "
                  f"{zp_str:>7s}  {res['status']:10s}")
