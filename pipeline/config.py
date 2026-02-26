"""Pipeline configuration: paths, telescope parameters, QC thresholds."""

import os
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
ARCHIVE_JSON = PROJECT_ROOT / "ohp-archive" / "ohp_archive.json"
SIMBAD_CACHE = PROJECT_ROOT / "ohp-archive" / "simbad_cache.json"
DATA_ROOT = Path(os.environ.get("OHP_DATA_ROOT", "/Users/hjmcc/Work/data"))

# ── Telescope parameters ──────────────────────────────────────────────────
TELESCOPES = {
    "T080": {
        "name": "80cm Cassegrain",
        "instrument": "SBIG STXL-6303",
        "pixel_scale": 0.484,       # arcsec/px (2x2 binned) — empirically determined
        "expected_shape": (1023, 1536),  # (NAXIS2, NAXIS1)
        "binning": "2x2",
        "nominal_temp": -25.0,       # °C — science CCD temperature
        "temp_tolerance": 5.0,       # °C — flag if bias temp differs by more
        "saturation": 65535,         # ADU (16-bit unsigned)
    },
    "T120": {
        "name": "1.2m Newton",
        "instrument": "Andor 1k×1k",
        "pixel_scale": 0.773,       # arcsec/px (2x2 binned)
        "expected_shape": (1024, 1024),
        "binning": "2x2",
        "nominal_temp": -60.0,
        "temp_tolerance": 5.0,
        "saturation": 65535,
    },
}

# ── QC thresholds ─────────────────────────────────────────────────────────
BIAS_SIGMA_CLIP = 3.0       # reject bias frames > N sigma from group median
FLAT_MIN_ADU = 10_000       # reject flats below this median signal
FLAT_MAX_ADU = 55_000       # reject flats above this median signal
FLAT_GRADIENT_MAX = 1.15    # reject if quadrant max/min ratio exceeds this
FLAT_SIGMA_CLIP = 3.0       # reject flat frames > N sigma from group median

# ── SEP source-extraction thresholds ─────────────────────────────────────
SEP_THRESH = 3.0            # detection threshold in σ above background
SEP_MINAREA = 5             # minimum connected pixels for a detection
FWHM_MIN_PX = 2.0           # minimum FWHM for "good" source (pixels)
FWHM_MAX_PX = 15.0          # maximum FWHM for "good" source (pixels)

# ── Science frame QC flagging thresholds ────────────────────────────────
QC_MIN_GOOD_SOURCES = 5     # fewer → LOW_SOURCES
QC_MAX_FWHM_ARCSEC = 6.0    # above → BAD_SEEING
QC_MAX_ELONGATION = 1.3      # above → ELONGATED (tracking problem)
QC_MAX_SATURATED_FRAC = 0.01 # above → SATURATED

# ── Astrometry ───────────────────────────────────────────────────────────
ASTROM_SEP_THRESH = 5.0             # σ above background (higher than QC)
ASTROM_MAX_SOURCES = 200            # brightest sources to use
ASTROM_MATCH_TOLERANCE_INIT = 200.0 # initial vote search range (pixels)
ASTROM_MATCH_TOLERANCE_FINAL = 5.0  # final match tolerance (pixels)
ASTROM_MAX_ITER = 5
ASTROM_SIGMA_CLIP = 3.0
ASTROM_MIN_MATCHES = 6              # minimum for valid 6-param fit
ASTROM_MAX_RMS_ARCSEC = 10.0        # above → FAILED
ASTROM_WARN_RMS_ARCSEC = 3.0        # above → HIGH_RMS warning
ASTROM_GAIA_MAG_LIMIT = {"T120": 18.0, "T080": 18.0}
ASTROM_GAIA_MAG_DEEP = 20.0         # fallback if < 20 stars at default
ASTROM_GAIA_QUERY_MARGIN = 2.5      # multiply FOV radius for cone search (generous for pointing errors)
ASTROM_KNOWN_ROTATION = {"T120": -90.0, "T080": 136.0}  # degrees
ASTROM_FLIPPED_PARITY = {"T120": False, "T080": True}  # True = det(CD) > 0

# ── Photometry ──────────────────────────────────────────────────────────
PHOT_APERTURE_MULT = 3.0         # aperture radius = mult × median FWHM (px)
PHOT_ANNULUS_INNER = 4.0         # sky annulus inner = mult × FWHM
PHOT_ANNULUS_OUTER = 7.0         # sky annulus outer = mult × FWHM
PHOT_ZP_SIGMA_CLIP = 2.5
PHOT_ZP_MIN_STARS = 5
PHOT_MATCH_RADIUS_ARCSEC = 5.0
PHOT_PS1_MAG_FAINT = 20.0
PHOT_PS1_MAG_BRIGHT = 13.0
PHOT_SEP_THRESH = 5.0
PHOT_FILTER_TO_PS1 = {
    "g'": "gmag", "g": "gmag", "G": "gmag",
    "r'": "rmag", "r": "rmag", "R": "rmag",
    "i'": "imag", "i": "imag",
    "V": "gmag", "B": "gmag",
}

# ── PSF (SExtractor + PSFex) ───────────────────────────────────────────
SEXTRACTOR_CMD = "sex"
PSFEX_CMD = "/Users/hjmcc/opt/astromatic/bin/psfex"
PSF_VIGNET_SIZE = 35
PSF_PSFVAR_DEGREES = 2
PSF_DETECT_THRESH = 5.0
PSF_MIN_SOURCES = 20
PSF_SAMPLE_MINSN = 20
PSF_SAMPLE_MAXELLIP = 0.3

# ── Pipeline metadata ────────────────────────────────────────────────────
PIPELINE_NAME = "OHP-M2-v1"
