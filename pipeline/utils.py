"""Utility functions: archive loading, FITS I/O, logging, directory setup."""

import json
import logging
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.nddata import CCDData
import astropy.units as u

from . import config

log = logging.getLogger("pipeline")


def setup_logging(verbose: bool = False):
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S"
    ))
    root = logging.getLogger("pipeline")
    root.setLevel(level)
    if not root.handlers:
        root.addHandler(handler)


# ── Archive JSON ──────────────────────────────────────────────────────────

def load_archive():
    """Load the archive JSON produced by scan_archive.py."""
    with open(config.ARCHIVE_JSON) as f:
        return json.load(f)


def get_cal_frames(archive, year: str, telescope: str, binning: str = "2x2"):
    """Return cal_index entry: {bias: [...], dark: [...], flat: {filter: [...]}}."""
    return archive["cal_index"].get(year, {}).get(telescope, {}).get(binning, {
        "bias": [], "dark": [], "flat": {}
    })


def get_science_frames(archive, year: str, telescope: str):
    """Return list of science file entries for given year/telescope."""
    frames = []
    for obj_name, files in archive["file_index"].items():
        for f in files:
            if f["date"].startswith(year) and f["telescope"] == telescope:
                frames.append(f)
    return frames


# ── FITS I/O ──────────────────────────────────────────────────────────────

def read_ccd(path: str, unit="adu") -> CCDData:
    """Read a FITS file as CCDData, handling BZERO/BSCALE unsigned-16 format."""
    hdu = fits.open(path, do_not_scale_image_hdu=False)
    data = hdu[0].data.astype(np.float32)
    header = hdu[0].header
    hdu.close()
    return CCDData(data, unit=u.Unit(unit), header=header)


def write_ccd(ccd: CCDData, path: Path, overwrite: bool = True):
    """Write CCDData to FITS as float32."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ccd.data = ccd.data.astype(np.float32)
    ccd.write(str(path), overwrite=overwrite)
    log.debug("Wrote %s", path)


# ── Output directories ───────────────────────────────────────────────────

def output_dir(year: str, telescope: str, subdir: str = "") -> Path:
    """Return output directory path, creating it if needed."""
    d = config.DATA_ROOT / year / telescope
    if subdir:
        d = d / subdir
    d.mkdir(parents=True, exist_ok=True)
    return d


def write_json(data, path: Path):
    """Write dict to JSON file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    log.debug("Wrote %s", path)
