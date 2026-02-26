# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Astronomical imaging data archive and reduction pipeline for the OHP (Observatoire de Haute-Provence) M2 student observing programme. Two main components:

- **`pipeline/`** — Python package for CCD reduction, astrometry, photometry, and PSF extraction
- **`ohp-archive/`** — Standalone archive scanner (pure stdlib, no astropy) that indexes 8 years of FITS data and generates HTML dashboards

GitHub Pages site served from `docs/` on `main` branch: https://hjmcc.github.io/m2-ohp-hjmcc/

## Commands

```bash
# Run individual pipeline stages
python -m pipeline scan         --year 2025 --telescope T080
python -m pipeline calibrate    --year 2025 --telescope T080
python -m pipeline qc           --year 2025 --telescope T080
python -m pipeline astrometry   --year 2025 --telescope T080
python -m pipeline photometry   --year 2025 --telescope T080
python -m pipeline psf          --year 2025 --telescope T080

# Run full pipeline (calibrate → qc → astrometry → photometry → psf)
python -m pipeline all --year 2025 --telescope T080

# Status overview
python -m pipeline status --year 2025

# Re-process already-completed frames
python -m pipeline photometry --year 2025 --telescope T120 --force

# Archive scanner (standalone, no astropy needed)
cd ohp-archive && python3 scan_archive.py
```

Omitting `--telescope` processes all configured telescopes. There is no test suite.

## Dependencies

Pipeline requires: `astropy`, `astroquery`, `scipy`, `sep`, `ccdproc`, `numpy`. Install manually (no requirements.txt).

PSF extraction requires external tools: **SExtractor** (`sex` on PATH) and **PSFex** (hardcoded at `/Users/hjmcc/opt/astromatic/bin/psfex`).

`$OHP_DATA_ROOT` sets the data directory (default: `/Users/hjmcc/Work/data`).

## Architecture

### Pipeline stages (sequential, each reads previous stage's output)

1. **scanner.py** — Reads FITS headers, classifies frame types (bias/flat/science), resolves target names via Simbad
2. **calibration.py** — Builds master bias/flat per night, applies to science frames → `reduced/`
3. **quality.py** — `sep` source extraction for FWHM, elongation, background stats → QC flags in headers + `qc/frame_stats.json`
4. **astrometry.py** — WCS via Gaia DR3 cross-matching with two-pass solver and multi-resolution offset voting → CD matrix + `ASTRSTAT` in headers
5. **photometry.py** — Aperture photometry + PS1 DR2 cross-match for per-frame zero-points → `PHOTZP` in headers
6. **psf.py** — SExtractor LDAC → PSFex for per-frame PSF models with degree-2 spatial variation

### Key design patterns

- **config.py is the single source** for all telescope parameters, thresholds, and paths. `config.TELESCOPES` dict keys (`T120`, `T080`) are used everywhere.
- Each stage has a `run_*(year, telescope, force)` batch entry point and writes both FITS header keywords and a `*_report.json` + `.csv`.
- External catalogue queries (Gaia, PS1, Simbad) are cached to local JSON files to survive re-runs. Pattern: check cache file → query if missing → write cache.
- `__main__.py` uses lazy imports (`from .photometry import run_photometry` inside command functions) to avoid loading unused stages.

### Telescope-specific quirks

- **T080** has flipped parity (Cassegrain mirror, `det(CD) > 0`) — handled via `config.ASTROM_FLIPPED_PARITY`. Pixel scale is 0.484"/px (empirically determined, not manufacturer spec). Inherently elongated PSF — `ELONGATED` QC flags on T080 are expected.
- **T120** has standard parity, 0.773"/px, rotation ~-90°.
- Many flats have `IMAGETYP='Light Frame'` (labeling error) — detected via OBJECT keyword containing "flat" tokens.

### Data layout

```
$OHP_DATA_ROOT/{year}/{telescope}/
  reduced/{night}/*.fits    # calibrated science (input to QC+)
  qc/frame_stats.json       # QC results
  astrometry/               # WCS solutions + gaia_cache/
  photometry/               # ZPs + ps1_cache/
  psf/                      # .psf models + mean_psf.fits
```

### Archive scanner (`ohp-archive/scan_archive.py`)

Intentionally **pure stdlib** (no astropy) so it can run at the telescope. Reads FITS headers with a minimal parser, handles T152 spectroscopy quirks (non-standard headers, classify from filename). Outputs self-contained HTML dashboards with embedded CSS/JS/data. `OUTPUT_DIR` is derived from `__file__` path.
