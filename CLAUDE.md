# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Astronomical imaging data archive and reduction pipeline for the OHP (Observatoire de Haute-Provence) M2 student observing programme. Three main components:

- **`pipeline/`** — Python package for CCD reduction, astrometry, photometry, and PSF extraction
- **`stack_all.py`** / **`check_phot.py`** — Batch SWarp stacking with shared WCS grids and photometric validation
- **`ohp-archive/`** — Standalone archive scanner (pure stdlib, no astropy) that indexes 8 years of FITS data and generates HTML dashboards

Full pipeline documentation: `docs/pipeline_description.md`

GitHub Pages site served from `docs/` on `main` branch: https://hjmcc.github.io/m2-ohp-hjmcc/

## Repository structure

```
pipeline/              # Python package — per-frame reduction pipeline
  config.py            # Single source for telescope params, thresholds, paths
  scanner.py           # FITS classification (bias/flat/science)
  calibration.py       # Master bias/flat, science frame reduction
  quality.py           # FWHM, elongation, QC flags via sep
  astrometry.py        # WCS via solve-field + Gaia DR3 validation
  photometry.py        # Aperture photometry + PS1 DR2 zero points
  psf.py               # SExtractor + PSFex PSF modelling
  stacking.py          # Target inventory, SCAMP-based stacking (legacy)
  __main__.py          # CLI entry point
stack_all.py           # Batch SWarp stacking — shared WCS grids, ZP=30
check_phot.py          # Photometric validation of stacks vs PS1
ohp-archive/           # Pure-stdlib archive scanner for telescope use
scripts/               # One-off utility scripts
  process_all.sh       # Batch reprocessing shell script
  make_phot_plots.py   # Photometry diagnostic plots
  resolve_t120.py      # T120 coordinate resolution
  run_scamp.py         # Manual SCAMP runs
docs/                  # Documentation and GitHub Pages
  pipeline_description.md  # Complete pipeline description
calibrations/          # Reference PSF models (T080, T120)
notebooks/             # Jupyter analysis notebooks
data/                  # Pipeline results summaries
```

## Commands

```bash
# Activate Python environment
micromamba activate astropy-3.12

# Run individual pipeline stages
python -m pipeline scan         --year 2025 --telescope T120
python -m pipeline calibrate    --year 2025 --telescope T120
python -m pipeline qc           --year 2025 --telescope T120
python -m pipeline astrometry   --year 2025 --telescope T120
python -m pipeline photometry   --year 2025 --telescope T120
python -m pipeline psf          --year 2025 --telescope T120

# Run full pipeline (calibrate → qc → astrometry → photometry → psf)
python -m pipeline all --year 2025 --telescope T120

# Status overview
python -m pipeline status --year 2025

# Re-process already-completed frames
python -m pipeline photometry --year 2025 --telescope T120 --force

# Batch stacking (all targets, shared WCS grids, ZP=30)
python stack_all.py
python stack_all.py --target M67
python stack_all.py --dry-run

# Validate stack photometry against PS1
python check_phot.py --plot
python check_phot.py --fix         # apply corrections

# Archive scanner (standalone, no astropy needed)
cd ohp-archive && python3 scan_archive.py
```

Omitting `--telescope` processes all configured telescopes. There is no test suite.

## Dependencies

Pipeline requires: `astropy`, `astroquery`, `scipy`, `sep`, `ccdproc`, `numpy`. Use `micromamba activate astropy-3.12`.

External tools: **solve-field** (on PATH), **SExtractor** (`sex` on PATH), **PSFex** (`/Users/hjmcc/opt/astromatic/bin/psfex`), **SWarp** (`/Users/hjmcc/opt/astromatic/bin/swarp`), **SCAMP** (`scamp` on PATH).

`$OHP_DATA_ROOT` sets the data directory (default: `/Users/hjmcc/Work/data`).

## Architecture

### Pipeline stages (sequential, each reads previous stage's output)

1. **scanner.py** — Reads FITS headers, classifies frame types (bias/flat/science), resolves target names via Simbad
2. **calibration.py** — Builds master bias/flat per night, applies to science frames → `reduced/`
3. **quality.py** — `sep` source extraction for FWHM (via `sep.flux_radius` half-light radius), elongation, background stats → QC flags in headers + `qc/frame_stats.json`
4. **astrometry.py** — WCS via `solve-field` plate solver + Gaia DR3 cross-match validation → CD matrix + `ASTRSTAT` in headers
5. **photometry.py** — Aperture photometry (`r = 3 × FWHM`) + PS1 DR2 cross-match → `PHOTZP` in headers (counts/sec convention: `m = -2.5 log10(flux/exptime) + ZP`)
6. **psf.py** — SExtractor LDAC → PSFex for per-frame PSF models with degree-2 spatial variation

### Stacking (stack_all.py)

- Processes all T120 targets (2018–2025) with SWarp
- Shared WCS grid per target group: all filters pixel-aligned for multi-band photometry
- Merged groups: Coma (Abell1656+border+NGC4874), M105_group, M38_group, Markarian (M84+M86)
- Flux scaling to ZP=30: `FLXSCALE = 10^(0.4*(30 - ZP)) / exptime` written into `.head` files
- Weight maps from normalised master flat-fields
- MEDIAN combine, LANCZOS3 resampling, background subtraction
- Output: `$OHP_DATA_ROOT/stacks/{group}/T120/{group}_{filter}.fits`

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
  master_bias/                # master bias frame
  master_flat/                # master flats per filter (normalised to median=1)
  reduced/{night}/*.fits      # calibrated science (input to QC+)
  qc/frame_stats.json         # QC results
  astrometry/                 # WCS solutions + gaia_cache/
  photometry/                 # ZPs + ps1_cache/
  psf/                        # .psf models

$OHP_DATA_ROOT/stacks/
  {group}/T120/
    {group}_{filter}.fits         # coadded image (ZP=30)
    {group}_{filter}.weight.fits  # weight map
    plots/                        # diagnostic photometry plots
  ps1_cache/                      # PS1 cache for stack validation
  phot_check_report.json          # photometric validation report
```

### Archive scanner (`ohp-archive/scan_archive.py`)

Intentionally **pure stdlib** (no astropy) so it can run at the telescope. Reads FITS headers with a minimal parser, handles T152 spectroscopy quirks (non-standard headers, classify from filename). Outputs self-contained HTML dashboards with embedded CSS/JS/data. `OUTPUT_DIR` is derived from `__file__` path.
