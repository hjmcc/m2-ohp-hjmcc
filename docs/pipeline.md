# OHP M2 Reduction Pipeline

Data reduction pipeline for the OHP (Observatoire de Haute-Provence) M2 student observing programme. Handles calibration, quality control, astrometry, photometric calibration, and PSF extraction for imaging data from the T080 (80cm) and T120 (1.2m) telescopes.

## Telescopes

| Telescope | Aperture | CCD | Pixel scale | FOV | Binning |
|-----------|----------|-----|-------------|-----|---------|
| T080 | 80cm Cassegrain | SBIG STXL-6303 | 0.484"/px | 12.4' × 8.3' | 2×2 |
| T120 | 1.2m Newton | Andor 1k×1k | 0.773"/px | 13.2' × 13.2' | 2×2 |

## Installation

### Dependencies

Python ≥ 3.10, plus:

```
astropy
astroquery
scipy
sep
photutils
```

For PSF extraction, the following external tools must be installed:
- **SExtractor** ≥ 2.28 (available as `sex` on the PATH)
- **PSFex** ≥ 3.24

### Data layout

The pipeline expects data under `$OHP_DATA_ROOT` (default: `/Users/hjmcc/Work/data`) with the structure:

```
{DATA_ROOT}/
  {year}/
    {telescope}/
      reduced/                  # calibrated science frames (input)
        {night}/
          *.fits / *.fit
      qc/                       # quality control outputs
        frame_stats.json
        frame_stats.csv
      astrometry/               # astrometric solutions
        astrom_report.json
        gaia_cache/
      photometry/               # photometric calibration outputs
        phot_report.json
        phot_report.csv
        ps1_cache/
      psf/                      # PSF models and stats
        psf_report.json
        psf_report.csv
        mean_psf.fits
        {night}/
          {frame_stem}.psf
```

## Usage

```bash
# Run individual stages
python -m pipeline qc         --year 2025 --telescope T080
python -m pipeline astrometry --year 2025 --telescope T080
python -m pipeline photometry --year 2025 --telescope T080
python -m pipeline psf        --year 2025 --telescope T080

# Run the full pipeline (calibrate → qc → astrometry → photometry → psf)
python -m pipeline all --year 2025 --telescope T080

# View pipeline status summary
python -m pipeline status --year 2025

# Re-process frames that already have results
python -m pipeline photometry --year 2025 --telescope T120 --force
```

If `--telescope` is omitted, all configured telescopes are processed.

## Pipeline stages

### 1. Calibration (`pipeline/calibration.py`)

Bias subtraction and flat-field division. Groups calibration frames by night, constructs master bias and master flat, applies to science frames. Output: `reduced/{night}/*.fits`.

### 2. Quality Control (`pipeline/quality.py`)

Per-frame source extraction with `sep` to measure:
- Source count, FWHM (median of good sources), elongation
- Background level and RMS
- QC flags: `OK`, `LOW_SOURCES`, `BAD_SEEING`, `ELONGATED`, `SATURATED`

Results written to FITS headers and `qc/frame_stats.json`.

### 3. Astrometry (`pipeline/astrometry.py`)

WCS calibration via Gaia DR3 cross-matching:
- Two-pass solver: pass 1 solves all frames; pass 2 retries failures using WCS from solved neighbours as initial guess
- Multi-resolution sky-coordinate offset voting (30" → 10" → 3" bins) for pre-alignment
- Progressive match tolerance: 30 → 10px (median offset), then 8 → 5px (6-parameter affine fit)
- Writes CD matrix, `CRVAL1/2`, `ASTRSTAT`, `ASTRRMS` to FITS headers
- Results: `astrometry/astrom_report.json`

### 4. Photometric Calibration (`pipeline/photometry.py`)

Per-frame zero-point via aperture photometry cross-matched to Pan-STARRS1 DR2:

1. Source extraction with `sep` — aperture radius = 3 × median FWHM
2. WCS pixel → sky coordinate transformation
3. PS1 DR2 query via VizieR (`II/349/ps1`) with local JSON file caching
4. Cross-match at 5" tolerance using `scipy.spatial.cKDTree`
5. Instrumental magnitude: `m_inst = -2.5 × log10(flux / exptime)`
6. Zero-point: iterative sigma-clipped median of `(m_cat − m_inst)`
7. Writes `PHOTZP`, `PHOTZPER`, `PHOTNSTR`, `PHOTCAT`, `PHOTSTAT` to FITS header

Only processes frames with WCS solutions (`ASTRSTAT` = `OK` or `HIGH_RMS`).

**Filter → PS1 band mapping:**

| Telescope | Filter | PS1 band | Notes |
|-----------|--------|----------|-------|
| T080 | g' | gmag | |
| T080 | r' | rmag | |
| T080 | i' | imag | |
| T080 | H_alpha | — | Skipped (no broadband equivalent) |
| T120 | R | rmag | |
| T120 | V | gmag | ZP absorbs colour term |
| T120 | B | gmag | ZP absorbs colour term |
| T120 | G | gmag | |

### 5. PSF Extraction (`pipeline/psf.py`)

Per-frame PSF models via SExtractor → PSFex:

1. SExtractor run with `CATALOG_TYPE FITS_LDAC`, extracting 35×35 pixel vignettes
2. PSFex modelling with `PIXEL_AUTO` basis, 25×25 pixel PSF size, degree-2 spatial variation
3. Statistics parsed from PSFex VOTABLE XML output
4. Per-frame `.psf` model files saved; intermediate LDAC catalogues cleaned up
5. Results: `psf/psf_report.json` and `psf/psf_report.csv`

A **mean PSF** per telescope is also available as `psf/mean_psf.fits` — the median stack of all centre-of-field PSF components, normalised to unit peak.

---

## 2025 Results

### Photometric Zero-Points

#### T080 (80cm Cassegrain)

97 of 114 frames calibrated (85%). 14 H_alpha frames skipped (no PS1 equivalent), 3 r' frames with too few matches.

| Filter | N frames | PS1 band | ZP (median) | ZP (mean) | ZP std | Median matches |
|--------|----------|----------|-------------|-----------|--------|----------------|
| g' | 38 | gmag | 20.64 | 20.72 | 0.97 | 17 |
| r' | 35 | rmag | 21.35 | 21.34 | 1.19 | 17 |
| i' | 24 | imag | 21.20 | 21.22 | 1.25 | 29 |

#### T120 (1.2m Newton)

213 of 359 frames calibrated (59%). Many crowded open-cluster fields have too few clean matches.

| Filter | N frames | PS1 band | ZP (median) | ZP (mean) | ZP std | Median matches |
|--------|----------|----------|-------------|-----------|--------|----------------|
| B | 37 | gmag | 25.34 | 25.55 | 1.66 | 26 |
| V | 40 | gmag | 25.91 | 25.95 | 1.35 | 30 |
| R | 112 | rmag | 25.00 | 25.00 | 1.49 | 7 |
| G | 24 | gmag | 24.68 | 24.82 | 1.87 | 6 |

**Notes on ZP scatter:**
- The 1–2 mag standard deviations are dominated by non-photometric conditions (frequent cloud cover at OHP) and crowded-field aperture contamination.
- In sparse fields observed under stable conditions, the frame-to-frame ZP scatter is 0.02–0.05 mag (e.g. `vulcano`: σ = 0.017, `M60`: σ = 0.020).
- Crowded open clusters (NGC 1912, NGC 2281) show large ZP scatter (> 1 mag) due to source blending in aperture photometry.
- The ZP values are approximate and should be used for relative comparisons within a night rather than absolute photometry.

### PSF Quality

#### T080

99 of 114 frames with successful PSF models (87%).

| Filter | N frames | FWHM (median) | FWHM (px) | Ellipticity | χ² |
|--------|----------|---------------|-----------|-------------|-----|
| g' | 32 | 3.81" | 7.9 | 0.049 | 1.31 |
| r' | 29 | 3.76" | 7.8 | 0.049 | 1.36 |
| i' | 24 | 3.74" | 7.7 | 0.053 | 1.43 |
| H_alpha | 14 | 3.35" | 6.9 | 0.147 | 1.63 |

#### T120

359 of 359 frames with successful PSF models (100%).

| Filter | N frames | FWHM (median) | FWHM (px) | Ellipticity | χ² |
|--------|----------|---------------|-----------|-------------|-----|
| B | 50 | 3.50" | 4.5 | 0.053 | 1.16 |
| V | 45 | 3.04" | 3.9 | 0.059 | 1.33 |
| R | 220 | 3.14" | 4.1 | 0.070 | 1.14 |
| G | 44 | 3.38" | 4.4 | 0.046 | 1.03 |

### Mean PSF

Median-stacked centre-of-field PSF models per telescope (normalised to unit peak):

| Telescope | N frames stacked | FWHM (px) | FWHM (arcsec) | EE50 radius | EE80 radius |
|-----------|-----------------|-----------|---------------|-------------|-------------|
| T080 | 91 | 4.10 | 1.98" | 1.08" | 2.00" |
| T120 | 359 | 4.36 | 3.37" | 1.73" | 2.79" |

The mean PSF FITS files (`psf/mean_psf.fits`) are 25×25 pixel arrays with headers containing `FWHM_PX`, `FWHM_AS`, `EE50_AS`, `EE80_AS`, `PIXSCALE`, and `NFRAMES`.

Note: The T080 per-frame PSFex FWHM values (~3.8") are systematically larger than the mean-PSF FWHM (1.98") because PSFex reports the FWHM of its spatially-varying model evaluated across the full field, while the mean PSF is the centre-of-field component only.

## Pipeline modules

| Module | Lines | Description |
|--------|-------|-------------|
| `config.py` | 99 | Paths, telescope parameters, QC and algorithm thresholds |
| `scanner.py` | 293 | Archive scanner — reads FITS headers, classifies frame types |
| `calibration.py` | 297 | Master bias/flat construction, science frame calibration |
| `quality.py` | 441 | Source extraction QC — FWHM, elongation, background stats |
| `astrometry.py` | 1317 | WCS calibration via Gaia DR3 cross-matching |
| `photometry.py` | 492 | Photometric zero-points via PS1 DR2 cross-matching |
| `psf.py` | 456 | PSF extraction via SExtractor + PSFex |
| `utils.py` | 91 | Shared utilities — logging, JSON I/O, path helpers |
| `__main__.py` | 330 | CLI entry point and command dispatch |
