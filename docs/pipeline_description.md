# OHP M2 Imaging Pipeline: Complete Description

This document describes the end-to-end reduction pipeline for the OHP (Observatoire de Haute-Provence) M2 student observing programme. The pipeline processes raw CCD frames from two telescopes through calibration, quality control, astrometric and photometric calibration, PSF modelling, and image stacking to produce science-ready multi-band coadded images with a common zero point of 30.0 mag.

## 1. Data and Instrumentation

### 1.1 Telescopes

| Parameter | T120 (1.2 m Newton) | T080 (80 cm Cassegrain) |
|-----------|---------------------|-------------------------|
| CCD | Andor 1k &times; 1k | SBIG STXL-6303 |
| Binning | 2&times;2 | 2&times;2 |
| Image size | 1024 &times; 1024 px | 1023 &times; 1536 px |
| Pixel scale | 0.773 &Prime;/px | 0.484 &Prime;/px |
| FOV | ~13.2&prime; &times; 13.2&prime; | ~8.3&prime; &times; 12.4&prime; |
| Nominal CCD temp | &minus;60 &deg;C | &minus;25 &deg;C |
| Saturation | 65 535 ADU | 65 535 ADU |
| Image rotation | &minus;90&deg; | +136&deg; |
| Parity | Standard (det(CD) < 0) | Flipped (det(CD) > 0) |

The T080 has flipped parity due to the Cassegrain mirror and an inherently elongated PSF; `ELONGATED` QC flags on T080 frames are expected and not indicative of tracking errors.

### 1.2 Data layout

Raw data are organised by observing run under `$OHP_DATA_ROOT` (default `/Users/hjmcc/Work/data`):

```
$OHP_DATA_ROOT/
  {year}/
    {telescope}/
      {night}/                  # raw frames (YYYYMMDD)
      master_bias/              # master bias (output of calibration)
      master_flat/              # master flats per filter (output)
      reduced/
        {night}/*.fits          # calibrated science frames
      qc/
        frame_stats.json        # per-frame quality metrics
      astrometry/
        gaia_cache/             # cached Gaia DR3 queries (JSON)
        astrom_report.json      # per-frame astrometry results
      photometry/
        ps1_cache/              # cached PS1 DR2 queries (JSON)
        phot_report.json        # per-frame photometry results
      psf/
        *.psf                   # PSFex models per frame
        psf_report.json         # per-frame PSF statistics
  stacks/
    {target}/
      T120/
        {target}_{filter}.fits          # coadded image
        {target}_{filter}.weight.fits   # weight map
        plots/                          # diagnostic photometry plots
    ps1_cache/                  # PS1 cache for stack validation
    phot_check_report.json      # photometric validation report
```

The archive spans 8 years (2018&ndash;2025) with ~2900 reduced T120 science frames.

### 1.3 Filters

Raw filter names (~50 variants) are normalised to 33 canonical filters. The broadband filters used for photometric calibration against Pan-STARRS1 are:

| Canonical filter | PS1 band | Raw variants |
|------------------|----------|--------------|
| g | gmag | g, G, g_Gunn, g Gunn, g' |
| R | rmag | R, r, R Cousins, R_Cousins, r_Gunn, r Gunn, r' |
| i | imag | i, I, i_Gunn, i Gunn, i' |
| V | gmag | V, V Cousins, V_Cousins, v_Gunn, v Gunn |
| B | gmag | B, B Cousins, B_Cousins |

Narrowband filters (H&alpha;, [O III], [S II]) and U-band have no PS1 counterpart and are stacked without photometric validation.


## 2. Pipeline Stages

The pipeline is invoked via `python -m pipeline {stage}` with `--year` and optional `--telescope` arguments. Stages are sequential; each reads the output of the preceding stage. The `all` command runs the full sequence: `calibrate` &rarr; `qc` &rarr; `astrometry` &rarr; `photometry` &rarr; `psf`.

### 2.1 Stage 0: Scanning (`pipeline/scanner.py`)

**Command:** `python -m pipeline scan --year 2025 --telescope T120`

The scanner walks `$OHP_DATA_ROOT/{year}/{telescope}/` recursively, reads FITS headers, and classifies each frame as bias, dark, flat, or science. It builds an archive JSON used by the calibration stage.

**Frame classification logic:**

1. **Flat detection override**: many OHP flats are mislabelled with `IMAGETYP='Light Frame'`. The scanner detects these by checking whether the filename starts with `flatdome`, `skyflat`, `flat_`, etc., or whether the `OBJECT` keyword contains the token `flat` (case-insensitive, whitespace/underscore-delimited).
2. **IMAGETYP parsing**: standard FITS keywords (`Flat Field`, `Bias Frame`, `Dark Frame`, `Light Frame`).
3. **Filename fallback**: if `IMAGETYP` is missing (common on T152 spectroscopy frames), the scanner classifies from filename prefixes (`bias`, `biais`, `offset` &rarr; bias; `flat`, `ff` &rarr; flat; `dark` &rarr; dark).

**Filter normalisation** maps ~50 raw filter strings to canonical names via a lookup table (`FILTER_MAP`) with case-insensitive fallback.

**Output:** `{year}/{telescope}/archive.json` containing:
- `cal_index`: bias, dark, and flat frames indexed by binning and filter.
- `file_index`: science frames indexed by object name.


### 2.2 Stage 1: Calibration (`pipeline/calibration.py`)

**Command:** `python -m pipeline calibrate --year 2025 --telescope T120`

Produces master bias frames, master flat-field frames, and bias-subtracted, flat-corrected science frames.

#### 2.2.1 Master bias

All bias frames for a given year/telescope/binning combination are combined:

1. **Input QC** (`quality.check_bias_frames`):
   - Reject frames with `EXPTIME > 0` (test exposures accidentally labelled as bias).
   - Flag CCD temperature mismatches: reject if `|CCD_TEMP - nominal_temp| > 5.0 °C`.
   - Compute median ADU level for each frame; sigma-clip outliers at 3.0&sigma; from the group median.
2. **Combination**: `ccdproc.combine(method='average', sigma_clip=True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=3)`.
3. **Output**: `master_bias/master_bias_{year}_{telescope}_{binning}.fits` with header keywords `NCOMBINE`, `COMBMETH`, `BIASQUAL`.

#### 2.2.2 Master flats

Flat frames are combined per filter:

1. **Input QC** (`quality.check_flat_frames`):
   - **Signal level**: reject frames with median ADU below 10 000 or above 55 000.
   - **Gradient**: divide the image into four quadrants, compute the ratio of the highest to lowest quadrant median; reject if ratio > 1.15 (non-uniform illumination beyond tolerance).
   - **Sigma clipping**: reject frames > 3.0&sigma; from the per-filter group median.
2. **Bias subtraction**: subtract the master bias from each accepted flat (if a master bias exists).
3. **Combination**: `ccdproc.combine(method='median', scale=1/median, sigma_clip=True, 3σ)` &mdash; each flat is scaled to its own median before combining.
4. **Normalisation**: the combined flat is divided by its median to give a flat field with median exactly 1.0. Pixel values directly represent relative sensitivity: vignetting and dust shadows produce values < 1.
5. **Output**: `master_flat/master_flat_{year}_{telescope}_{binning}_{filter}.fits`.

#### 2.2.3 Science frame reduction

Each science frame is processed individually:

1. **Quick QC rejects**: frames with `EXPTIME = 0`, `OBJECT` containing calibration keywords (`bias`, `dark`, `flat`, `test`), or shape mismatch with the telescope's expected dimensions.
2. **Bias subtraction**: `ccdproc.subtract_bias(ccd, master_bias)` if a master bias is available. Sets `BIASCORR` header keyword.
3. **Flat-field correction**: if the frame's filter matches an available master flat, `ccdproc.flat_correct(ccd, master_flat)`. Sets `FLATCORR` to the flat's filter name. If no flat is available, `FLATCORR = 'NONE'` and the frame is bias-only corrected.
4. **Output**: `reduced/{night}/{filename}` preserving the night subdirectory structure.

**Report:** `qc/qc_report.json` summarising bias quality, flat acceptance counts per filter, and science frame reduction statistics.


### 2.3 Stage 2: Quality Control (`pipeline/quality.py`)

**Command:** `python -m pipeline qc --year 2025 --telescope T120`

Measures image quality metrics for every reduced science frame using the `sep` source extraction library (a C implementation of SExtractor algorithms).

#### Per-frame measurements (`measure_frame_stats`):

1. **Image statistics**: median, MAD-based sigma (`MAD × 1.4826`), min/max, saturation fraction (fraction of pixels > 90% of the saturation level).
2. **Source extraction**: `sep.extract()` at 3.0&sigma; above the locally-estimated background, with a minimum connected area of 5 pixels.
3. **FWHM measurement**: for each detected source, the half-light radius is computed via `sep.flux_radius(data_sub, x, y, 6.0*a, 0.5, subpix=5)` &mdash; this finds the radius enclosing 50% of the total flux within a circular aperture of radius `6a` (where `a` is the semi-major axis). The FWHM is defined as twice the half-light radius: `FWHM = 2 × r_half`.
4. **Source filtering**: keep sources that are unflagged, unsaturated, have `2.0 ≤ FWHM ≤ 15.0` px, and finite elongation. The reported FWHM is the median over all accepted sources.
5. **Elongation**: `a / b` (semi-major / semi-minor axis ratio) for each source; the median is reported.

#### QC flags

Each frame receives a comma-separated list of flags (or `OK` if none apply):

| Flag | Condition | Threshold |
|------|-----------|-----------|
| `LOW_SOURCES` | Too few good sources | n_good < 5 |
| `BAD_SEEING` | Poor seeing | FWHM > 6.0&Prime; |
| `ELONGATED` | Tracking or guiding issue | elongation > 1.3 |
| `SATURATED` | Excessive saturation | saturated fraction > 1% |

**Output:**
- `qc/frame_stats.json`: JSON array of per-frame statistics.
- `qc/frame_stats.csv`: CSV export.
- Console summary table: per-filter FWHM, sky level, source counts, flag breakdown.


### 2.4 Stage 3: Astrometric Calibration (`pipeline/astrometry.py`)

**Command:** `python -m pipeline astrometry --year 2025 --telescope T120`

Solves the World Coordinate System (WCS) for every reduced frame by plate-solving with `astrometry.net`'s `solve-field` and validating against Gaia DR3. This provides a per-frame WCS written into the FITS headers. A subsequent multi-frame refinement with SCAMP (Section 2.7) produces the `.head` files used for stacking.

#### 2.4.1 Initial pointing estimate

Before calling the solver, the pipeline obtains an approximate sky position:

1. Parse `OBJCTRA` / `OBJCTDEC` from the FITS header (sexagesimal format).
2. If the header lacks coordinates (common on T080), resolve the `OBJECT` name via a pre-loaded Simbad cache (`simbad_cache.json`, 330 entries). Object names are cleaned (strip telescope prefixes, normalise known aliases like `Coma_N4874` &rarr; `NGC 4874`).
3. As a further fallback, pre-compute per-object coordinates from the median `OBJCTRA`/`OBJCTDEC` across all frames of the same object.
4. If all methods fail, the solver attempts a blind solve (no position hint).

#### 2.4.2 Plate solving with `solve-field`

Each frame is solved using `astrometry.net`'s `solve-field`:

- **Pixel scale bounds**: 0.6&times; to 1.4&times; the nominal pixel scale (arcsec/px).
- **Image downsampling**: factor 2 for speed.
- **CPU limit**: 120 seconds per frame.
- **Search centre**: if an approximate position is available, the solver is given `--ra`, `--dec`, `--radius 5°` to constrain the search. Otherwise, a blind solve is attempted.

The solver outputs a `.wcs` file containing the solved WCS, from which the pipeline extracts the CD matrix, rotation angle, and pixel scale.

#### 2.4.3 Gaia DR3 validation

Each solve-field solution is validated by cross-matching against Gaia DR3:

1. **Gaia query**: cone search centred on the solved WCS centre, radius = 2.5 &times; the FOV diagonal. Magnitude limit G < 18 mag (or G < 20 mag if fewer than 20 stars returned). Positions are proper-motion-corrected from the Gaia epoch (J2016.0) to J2025.2, using per-star proper motions (`pmra`, `pmdec`).
2. **Source extraction**: `sep` at 5.0&sigma; on the science frame, keeping the 200 brightest unflagged, unsaturated, non-elongated (`a/b < 3`) sources.
3. **Cross-match**: project extracted sources to sky coordinates via the solved WCS, match against Gaia positions with a 2&Prime; tolerance using a KD-tree.
4. **RMS**: the astrometric RMS is computed as `sqrt(mean(dRA² × cos²(Dec) + dDec²))` over all matched pairs, in arcseconds.

#### 2.4.4 Quality classification

Each frame receives a status based on the Gaia-validated RMS:

| Status | Condition | Interpretation |
|--------|-----------|----------------|
| `OK` | RMS &le; 3.0&Prime; | Good astrometric solution |
| `HIGH_RMS` | 3.0&Prime; < RMS &le; 10.0&Prime; | Usable but degraded (poor seeing, sparse field) |
| `FAILED` | RMS > 10.0&Prime; or solve-field did not converge | No reliable WCS |

Frames with `ASTRSTAT = FAILED` are excluded from subsequent photometric calibration and stacking. Frames with `HIGH_RMS` are included in photometry but may be rejected during stacking.

#### 2.4.5 WCS output

The following keywords are written to each FITS header:

| Keyword | Content |
|---------|---------|
| `CTYPE1/2` | `RA---TAN` / `DEC--TAN` |
| `CRPIX1/2` | Reference pixel (image centre) |
| `CRVAL1/2` | Reference RA / Dec (ICRS, degrees) |
| `CD1_1/1_2/2_1/2_2` | CD rotation/scale matrix |
| `EQUINOX` | 2000.0 |
| `RADESYS` | ICRS |
| `ASTRRMS` | Astrometric RMS vs Gaia (arcsec) |
| `ASTRNMAT` | Number of matched Gaia sources |
| `ASTRSTAT` | `OK`, `HIGH_RMS`, or `FAILED` |

**Reports:** `astrometry/astrom_report.json` and `.csv` containing per-frame results with filename, object, filter, status, RMS, number of matches, rotation, and pixel scale.

**Caching:** Gaia queries are cached as JSON files in `astrometry/gaia_cache/`, keyed by rounded (RA, Dec, radius, mag_limit). This allows re-runs without repeated TAP queries.


### 2.5 Stage 4: Photometric Calibration (`pipeline/photometry.py`)

**Command:** `python -m pipeline photometry --year 2025 --telescope T120`

Determines a per-frame photometric zero point by cross-matching aperture photometry against the Pan-STARRS1 DR2 catalogue.

#### 2.5.1 Source extraction and aperture photometry

For each frame with a valid astrometric solution (`ASTRSTAT` = `OK` or `HIGH_RMS`):

1. **Background subtraction**: `sep.Background()` computes a smooth sky model; the background-subtracted image is used for photometry.
2. **Source detection**: `sep.extract()` at 5.0&sigma; with minimum area 5 pixels.
3. **Source filtering**: keep unflagged sources below 90% of saturation, with aspect ratio `a/b < 3.0`. Require at least 5 sources.
4. **Aperture photometry**: `sep.sum_circle(data_sub, x, y, r_aper, err=bkg.globalrms, subpix=5)` where the aperture radius is:

$$r_{\rm aper} = 3.0 \times \mathrm{FWHM}_{\rm px}$$

The FWHM is taken from the QC stage; if unavailable, a default of 4.0 px is used.

5. **Instrumental magnitudes**: normalised to counts per second:

$$m_{\rm inst} = -2.5 \log_{10}\!\left(\frac{f}{\Delta t}\right)$$

where $f$ is the total flux in the aperture (ADU) and $\Delta t$ is the exposure time (seconds).

#### 2.5.2 PS1 cross-match

1. **PS1 query**: VizieR TAP query to Pan-STARRS1 DR2 (`II/349/ps1`) in a cone of radius 0.75 &times; the FOV, limited to stars with 13.0 < mag < 20.0 in the relevant PS1 band (up to 5000 rows).
2. **Cross-match**: KD-tree nearest-neighbour in `(RA × cos(Dec), Dec)` space with a 5.0&Prime; match radius.
3. **Duplicate removal**: if multiple extracted sources match the same PS1 star, only the closest match is kept.

#### 2.5.3 Zero-point fitting

The zero point is defined as:

$$\mathrm{ZP} = \mathrm{median}(m_{\rm PS1} - m_{\rm inst})$$

computed via iterative sigma clipping (5 iterations, 2.5&sigma; threshold using MAD-based robust standard deviation). Minimum 5 stars required.

**Output keywords written to FITS header:**

| Keyword | Content |
|---------|---------|
| `PHOTZP` | Zero point (mag) |
| `PHOTZPER` | Zero point uncertainty (mag) |
| `PHOTNSTR` | Number of stars used |
| `PHOTNMAT` | Number of cross-matches |
| `PHOTCAT` | Catalogue and band (e.g., `PS1_gmag`) |
| `PHOTSTAT` | `OK`, `SKIPPED`, or failure reason |
| `PHOTAPER` | Aperture radius used (pixels) |

**Magnitude convention:** the calibrated magnitude of a source is:

$$m = -2.5 \log_{10}\!\left(\frac{f}{\Delta t}\right) + \mathrm{PHOTZP}$$

**Caching:** PS1 queries are cached as JSON in `photometry/ps1_cache/`.

**Reports:** `photometry/phot_report.json` and `.csv`.


### 2.6 Stage 5: PSF Modelling (`pipeline/psf.py`)

**Command:** `python -m pipeline psf --year 2025 --telescope T120`

Builds a spatially-varying PSF model for each reduced science frame using SExtractor and PSFex.

#### 2.6.1 SExtractor

SExtractor is run on each frame to produce an LDAC (Leiden Data Analysis Centre) catalogue containing postage stamps (`VIGNET`) of detected sources. Key parameters:

| Parameter | Value |
|-----------|-------|
| Detection threshold | 5.0&sigma; |
| Vignet size | 35 &times; 35 pixels |
| Background mesh | 64 pixels |
| Deblending | 32 thresholds, min contrast 0.005 |
| Photometry aperture | 3.0 &times; FWHM / pixel_scale |

#### 2.6.2 PSFex

PSFex ingests the LDAC catalogue and fits a spatially-varying PSF model:

| Parameter | Value |
|-----------|-------|
| Basis | PIXEL_AUTO |
| PSF size | 25 &times; 25 pixels |
| Spatial variation | Degree 2 polynomial in (X, Y) |
| Star selection | FWHM 2&ndash;10&Prime;, ellipticity < 0.3, S/N > 20 |

The output is a `.psf` binary file that encodes the PSF as a set of basis images whose coefficients vary as a quadratic function of position across the detector. This captures focus variations, coma, and other field-dependent aberrations.

**Output statistics** (parsed from PSFex XML): mean FWHM, mean ellipticity, mean &chi;&sup2;, number of stars loaded and accepted.

**Reports:** `psf/psf_report.json` and `.csv`.


### 2.7 Stage 6: SCAMP Astrometric Refinement (`scripts/run_scamp.py`)

**Command:** `python scripts/run_scamp.py [--target TARGET] [--years YEAR ...] [--force] [--dry-run]`

SCAMP performs a multi-frame astrometric refinement against Gaia DR3, producing external `.head` header files that replace the per-frame WCS solutions from `solve-field`. These `.head` files are the WCS used by SWarp for stacking (Section 3).

#### 2.7.1 Why SCAMP?

The per-frame `solve-field` solutions (Section 2.4) are independent — each frame is solved in isolation. SCAMP instead solves all frames of a given target/filter group simultaneously, enforcing a consistent instrument model across exposures. This yields:

- Sub-pixel relative alignment between frames of the same group (critical for median stacking).
- A single distortion model shared across the instrument, reducing per-frame noise in the WCS.
- Direct calibration against Gaia DR3 proper-motion-corrected positions.

#### 2.7.2 SExtractor LDAC catalogue generation

For each input frame, SExtractor produces a FITS_LDAC catalogue containing windowed centroid positions (`XWIN_IMAGE`, `YWIN_IMAGE`) and associated errors. Multi-extension FITS files (from `ccdproc`) are first reduced to single-HDU files to avoid confusing SCAMP. Key SExtractor parameters:

| Parameter | Value |
|-----------|-------|
| Detection threshold | 5.0σ |
| Catalogue type | FITS_LDAC |
| Background mesh | 64 pixels |
| Deblending | 32 thresholds, min contrast 0.005 |
| Seeing FWHM | 3.0″ |

The output parameters include `XWIN_IMAGE`, `YWIN_IMAGE`, windowed position errors, `FLUX_AUTO`, `MAG_AUTO`, `FLAGS`, `FLUX_RADIUS`, `ELONGATION`, and `SNR_WIN`.

#### 2.7.3 SCAMP multi-frame solution

Frames are grouped by target and filter (filters are not mixed). SCAMP is invoked with:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `ASTREF_CATALOG` | GAIA-DR3 | External reference catalogue |
| `DISTORT_DEGREES` | 1 | First-order (affine) distortion model |
| `STABILITY_TYPE` | INSTRUMENT | Assume distortion is stable across exposures |
| `CROSSID_RADIUS` | 3.0″ | Maximum match radius for source identification |
| `POSITION_MAXERR` | 2.0″ | Maximum allowed positional error |
| `POSANGLE_MAXERR` | 2.0° | Maximum allowed position angle error |
| `PIXSCALE_MAXERR` | 1.05 | Maximum allowed pixel scale ratio error |
| `SN_THRESHOLDS` | 5.0, 20.0 | S/N for detection and high-quality stars |
| `ASTREFMAG_LIMITS` | 10.0, 19.0 | Reference catalogue magnitude range |
| `MOSAIC_TYPE` | UNCHANGED | No mosaic grouping |
| `SOLVE_PHOTOM` | N | Photometric solution handled separately |

#### 2.7.4 Output `.head` files

SCAMP writes one `.head` file per input frame, saved alongside the original reduced FITS file (e.g., `reduced/20250123/M67_r_001.head`). Each `.head` file contains a complete FITS header with the refined WCS (CTYPE, CRPIX, CRVAL, CD matrix, PV distortion terms).

These `.head` files are read by SWarp as "external header" overrides — SWarp uses the `.head` WCS in preference to the WCS embedded in the FITS header itself. The stacking script (Section 3.4) appends `FLXSCALE` to these `.head` files before invoking SWarp.

#### 2.7.5 CD matrix validation

Before stacking, each `.head` file is validated by checking the CD matrix determinant against the expected value for the T120:

$$|\det(\mathbf{CD})| \approx \left(\frac{0.770}{3600}\right)^2 \approx 4.57 \times 10^{-8} \;\mathrm{deg}^2$$

Frames with $|\det(\mathbf{CD})|$ deviating by more than a factor of 10 from this value are rejected as corrupt SCAMP solutions (typically caused by insufficient reference stars or failed convergence).

**Report:** `scamp_report.json` in the project root, containing per-group status, number of `.head` files produced, and external RMS values.


## 3. Image Stacking (`scripts/stack_all.py`)

**Command:** `python scripts/stack_all.py [--target TARGET] [--force] [--dry-run]`

The stacking script produces deep, pixel-aligned, photometrically calibrated coadded images from all T120 data (2018&ndash;2025) using SWarp. All stacks are placed on a common zero point of ZP = 30.0.

### 3.1 Target inventory

The script calls `pipeline.stacking.build_target_inventory()` which scans all reduced frames, normalises target names (via `TARGET_ALIASES`), filters out calibration frames and solar system objects, and groups frames by canonical target name.

**Merged groups** combine spatially overlapping targets onto a single WCS grid:

| Group name | Constituent targets |
|------------|---------------------|
| Coma | Abell1656, Abell1656_border, NGC4874 |
| M105_group | M105, NGC3384 |
| M38_group | M38, NGC1912 |
| Markarian | M84, M86 |

**Skipped targets** (too few frames or single-band): J0254, Gaia18aak, IC10, NGC628, Abell1228, NGC4636, NGC4278, NGC5322.

### 3.2 Frame selection

For each target group:

1. **Require `.head` files**: only frames with SCAMP astrometric solutions (`.head` external header files) are included.
2. **Spatial outlier rejection**: frames whose `CRVAL` position (from the `.head` file) is more than 0.5&deg; from the group median are rejected.
3. **CD matrix validation**: the `.head` file's CD matrix determinant is checked against the expected value of `(0.770/3600)^2`. Files with `|det(CD)|` deviating by more than a factor of 10 from expected are rejected (catches corrupt SCAMP solutions).
4. **Per-filter ZP outlier rejection**: within each filter group, frames whose PHOTZP is more than 2.0&sigma; (MAD-based) from the filter median are rejected.
5. **Minimum**: at least 2 frames per filter are required.

### 3.3 Shared WCS grid

A common pixel grid is computed across **all filters** for each target group:

1. For each surviving frame, the four corner sky coordinates are computed from the `.head` WCS and the frame dimensions.
2. The bounding box in (RA, Dec) is computed, with RA wrapping handled near 0&deg;/360&deg;.
3. The grid centre is the midpoint of the bounding box; the grid size is:

$$n_x = \left\lceil \frac{\Delta\alpha \cos\delta}{\theta_{\rm px}} \right\rceil + 40, \qquad n_y = \left\lceil \frac{\Delta\delta}{\theta_{\rm px}} \right\rceil + 40$$

where &theta;<sub>px</sub> = 0.770&Prime;/px is the output pixel scale and 40 pixels of padding are added on each axis.

This ensures that pixel (*i*, *j*) maps to the same sky position in every filter band &mdash; a requirement for multi-band photometry (e.g., SExtractor dual-mode).

### 3.4 Photometric flux scaling

Each frame is scaled to a target zero point of ZP = 30.0 using the `FLXSCALE` keyword. Because the pipeline's zero points are defined in counts per second (see Section 2.5.1), the flux scale factor must also account for the exposure time:

$$\mathrm{FLXSCALE} = \frac{10^{0.4\,(30 - \mathrm{PHOTZP})}}{t_{\rm exp}}$$

This ensures that in the stacked image, the calibrated magnitude of a source is:

$$m = -2.5 \log_{10}(f_{\rm stack}) + 30.0$$

without needing to know the exposure time. The `FLXSCALE` value is appended to each frame's `.head` file (before the `END` card) so that SWarp reads it via the `FSCALE_KEYWORD` mechanism.

### 3.5 Weight maps

For each input frame, the master flat-field is used as a weight map:

- The `FLATCORR` header keyword in the reduced frame identifies which filter's flat was applied.
- The flat file at `master_flat/master_flat_{year}_T120_2x2_{flatcorr}.fits` is symlinked as `{stem}.weight.fits` in the SWarp working directory.
- SWarp is invoked with `-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX .weight.fits`.

Since the master flats are normalised to median 1.0, pixel values directly encode relative sensitivity: vignetted or dusty regions have values < 1 and are correspondingly down-weighted in the stack.

For frames lacking a flat (`FLATCORR = 'NONE'`), a uniform weight of 1.0 is used instead.

### 3.6 SWarp parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `COMBINE_TYPE` | MEDIAN | Robust to outliers, cosmic rays |
| `SUBTRACT_BACK` | Y | Remove residual sky gradients |
| `BACK_SIZE` | 256 | Background mesh size (pixels) |
| `RESAMPLING_TYPE` | LANCZOS3 | High-quality flux-conserving interpolation |
| `FSCALE_KEYWORD` | FLXSCALE | Photometric scaling keyword |
| `FSCALASTRO_TYPE` | NONE | No automatic flux scaling |
| `PIXEL_SCALE` | 0.770 | Output pixel scale (arcsec/px) |
| `PIXELSCALE_TYPE` | MANUAL | Use exact specified scale |
| `CENTER_TYPE` | MANUAL | Use computed grid centre |
| `PROJECTION_TYPE` | TAN | Gnomonic (tangent-plane) projection |

### 3.7 Output headers

Each stack header contains:

| Keyword | Content |
|---------|---------|
| `NCOMBINE` | Number of input frames |
| `FILTER` | Filter band |
| `OBJECT` | Target group name |
| `PHOTZP` | 30.0 (target zero point) |
| `COMBTYPE` | MEDIAN |
| `BACKSUB` | Y |
| `BACKSIZE` | 256 |
| `PIXSCALE` | 0.770 |
| `INP0001`&ndash;`INPnnnn` | Input frame filenames |

### 3.8 Output files

Stacks are written to `$OHP_DATA_ROOT/stacks/{group}/T120/`:

- `{group}_{filter}.fits` &mdash; coadded image.
- `{group}_{filter}.weight.fits` &mdash; combined weight map.

The current dataset produces **131 stacks** across **40 target groups** in 8 filter bands (B, V, R, g, i, H&alpha;, [O III], [S II]). Note that all r-band variants (r, r_Gunn, R, R_Cousins) are merged into a single canonical filter &ldquo;R&rdquo; and calibrated against PS1 rmag.


### 3.9 Per-stack PSF extraction (`scripts/extract_stack_psfs.py`)

**Command:** `python scripts/extract_stack_psfs.py [--target TARGET] [--force]`

After stacking, SExtractor and PSFex are run on each coadded image to produce a PSF model specific to that stack. This captures the effective PSF of the coadd, which includes contributions from the seeing of all input frames, astrometric alignment residuals, and Lanczos-3 resampling.

1. **SExtractor** extracts an LDAC catalogue from the stack, using the weight map for source weighting. The background mesh is set to 128 pixels (larger than the per-frame value of 64, appropriate for the larger stack images).
2. **PSFex** fits a spatially varying PSF model (degree-2 polynomial, 25 &times; 25 pixel basis).
3. The centre-of-field component is extracted, normalised to sum = 1, and saved as a simple 2D FITS image: `{group}_{filter}_psf.fits`.

The PSF header contains `FWHM_PX`, `FWHM_AS`, `EE50_AS`, `EE80_AS`, `ELLIPT`, `NSTARS`, and `CHI2`.

**Report:** `stacks/stack_psf_report.json` with per-stack PSF statistics.


## 4. Photometric Validation (`scripts/check_phot.py`)

**Command:** `python scripts/check_phot.py [--fix] [--plot] [--target TARGET]`

Validates the photometric calibration of all stacks by independently measuring zero points against Pan-STARRS1 DR2, using exactly the same methodology as the pipeline's per-frame photometry.

### 4.1 Measurement methodology

The validation script deliberately replicates the pipeline's photometric methodology to ensure any discrepancy reflects a genuine calibration error rather than a measurement inconsistency.

For each stack with a PS1-compatible filter:

1. **FWHM measurement**: identical to `pipeline/quality.py` &mdash; `sep.flux_radius()` half-light radius method (radius enclosing 50% of the flux within `6a`; FWHM = 2 &times; r_half). The FWHM of the stack PSF is typically broader than individual frames due to the combination of frames with different seeing.
2. **Aperture photometry**: identical to `pipeline/photometry.py` &mdash; `sep.sum_circle()` with aperture radius = 3.0 &times; FWHM (the same `PHOT_APERTURE_MULT` multiplier used in the pipeline). Zero-coverage regions (pixel value = 0) are masked before source extraction.
3. **Instrumental magnitudes**: because the stacks are flux-scaled (Section 3.4), pixel values are already in calibrated flux units:

$$m_{\rm inst} = -2.5 \log_{10}(f_{\rm stack})$$

No exposure-time division is needed. This contrasts with the per-frame pipeline convention ($m = -2.5 \log_{10}(f / \Delta t) + \mathrm{ZP}$), because the FLXSCALE factor has already absorbed both the ZP correction and the exposure-time normalisation.

4. **PS1 cross-match**: cone search from the stack centre, cross-matched via KD-tree with 3.0&Prime; tolerance.
5. **Zero-point fit**: iterative sigma-clipped median of $(m_{\rm PS1} - m_{\rm inst})$, using the same `compute_zeropoint()` function from `pipeline.photometry` (5 iterations, 2.5&sigma; MAD-based clipping, minimum 5 stars).

#### Methodology verification

To confirm that the validation and pipeline measurements are consistent, the methodology was tested on individual reduced frames before stacking. For a single frame with known pipeline zero point (e.g., PHOTZP = 24.50, EXPTIME = 2.0 s), the validation script &mdash; which does not divide by exposure time &mdash; measured ZP = 25.25, a difference of exactly $2.5 \times \log_{10}(2.0) = 0.75$ mag. This confirmed that the only difference between the two measurements is the counts-per-second normalisation convention used by the pipeline, and motivated the FLXSCALE formula in Section 3.4 (which divides by $t_{\rm exp}$).

### 4.2 Validation results

Of 131 total stacks, 117 broadband stacks were validated (narrowband and stacks with too few PS1 matches excluded):

| Metric | Value |
|--------|-------|
| Median ZP offset | &minus;0.001 mag |
| Standard deviation | 0.20 mag |
| Median RMS (per stack) | 0.82 mag |
| Range | &minus;0.39 to +1.65 mag |

The residual scatter is dominated by colour terms (particularly B &rarr; gmag and V &rarr; gmag conversions) and crowded-field photometry in clusters.

### 4.3 Correction mode

With `--fix`, the script applies a multiplicative correction to bring each stack to exactly ZP = 30.0:

$$f_{\rm corrected} = f_{\rm stack} \times 10^{0.4 \times \mathrm{offset}}$$

The weight map is correspondingly adjusted by dividing by the factor squared. Corrections smaller than 0.01 mag are skipped. Updated header keywords: `ZPOFFSET`, `ZPFACTOR`, `ZPMEAS`, `ZPNSTAR`, `ZPRMS`.

### 4.4 Diagnostic plots

With `--plot`, a two-panel diagnostic plot is saved for each stack:
- **Left panel**: PS1 catalogue magnitude vs calibrated instrumental magnitude (1:1 plot).
- **Right panel**: residuals (PS1 &minus; calibrated) vs magnitude, with &plusmn;RMS lines.

Plots are saved to `{group}/T120/plots/{group}_{filter}_phot.png`.


## 5. External Dependencies

| Tool | Version/Path | Purpose |
|------|-------------|---------|
| Python | 3.12 (`micromamba activate astropy-3.12`) | Runtime |
| astropy | &mdash; | FITS I/O, WCS, coordinates |
| ccdproc | &mdash; | Bias/flat combination and correction |
| sep | &mdash; | Source extraction and aperture photometry |
| numpy, scipy | &mdash; | Numerical computation |
| astroquery | &mdash; | Simbad and Gaia TAP queries |
| solve-field | On PATH | Astrometric plate solving |
| SExtractor | `sex` on PATH | LDAC catalogue generation for PSFex and SCAMP |
| PSFex | `/Users/hjmcc/opt/astromatic/bin/psfex` | PSF modelling |
| SWarp | `/Users/hjmcc/opt/astromatic/bin/swarp` | Image resampling and coaddition |
| SCAMP | `scamp` on PATH | Multi-frame astrometric refinement (Section 2.7, `scripts/run_scamp.py`) |


## 6. Catalogue Caching

All external catalogue queries are cached to local JSON files to survive re-runs and avoid redundant network access:

| Catalogue | Cache location | Key format |
|-----------|---------------|------------|
| Simbad | `ohp-archive/simbad_cache.json` | Object name |
| Gaia DR3 | `{year}/{telescope}/astrometry/gaia_cache/` | `gaia_{ra}_{dec}_{radius}_{maglim}.json` |
| PS1 DR2 (per-frame) | `{year}/{telescope}/photometry/ps1_cache/` | `ps1_{ra}_{dec}_{radius}_{band}.json` |
| PS1 DR2 (stacks) | `stacks/ps1_cache/` | Same format |


## 7. Quick Reference: Running the Full Pipeline

```bash
# Activate environment
micromamba activate astropy-3.12

# Full pipeline for one year/telescope
python -m pipeline all --year 2025 --telescope T120

# Or run stages individually
python -m pipeline scan        --year 2025 --telescope T120
python -m pipeline calibrate   --year 2025 --telescope T120
python -m pipeline qc          --year 2025 --telescope T120
python -m pipeline astrometry  --year 2025 --telescope T120
python -m pipeline photometry  --year 2025 --telescope T120
python -m pipeline psf         --year 2025 --telescope T120

# Stack all targets
python scripts/stack_all.py

# Stack one target
python scripts/stack_all.py --target M67

# Validate photometry of stacks
python scripts/check_phot.py --plot

# Apply photometric corrections if needed
python scripts/check_phot.py --fix

# Re-process with --force to overwrite existing results
python -m pipeline photometry --year 2025 --telescope T120 --force
python scripts/stack_all.py --force
```
