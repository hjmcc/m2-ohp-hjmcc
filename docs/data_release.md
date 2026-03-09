# OHP M2 T120 Stacked Image Data Release

## Overview

This data release provides deep, multi-band stacked images from the OHP (Observatoire de Haute-Provence) 1.2 m Newton telescope (T120), covering 8 years of M2 student observing campaigns (2018--2025). Raw CCD frames have been processed through a fully automated pipeline (bias subtraction, flat-fielding, astrometric calibration, photometric calibration) and combined into 131 coadded images across 40 target fields.

All stacks are placed on a common photometric zero point of **ZP = 30.0** and validated against Pan-STARRS1 DR2. Within each target group, all filter bands share a common pixel grid, enabling direct pixel-to-pixel multi-band analysis.


## Data Summary

| Property | Value |
|----------|-------|
| Telescope | OHP 1.2 m Newton (T120) |
| CCD | Andor 1k x 1k, 2x2 binning |
| Output pixel scale | 0.770 arcsec/px |
| Projection | Gnomonic (TAN) |
| Coordinate system | ICRS, J2000 |
| Photometric zero point | 30.0 mag |
| Number of target groups | 40 |
| Number of stacks | 131 (118 broadband + 13 narrowband) |
| Total input frames | 1268 (from 2894 reduced T120 frames) |
| Observing years | 2018, 2019, 2020, 2021, 2023, 2024, 2025 |
| Combination method | Median |
| Resampling | Lanczos-3 |


## Filter Inventory

| Filter | N stacks | PS1 calibration band |
|--------|----------|---------------------|
| B | 24 | gmag |
| V | 21 | gmag |
| R | 11 | rmag |
| g | 26 | gmag |
| r | 28 | rmag |
| i | 10 | imag |
| Halpha | 9 | -- (no PS1 match) |
| [O III] | 2 | -- |


## Target Catalogue

Targets are grouped by sky position. Spatially overlapping fields are merged into single groups with a shared WCS grid.

### Merged groups

| Group name | Constituent targets | Description |
|------------|---------------------|-------------|
| Coma | Abell1656, Abell1656_border, NGC4874 | Coma galaxy cluster |
| M105_group | M105, NGC3384 | Leo I group |
| M38_group | M38, NGC1912 | Open cluster |
| Markarian | M84, M86 | Markarian's Chain (Virgo) |

### Full target list

| Target | Filters | N frames | Type |
|--------|---------|----------|------|
| Abell1703 | B, V, g, i, r | 11 | Galaxy cluster (lensing) |
| Abell2744 | B, V, g, i, r | 96 | Galaxy cluster (Pandora) |
| Aspirateur | B, V, r | 12 | Galaxy |
| AWM7 | B, R | 20 | Galaxy cluster |
| Coma | g, r | 178 | Galaxy cluster |
| M1 | B, Ha, [O III] | 8 | Supernova remnant (Crab) |
| M13 | B, V, r | 13 | Globular cluster |
| M38_group | B, V, g, r | 46 | Open cluster |
| M48 | B, V, g, i, r | 72 | Open cluster |
| M49 | B, g, i, r | 33 | Elliptical galaxy (Virgo) |
| M51 | B, Ha, [O III], V, g, r | 44 | Spiral galaxy (Whirlpool) |
| M57 | B, Ha, V, r | 17 | Planetary nebula (Ring) |
| M60 | Ha, g, r | 26 | Elliptical galaxy (Virgo) |
| M63 | B, V, r | 10 | Spiral galaxy (Sunflower) |
| M64 | B, V, r | 9 | Spiral galaxy (Black Eye) |
| M67 | B, V, g, i, r | 107 | Open cluster |
| M81 | B, Ha, V, g, i, r | 93 | Spiral galaxy (Bode's) |
| M82 | B, Ha, V, g, i, r | 59 | Starburst galaxy |
| M87 | B, g, i, r | 31 | Elliptical galaxy (Virgo A) |
| M91 | Ha, g, i, r | 24 | Barred spiral galaxy |
| M92 | B, V, r | 15 | Globular cluster |
| M97 | B, V, r | 14 | Planetary nebula (Owl) |
| M99 | B, V, r | 12 | Spiral galaxy |
| M104 | B, R, V | 13 | Spiral galaxy (Sombrero) |
| M105_group | g, i, r | 32 | Elliptical galaxy + NGC3384 |
| M108 | B, R, V | 11 | Edge-on spiral galaxy |
| Markarian | g, r | 40 | Markarian's Chain |
| NGC2158 | B, V, r | 7 | Open cluster |
| NGC2301 | R, g | 16 | Open cluster |
| NGC2420 | R, V, g | 25 | Open cluster |
| NGC2539 | B, R, V | 25 | Open cluster |
| NGC2768 | R, g | 19 | Lenticular galaxy |
| NGC3377 | R, g | 17 | Elliptical galaxy |
| NGC3718 | Ha, R, g | 39 | Peculiar galaxy |
| NGC4008 | g, r | 10 | Elliptical galaxy |
| NGC4365 | R, g | 14 | Elliptical galaxy |
| NGC4494 | g, r | 14 | Elliptical galaxy |
| NGC5557 | g, r | 16 | Elliptical galaxy |
| NGC5566 | Ha, R, g | 10 | Barred spiral galaxy |
| NGC5813 | g, r | 10 | Elliptical galaxy |


## File Format

### Stacked images

Each stack is a single-extension FITS file with 32-bit floating-point pixel values. The WCS is a standard gnomonic (TAN) projection with a CD rotation matrix.

**File naming:** `{target}_{filter}.fits`

**Key header keywords:**

| Keyword | Description |
|---------|-------------|
| `OBJECT` | Target group name |
| `FILTER` | Filter band |
| `PHOTZP` | Photometric zero point (= 30.0) |
| `NCOMBINE` | Number of input frames |
| `COMBTYPE` | Combination method (MEDIAN) |
| `BACKSUB` | Background subtracted (Y) |
| `BACKSIZE` | Background mesh size (256 px) |
| `PIXSCALE` | Output pixel scale (0.770 arcsec/px) |
| `FWHM_PX` | Stack PSF FWHM in pixels |
| `FWHM_AS` | Stack PSF FWHM in arcsec |
| `PHOTNSTR` | Number of stars used in ZP validation |
| `PHOTRMS` | ZP validation RMS (mag) |
| `CTYPE1/2` | `RA---TAN` / `DEC--TAN` |
| `CRPIX1/2` | Reference pixel |
| `CRVAL1/2` | Reference sky position (degrees, ICRS) |
| `CD1_1` ... `CD2_2` | CD rotation/scale matrix |
| `INP0001`--`INPnnnn` | Input frame filenames |

**Photometric convention:** the calibrated magnitude of a source measured by simple aperture photometry on a stack is:

```
m = -2.5 * log10(flux_ADU) + 30.0
```

No exposure-time normalisation is needed; the flux scaling is already applied.

**Units:** pixel values are in flux units scaled to ZP = 30.0. Specifically, a source with magnitude *m* produces a total aperture flux of `10^(0.4 * (30 - m))` ADU. The background has been subtracted by SWarp (mesh size 256 px).

### Weight maps

Each stack has a corresponding `{target}_{filter}.weight.fits` file. Weight values encode the effective relative sensitivity at each pixel, derived from the normalised master flat-fields of the input frames and the SWarp combination process. Pixels outside the frame overlap region have weight = 0.

**Masking bad pixels:** regions with `weight = 0` lie outside the frame overlap area and contain no data. Use the weight map to mask these regions:

```python
weight = fits.getdata('M67_r.weight.fits')
mask = weight == 0  # True where there is no data
```

### PSF models

Each stack has a corresponding PSF image: `{target}_{filter}_psf.fits`. These are extracted directly from each stack using SExtractor + PSFex, so they represent the actual effective PSF of that specific coadded image (accounting for the seeing, astrometric alignment residuals, and resampling of all input frames).

**File naming:** `{target}_{filter}_psf.fits`

| Property | Value |
|----------|-------|
| Array size | 25 x 25 pixels |
| Pixel scale | 0.770 arcsec/px |
| Normalisation | Sum = 1.0 |
| Spatial variation | Centre-of-field component |
| Median FWHM across stacks | 2.83 arcsec (range 1.69--5.59 arcsec) |

**PSF header keywords:**

| Keyword | Description |
|---------|-------------|
| `FWHM_PX` | PSF FWHM in pixels |
| `FWHM_AS` | PSF FWHM in arcsec |
| `EE50_AS` | 50% encircled energy radius (arcsec) |
| `EE80_AS` | 80% encircled energy radius (arcsec) |
| `ELLIPT` | Mean ellipticity (0 = circular) |
| `NSTARS` | Number of stars used for PSF model |

**When to use:** the per-stack PSF is the correct PSF for any analysis on that stack image: PSF photometry, PSF-convolved model fitting (Sérsic profiles with PetroFit or Galfit), deconvolution, and aperture corrections.

A global mean PSF (`T120_mean_psf.fits`) is also provided in the stacks root directory. It is the median of 2754 per-frame PSF models across all years (2018--2025), with FWHM = 3.44 arcsec. Use it only when a generic T120 PSF is needed (e.g., simulation inputs) -- for science on individual stacks, always use the per-stack PSF.

### Pixel alignment

Within each target group, all filter stacks share identical `NAXIS1`, `NAXIS2`, `CRVAL1`, `CRVAL2`, `CRPIX1`, `CRPIX2`, and CD matrix values. Pixel (i, j) maps to the same sky position in every band. This enables:

- Direct RGB colour composites without resampling
- SExtractor dual-image mode (detect in one band, measure in all)
- Pixel-by-pixel colour maps and SED fitting


## Data Quality

### Photometric calibration

All broadband stacks (B, V, R, g, r, i) are calibrated against Pan-STARRS1 DR2 using aperture photometry with radius = 3 x FWHM and iterative sigma-clipped median zero-point fitting.

| Metric | Value |
|--------|-------|
| Stacks validated | 117 / 118 broadband |
| Median ZP offset from 30.0 | +0.001 mag |
| Standard deviation of offsets | 0.20 mag |
| Median internal RMS | 0.82 mag |

The internal RMS is dominated by colour terms (particularly B and V mapped to PS1 gmag) and crowded-field photometry in dense clusters. For the g, r, and i bands (native PS1 matches), typical per-stack RMS is 0.5--0.7 mag.

Three M64 stacks failed validation due to insufficient PS1 matches (crowded field with a dominant galaxy). Narrowband stacks (Halpha, [O III]) have no PS1 counterpart and are not photometrically validated.

### Astrometric calibration

Input frame WCS solutions are obtained via `astrometry.net` solve-field and validated against Gaia DR3. Only frames with valid SCAMP-refined `.head` files are included in the stacking. Frames with aberrant CD matrices (determinant deviating >10x from expected) are rejected.

### Known limitations

- **Colour terms**: B and V magnitudes are mapped to PS1 gmag without colour-term correction. For precise photometry in these bands, users should derive their own colour transformations.
- **Crowded fields**: aperture photometry in dense clusters (M13, M67, M92) is affected by blending. PSF photometry is recommended for these targets.
- **Narrowband stacks**: Halpha, [O III], and [S II] stacks are combined without photometric flux scaling for frames lacking a broadband zero point. Their relative calibration across frames relies on exposure-time normalisation only.
- **M81 B-band**: the M81_B stack has a +1.65 mag ZP offset, likely due to poor input frame ZPs in this band/field combination.
- **No T080 stacks**: only T120 data are included in this release. T080 data are processed through the per-frame pipeline but not stacked.
- **Year 2022**: no T120 observations were taken in 2022.


## Processing Summary

Each input frame passed through the following stages before stacking:

1. **Bias subtraction** using a master bias (averaged, 3-sigma-clipped combination of all bias frames for the year).
2. **Flat-field correction** using a master flat (median-combined, normalised to median = 1.0) matched by filter.
3. **Quality control**: FWHM, elongation, saturation fraction, and source count measured with `sep`. Frames flagged for bad seeing (>6 arcsec), low source count (<5), or excessive saturation (>1%) are noted but not automatically excluded from stacking.
4. **Astrometric calibration**: WCS solved via `astrometry.net` solve-field, validated against Gaia DR3 (proper-motion corrected to J2025.2). Typical RMS: 0.3--1.0 arcsec.
5. **Photometric calibration**: aperture photometry (r = 3 x FWHM) cross-matched against Pan-STARRS1 DR2 (13 < mag < 20). Zero point fitted as sigma-clipped median. Magnitudes are in counts per second: `m = -2.5 log10(flux/exptime) + ZP`.
6. **SCAMP astrometric refinement**: multi-frame astrometric solution against Gaia DR3, producing `.head` external header files with refined WCS.

Stacking was performed with SWarp using:

- **Flux scaling**: `FLXSCALE = 10^(0.4 * (30 - ZP)) / exptime`, written into `.head` files. This accounts for the pipeline's counts-per-second magnitude convention.
- **Weight maps**: normalised master flat-fields used as pixel-level weight maps.
- **Shared grid**: a common (RA, Dec) bounding box computed from all frames across all filters, with 40-pixel padding.
- **Outlier rejection**: spatial outliers (>0.5 deg from median), ZP outliers (>2-sigma MAD), and corrupt SCAMP solutions (CD matrix validation) removed before stacking.


## Image Quality

### Seeing

The effective FWHM of each stack is measured by PSFex directly on the coadded image and recorded in the PSF file header and in `stack_psf_report.json`. The stack header also contains the FWHM (`FWHM_AS` keyword). The FWHM reflects the combined effect of atmospheric seeing across all input frames (which may span multiple nights and years), astrometric alignment residuals, and resampling.

| Metric | Value |
|--------|-------|
| Median stack FWHM (PSFex) | 2.83 arcsec |
| Best seeing | 1.69 arcsec |
| Worst seeing | 5.59 arcsec |

### Depth

The depth varies by target and filter, depending on the number of input frames and their exposure times. As a rough guide, stacks with NCOMBINE > 10 typically reach surface brightnesses of ~26--27 mag/arcsec² (1σ per pixel) in the g and r bands. Point-source limiting magnitudes can be estimated from the weight maps and background noise.


## Quickstart: Python Examples

### Reading a stack

```python
from astropy.io import fits
from astropy.wcs import WCS

hdu = fits.open('M67/T120/M67_r.fits')[0]
data = hdu.data       # 2D numpy array, flux in ZP=30 units
header = hdu.header
wcs = WCS(header)

print(f"Target: {header['OBJECT']}, Filter: {header['FILTER']}")
print(f"FWHM: {header['FWHM_AS']:.1f} arcsec")
print(f"N frames combined: {header['NCOMBINE']}")
```

### Converting pixel values to magnitudes

```python
import numpy as np

# Magnitude of a source with total flux f (ADU) in the aperture:
flux = 1234.5  # from aperture photometry
mag = -2.5 * np.log10(flux) + 30.0

# Surface brightness (mag/arcsec²):
pixel_scale = 0.770  # arcsec/px
pixel_area = pixel_scale**2  # arcsec²
sb = -2.5 * np.log10(data / pixel_area) + 30.0
```

### Aperture photometry with `sep`

```python
import sep
from astropy.io import fits

data = fits.getdata('M67/T120/M67_r.fits').astype(float)
weight = fits.getdata('M67/T120/M67_r.weight.fits')

# Mask zero-weight pixels
data[weight == 0] = 0.0

# Background subtraction (SWarp already subtracted, but local refinement helps)
bkg = sep.Background(data)
data_sub = data - bkg

# Source detection
sources = sep.extract(data_sub, thresh=3.0, err=bkg.globalrms)

# Aperture photometry (r = 3 * FWHM)
fwhm_px = fits.getheader('M67/T120/M67_r.fits')['FWHM_PX']
r_aper = 3.0 * fwhm_px
flux, fluxerr, flag = sep.sum_circle(data_sub, sources['x'], sources['y'],
                                      r_aper, err=bkg.globalrms)

# Calibrated magnitudes
mag = -2.5 * np.log10(flux) + 30.0
```

### Making an RGB colour composite

```python
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# All bands are pixel-aligned — no reprojection needed
r = fits.getdata('M51/T120/M51_r.fits')
g = fits.getdata('M51/T120/M51_g.fits')
b = fits.getdata('M51/T120/M51_B.fits')

# Asinh stretch for display
def stretch(img, Q=10, minimum=0):
    img = np.clip(img - minimum, 0, None)
    return np.arcsinh(Q * img / img.max()) / np.arcsinh(Q)

rgb = np.stack([stretch(r), stretch(g), stretch(b)], axis=-1)
plt.imshow(rgb, origin='lower')
plt.show()
```

### SExtractor dual-image mode

Because all bands share the same pixel grid, you can detect sources in one band and measure fluxes in all bands:

```bash
# Detect in r-band, measure in g-band
sex M67_r.fits,M67_g.fits -c default.sex \
    -CATALOG_NAME M67_g_dualmode.cat \
    -WEIGHT_IMAGE M67_r.weight.fits,M67_g.weight.fits \
    -WEIGHT_TYPE MAP_WEIGHT,MAP_WEIGHT \
    -MAG_ZEROPOINT 30.0
```

### Using the PSF for model fitting

```python
from astropy.io import fits

# Load the PSF for this specific stack
psf = fits.getdata('M49/T120/M49_r_psf.fits')  # 25x25, sum=1.0
psf_hdr = fits.getheader('M49/T120/M49_r_psf.fits')
print(f"PSF FWHM: {psf_hdr['FWHM_AS']:.2f} arcsec, "
      f"ellipticity: {psf_hdr['ELLIPT']:.3f}")

# Example with PetroFit:
from petrofit.models import PSFConvolvedModel2D
from astropy.modeling.models import Sersic2D

sersic = Sersic2D(amplitude=1, r_eff=10, n=4, x_0=50, y_0=50)
psf_model = PSFConvolvedModel2D(sersic, psf=psf, oversample=4)
```

### Sky coordinate lookup

```python
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

wcs = WCS(fits.getheader('M87/T120/M87_r.fits'))

# Sky → pixel
coord = SkyCoord(ra=187.7059, dec=12.3911, unit='deg')
x, y = wcs.world_to_pixel(coord)

# Pixel → sky
ra, dec = wcs.pixel_to_world_values(100, 200)
```


## File Locations

```
$OHP_DATA_ROOT/stacks/
  T120_mean_psf.fits                            # global mean PSF (25x25, sum=1)
  phot_check_report.json                        # photometric validation results
  stack_psf_report.json                         # per-stack PSF extraction results
  stack_all_report.json                         # stacking log
  {target}/T120/
    {target}_{filter}.fits                      # science image
    {target}_{filter}.weight.fits               # weight map
    {target}_{filter}_psf.fits                  # PSF model for this stack (25x25, sum=1)
    plots/
      {target}_{filter}_phot.png                # photometric validation plot
```

Total data volume: 131 science images + 131 weight maps + 131 PSF models + 1 global mean PSF.


## Per-Frame Pipeline Data

In addition to the stacked images, the full per-frame pipeline output is available for advanced analysis:

```
$OHP_DATA_ROOT/{year}/T120/
  reduced/{night}/*.fits                        # calibrated individual frames
  reduced/{night}/*.head                        # SCAMP-refined WCS
  master_bias/master_bias_*.fits                # master bias
  master_flat/master_flat_*_{filter}.fits       # master flat per filter
  qc/frame_stats.json                           # per-frame quality metrics
  astrometry/astrom_report.json                 # per-frame astrometry results
  photometry/phot_report.json                   # per-frame photometry + ZPs
  psf/{night}/*.psf                             # per-frame PSFex models
  psf/psf_report.json                           # per-frame PSF statistics
```

### Per-frame PSF models

Each reduced frame has a corresponding PSFex model (`.psf` file) stored as a FITS binary table. The PSF model encodes a spatially varying PSF as a set of basis images whose coefficients vary as a degree-2 polynomial in detector position. To extract the centre-of-field PSF:

```python
from astropy.io import fits
with fits.open('2025/T120/psf/20250317/M67_-0003r60.psf') as hdul:
    psf_mask = hdul[1].data['PSF_MASK'][0]  # shape (6, 25, 25)
    psf_centre = psf_mask[0]                 # centre-of-field component
    psf_centre /= psf_centre.sum()           # normalise
```

The 6 basis images correspond to the polynomial terms: 1, x, y, x², xy, y². To evaluate the PSF at a specific detector position, combine the basis images with the polynomial evaluated at that position.

Key per-frame PSF statistics (from `psf_report.json`):

| Field | Description |
|-------|-------------|
| `fwhm_px` | PSF FWHM in pixels |
| `fwhm_arcsec` | PSF FWHM in arcsec |
| `ellipticity` | PSF elongation (0 = circular) |
| `chi2` | PSFex goodness of fit |
| `n_stars_accepted` | Stars used for PSF model |


## How to Cite

If you use these data, please cite the OHP M2 observing programme and acknowledge the Observatoire de Haute-Provence. Pipeline documentation is available at `docs/pipeline_description.md` in the repository.
