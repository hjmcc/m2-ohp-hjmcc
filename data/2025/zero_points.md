# Approximate Photometric Zero Points — 2025 OHP Run

Zero points derived from PS1 DR2 cross-matching of all science frames from the
March 2025 run. Conditions were **not photometric** on most nights, so these
are rough estimates only (typical scatter ~1–2 mag frame-to-frame).

Convention: `mag = ZP - 2.5 * log10(counts / exptime)`

where `counts` is the total aperture flux in ADU and `exptime` is in seconds.

## T080 (80 cm Cassegrain)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| g'     | 20.6     | 1.0           | 38       |
| r'     | 21.4     | 1.2           | 35       |
| i'     | 21.2     | 1.3           | 24       |

## T120 (1.2 m Newton)

| Filter | ZP (mag) | Scatter (mag) | N frames |
|--------|----------|---------------|----------|
| B      | 25.3     | 1.7           | 37       |
| V      | 25.9     | 1.4           | 40       |
| R      | 25.0     | 1.5           | 112      |
| g      | 24.7     | 1.9           | 24       |

ZP values are medians across all calibrated frames. Scatter is the standard
deviation. Use these for order-of-magnitude flux estimates; do not rely on
them for precision photometry.
