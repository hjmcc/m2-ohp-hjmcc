# OHP M2 Archive

Scanner and database for 8 years (2018-2025) of OHP M2 student observing data
across four telescopes (T120, T080, IRIS, T152). Pure stdlib Python — no astropy
or other dependencies required.

## Files

| File | Description |
|------|-------------|
| `scan_archive.py` | Main scanner script. Reads FITS headers, resolves targets via Simbad, and produces all output files. Run with `python3 scan_archive.py`. |
| `ohp_archive.csv` | Flat CSV database — one row per FITS file with header metadata (object, filter, exposure, coordinates, etc.). |
| `ohp_archive.json` | Aggregated JSON with per-year summaries, master target catalogue, filter statistics, and per-object file index. |
| `archive_summary.html` | Self-contained HTML dashboard (embedded CSS/JS/data) with local file paths. Includes sortable target tables with click-to-expand file lists, filter usage, and T152 spectroscopy log. |
| `ohp-m2-archive.html` | Web-deployable version of the dashboard. File paths are relative URLs with clickable download links, intended for `https://ohp.ias.universite-paris-saclay.fr/archive/`. |
| `simbad_cache.json` | Cached Simbad TAP query results to avoid re-querying known targets on each run. |
