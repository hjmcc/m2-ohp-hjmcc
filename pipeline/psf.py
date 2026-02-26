"""PSF modeling via SExtractor + PSFex: per-frame PSF extraction and quality stats."""

import csv
import json
import logging
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from astropy.io import fits

from . import config

log = logging.getLogger("pipeline")


# ── SExtractor config generation ────────────────────────────────────────

def _write_sex_config(tmpdir, catalog_path, telescope, fwhm_arcsec=None):
    """Write SExtractor .sex and .param files for LDAC output.

    Returns path to the .sex config file.
    """
    tel = config.TELESCOPES[telescope]
    pixel_scale = tel["pixel_scale"]
    saturation = tel["saturation"]
    vig = config.PSF_VIGNET_SIZE

    if fwhm_arcsec is None or not np.isfinite(fwhm_arcsec):
        fwhm_arcsec = 3.0  # reasonable default
    seeing_fwhm = fwhm_arcsec

    # .param file
    param_path = Path(tmpdir) / "default.param"
    param_path.write_text(
        "X_IMAGE\n"
        "Y_IMAGE\n"
        "FLUX_APER(1)\n"
        "FLUXERR_APER(1)\n"
        "MAG_AUTO\n"
        "FLUX_RADIUS\n"
        "FLAGS\n"
        "ELONGATION\n"
        "SNR_WIN\n"
        f"VIGNET({vig},{vig})\n"
    )

    # .sex config
    sex_path = Path(tmpdir) / "default.sex"
    sex_path.write_text(
        f"CATALOG_NAME     {catalog_path}\n"
        f"CATALOG_TYPE     FITS_LDAC\n"
        f"PARAMETERS_NAME  {param_path}\n"
        f"\n"
        f"DETECT_TYPE      CCD\n"
        f"DETECT_THRESH    {config.PSF_DETECT_THRESH}\n"
        f"ANALYSIS_THRESH  {config.PSF_DETECT_THRESH}\n"
        f"\n"
        f"FILTER           N\n"
        f"\n"
        f"DEBLEND_NTHRESH  32\n"
        f"DEBLEND_MINCONT  0.005\n"
        f"\n"
        f"CLEAN            Y\n"
        f"CLEAN_PARAM      1.0\n"
        f"\n"
        f"PHOT_APERTURES   {seeing_fwhm / pixel_scale * 3.0:.1f}\n"
        f"\n"
        f"SATUR_LEVEL      {saturation}\n"
        f"PIXEL_SCALE      {pixel_scale}\n"
        f"SEEING_FWHM      {seeing_fwhm}\n"
        f"\n"
        f"BACK_SIZE        64\n"
        f"BACK_FILTERSIZE  3\n"
        f"\n"
        f"WEIGHT_TYPE      NONE\n"
        f"CHECKIMAGE_TYPE  NONE\n"
        f"VERBOSE_TYPE     QUIET\n"
    )
    return str(sex_path)


def _write_psfex_config(tmpdir, psf_dir, telescope):
    """Write PSFex .psfex config file.

    Returns path to the .psfex config file.
    """
    tel = config.TELESCOPES[telescope]
    vig = config.PSF_VIGNET_SIZE

    psfex_path = Path(tmpdir) / "default.psfex"
    psfex_path.write_text(
        f"BASIS_TYPE       PIXEL_AUTO\n"
        f"PSF_SIZE         25,25\n"
        f"PSF_SAMPLING     0.0\n"
        f"\n"
        f"PSFVAR_KEYS      X_IMAGE,Y_IMAGE\n"
        f"PSFVAR_GROUPS    1,1\n"
        f"PSFVAR_DEGREES   {config.PSF_PSFVAR_DEGREES}\n"
        f"\n"
        f"SAMPLE_AUTOSELECT Y\n"
        f"SAMPLE_FWHMRANGE  2.0,10.0\n"
        f"SAMPLE_VARIABILITY 0.2\n"
        f"SAMPLE_MINSN     {config.PSF_SAMPLE_MINSN}\n"
        f"SAMPLE_MAXELLIP  {config.PSF_SAMPLE_MAXELLIP}\n"
        f"\n"
        f"PSF_DIR          {psf_dir}\n"
        f"\n"
        f"WRITE_XML        Y\n"
        f"XML_NAME         {Path(tmpdir) / 'psfex.xml'}\n"
        f"\n"
        f"CHECKPLOT_TYPE   NONE\n"
        f"CHECKIMAGE_TYPE  NONE\n"
        f"VERBOSE_TYPE     QUIET\n"
    )
    return str(psfex_path)


# ── Subprocess wrappers ─────────────────────────────────────────────────

def run_sextractor(filepath, telescope, tmpdir, fwhm_arcsec=None):
    """Run SExtractor on a FITS frame, producing LDAC catalogue.

    Returns path to LDAC catalogue, or None on failure.
    """
    ldac_path = str(Path(tmpdir) / "ldac.cat")
    sex_config = _write_sex_config(tmpdir, ldac_path, telescope, fwhm_arcsec)

    cmd = [config.SEXTRACTOR_CMD, str(filepath), "-c", sex_config]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except FileNotFoundError:
        log.error("SExtractor not found: %s", config.SEXTRACTOR_CMD)
        return None
    except subprocess.TimeoutExpired:
        log.error("SExtractor timed out on %s", filepath)
        return None

    if proc.returncode != 0:
        log.warning("SExtractor failed on %s: %s", filepath, proc.stderr.strip())
        return None

    if not Path(ldac_path).exists():
        log.warning("SExtractor produced no catalogue for %s", filepath)
        return None

    return ldac_path


def run_psfex(ldac_path, telescope, psf_output_dir, tmpdir):
    """Run PSFex on an LDAC catalogue, producing a .psf model.

    Returns dict with PSF stats parsed from XML, or None on failure.
    """
    psf_output_dir = Path(psf_output_dir)
    psf_output_dir.mkdir(parents=True, exist_ok=True)

    psfex_config = _write_psfex_config(tmpdir, str(psf_output_dir), telescope)
    xml_path = Path(tmpdir) / "psfex.xml"

    cmd = [config.PSFEX_CMD, ldac_path, "-c", psfex_config]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except FileNotFoundError:
        log.error("PSFex not found: %s", config.PSFEX_CMD)
        return None
    except subprocess.TimeoutExpired:
        log.error("PSFex timed out on %s", ldac_path)
        return None

    if proc.returncode != 0:
        log.warning("PSFex failed: %s", proc.stderr.strip())
        return None

    # Parse PSFex XML for statistics
    stats = _parse_psfex_xml(xml_path, telescope)
    return stats


def _parse_psfex_xml(xml_path, telescope):
    """Parse PSFex output XML for FWHM, ellipticity, chi2, star counts.

    Returns dict or None.
    """
    xml_path = Path(xml_path)
    if not xml_path.exists():
        log.warning("PSFex XML not found: %s", xml_path)
        return None

    try:
        tree = ET.parse(xml_path)
    except ET.ParseError as exc:
        log.warning("Failed to parse PSFex XML: %s", exc)
        return None

    root = tree.getroot()

    stats = {}

    # Find PSF_Fields TABLE anywhere in the tree (PSFex nests it under
    # RESOURCE ID="MetaData", not type="results")
    for table in root.iter("TABLE"):
        if table.get("name") != "PSF_Fields":
            continue

        # Collect FIELD names in order
        fields = [f.get("name", "") for f in table.iter("FIELD")]

        # Read first data row
        for tr in table.iter("TR"):
            tds = [td.text for td in tr.iter("TD")]
            if len(tds) != len(fields):
                continue
            row = dict(zip(fields, tds))
            try:
                if row.get("FWHM_Mean"):
                    stats["fwhm_px"] = float(row["FWHM_Mean"])
                if row.get("Ellipticity_Mean"):
                    stats["ellipticity"] = float(row["Ellipticity_Mean"])
                if row.get("Chi2_Mean"):
                    stats["chi2"] = float(row["Chi2_Mean"])
                if row.get("NStars_Loaded_Mean"):
                    stats["n_stars_loaded"] = int(float(row["NStars_Loaded_Mean"]))
                if row.get("NStars_Accepted_Mean"):
                    stats["n_stars_accepted"] = int(float(row["NStars_Accepted_Mean"]))
            except (ValueError, TypeError):
                pass
            break  # only need first row

    # Convert FWHM from pixels to arcseconds
    if "fwhm_px" in stats:
        pixel_scale = config.TELESCOPES[telescope]["pixel_scale"]
        stats["fwhm_arcsec"] = round(stats["fwhm_px"] * pixel_scale, 2)
        stats["fwhm_px"] = round(stats["fwhm_px"], 2)

    if "ellipticity" in stats:
        stats["ellipticity"] = round(stats["ellipticity"], 3)
    if "chi2" in stats:
        stats["chi2"] = round(stats["chi2"], 2)

    return stats if stats else None


# ── Per-frame orchestrator ──────────────────────────────────────────────

def process_frame(filepath, telescope, psf_dir, qc_stats=None, force=False):
    """Run SExtractor + PSFex on a single frame.

    Parameters
    ----------
    filepath : Path
        Reduced science frame.
    telescope : str
    psf_dir : Path
        Base output directory for PSF models.
    qc_stats : dict or None
        QC entry for this frame (for FWHM, n_good).
    force : bool
        Reprocess even if .psf model exists.

    Returns
    -------
    dict with status and PSF statistics.
    """
    filepath = Path(filepath)
    night = filepath.parent.name
    stem = filepath.stem

    psf_night_dir = psf_dir / night
    psf_model_path = psf_night_dir / f"{stem}.psf"

    hdr = fits.getheader(str(filepath))
    obj = hdr.get("OBJECT", "")
    filt = hdr.get("FILTER", "")

    result = {
        "filename": filepath.name,
        "night": night,
        "object": obj,
        "filter": filt,
    }

    # Skip if already done
    if not force and psf_model_path.exists():
        result["status"] = "SKIPPED"
        return result

    # Check minimum sources from QC
    if qc_stats and qc_stats.get("n_good", 0) < config.PSF_MIN_SOURCES:
        result["status"] = "FEW_SOURCES"
        result["n_good"] = qc_stats.get("n_good", 0)
        return result

    # Get FWHM estimate
    fwhm_arcsec = None
    if qc_stats and qc_stats.get("fwhm_arcsec") is not None:
        fwhm_arcsec = qc_stats["fwhm_arcsec"]

    with tempfile.TemporaryDirectory(prefix="psfex_") as tmpdir:
        # Step 1: SExtractor → LDAC
        ldac_path = run_sextractor(filepath, telescope, tmpdir,
                                   fwhm_arcsec=fwhm_arcsec)
        if ldac_path is None:
            result["status"] = "SEX_FAILED"
            return result

        # Step 2: PSFex → .psf model + XML stats
        psf_night_dir.mkdir(parents=True, exist_ok=True)
        stats = run_psfex(ldac_path, telescope, str(psf_night_dir), tmpdir)
        if stats is None:
            result["status"] = "PSFEX_FAILED"
            return result

        # Check if .psf file was created (PSFex names it from the LDAC filename)
        ldac_stem = Path(ldac_path).stem
        generated_psf = psf_night_dir / f"{ldac_stem}.psf"
        if generated_psf.exists():
            # Rename to match the science frame
            if generated_psf != psf_model_path:
                shutil.move(str(generated_psf), str(psf_model_path))
        elif not psf_model_path.exists():
            # PSFex may not have created a model
            result["status"] = "NO_PSF_MODEL"
            return result

    result["status"] = "OK"
    result.update(stats)
    return result


# ── Batch runner ────────────────────────────────────────────────────────

def run_psf(year, telescope, force=False):
    """Run PSF extraction on all reduced frames for a telescope/year.

    Writes psf/psf_report.json and psf_report.csv.
    """
    base = config.DATA_ROOT / year / telescope
    reduced_dir = base / "reduced"
    psf_dir = base / "psf"

    if not reduced_dir.exists():
        log.warning("No reduced directory: %s", reduced_dir)
        return

    # Check if output already exists
    json_out = psf_dir / "psf_report.json"
    if json_out.exists() and not force:
        log.info("PSF report already exists: %s (use --force to rerun)", json_out)
        print(f"  PSF report already exists for {telescope} — use --force to rerun")
        return

    # Load QC frame stats
    frame_stats_file = base / "qc" / "frame_stats.json"
    qc_lookup = {}
    if frame_stats_file.exists():
        with open(frame_stats_file) as f:
            for entry in json.load(f):
                qc_lookup[entry["filename"]] = entry

    # Collect all reduced FITS files
    fits_files = sorted(reduced_dir.rglob("*.fits")) + sorted(reduced_dir.rglob("*.fit"))
    if not fits_files:
        log.warning("No reduced FITS files in %s", reduced_dir)
        return

    print(f"\n  {telescope}: PSF extraction for {len(fits_files)} frames")

    results = []
    for i, fpath in enumerate(fits_files, 1):
        qc = qc_lookup.get(fpath.name)
        try:
            result = process_frame(fpath, telescope, psf_dir,
                                   qc_stats=qc, force=force)
        except Exception as exc:
            log.error("PSF extraction failed for %s: %s", fpath.name, exc)
            result = {
                "filename": fpath.name,
                "night": fpath.parent.name,
                "status": f"ERROR: {exc}",
            }
        results.append(result)

        if i % 20 == 0:
            n_ok = sum(1 for r in results if r.get("status") == "OK")
            print(f"    [{i}/{len(fits_files)}] {n_ok} models built so far...")

    # Write reports
    psf_dir.mkdir(parents=True, exist_ok=True)

    with open(json_out, "w") as f:
        json.dump(results, f, indent=2, default=str)
    log.info("Wrote %s (%d frames)", json_out, len(results))

    csv_out = psf_dir / "psf_report.csv"
    if results:
        fieldnames = list(results[0].keys())
        for r in results:
            for k in r:
                if k not in fieldnames:
                    fieldnames.append(k)
        with open(csv_out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow({k: row.get(k, "") for k in fieldnames})
        log.info("Wrote %s", csv_out)

    # Console summary
    _print_summary(results, telescope)


def _print_summary(results, telescope):
    """Print PSF extraction summary."""
    n_ok = sum(1 for r in results if r.get("status") == "OK")
    n_total = len(results)
    status_counts = {}
    for r in results:
        s = r.get("status", "?")
        status_counts[s] = status_counts.get(s, 0) + 1

    print(f"\n  ── {telescope} PSF Summary: {n_ok}/{n_total} models built ──")

    ok_results = [r for r in results if r.get("status") == "OK"]
    if ok_results:
        fwhms = [r["fwhm_arcsec"] for r in ok_results
                 if r.get("fwhm_arcsec") is not None]
        ellips = [r["ellipticity"] for r in ok_results
                  if r.get("ellipticity") is not None]
        chi2s = [r["chi2"] for r in ok_results
                 if r.get("chi2") is not None]
        n_stars = [r["n_stars_accepted"] for r in ok_results
                   if r.get("n_stars_accepted") is not None]

        if fwhms:
            print(f"  FWHM:  median={np.median(fwhms):.2f}\"  "
                  f"range=[{min(fwhms):.2f}\", {max(fwhms):.2f}\"]")
        if ellips:
            print(f"  Ellip: median={np.median(ellips):.3f}  "
                  f"range=[{min(ellips):.3f}, {max(ellips):.3f}]")
        if chi2s:
            print(f"  Chi2:  median={np.median(chi2s):.2f}")
        if n_stars:
            print(f"  Stars: median={int(np.median(n_stars))}  "
                  f"range=[{min(n_stars)}, {max(n_stars)}]")

    # Status breakdown
    non_ok = {k: v for k, v in status_counts.items() if k != "OK"}
    if non_ok:
        print(f"\n  Status breakdown:")
        for s, cnt in sorted(non_ok.items(), key=lambda x: -x[1]):
            print(f"    {s}: {cnt}")
    print()
