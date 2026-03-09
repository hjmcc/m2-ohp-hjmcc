#!/usr/bin/env python3
"""Extract PSF models from stacked images using SExtractor + PSFex.

For each stack, runs SExtractor → LDAC → PSFex → .psf model, then
extracts the centre-of-field PSF into a simple 2D FITS image that
students can use directly.

Output per stack:
    {target}_{filter}_psf.fits   — 2D PSF image (25x25 px, sum=1)

Also writes stack_psf_report.json to the stacks root.

Usage:
    python scripts/extract_stack_psfs.py
    python scripts/extract_stack_psfs.py --target M67
    python scripts/extract_stack_psfs.py --force
    python scripts/extract_stack_psfs.py --dry-run
"""

import argparse
import json
import logging
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pipeline import config
from pipeline.psf import run_sextractor, run_psfex, _parse_psfex_xml

logging.basicConfig(level=logging.WARNING, format="%(message)s")
log = logging.getLogger("pipeline")
log.setLevel(logging.INFO)

PIXEL_SCALE = 0.770  # arcsec/px (stack pixel scale)
TELESCOPE = "T120"


def extract_psf_image(psf_path):
    """Read a PSFex .psf file and return the centre-of-field PSF as a 2D array.

    Returns (psf_2d, n_basis) or (None, 0) on failure.
    """
    try:
        with fits.open(str(psf_path)) as hdul:
            if len(hdul) < 2 or "PSF_MASK" not in hdul[1].columns.names:
                return None, 0
            mask = hdul[1].data["PSF_MASK"][0]  # (n_basis, ny, nx)
            psf_centre = mask[0].copy()
            s = psf_centre.sum()
            if s > 0:
                psf_centre /= s
            return psf_centre, mask.shape[0]
    except Exception as e:
        log.warning("Failed to read PSF model %s: %s", psf_path, e)
        return None, 0


def process_stack(stack_path, weight_path, output_psf_path, force=False):
    """Run SExtractor + PSFex on a single stack image.

    Returns dict with results.
    """
    stack_path = Path(stack_path)
    output_psf_path = Path(output_psf_path)
    stem = stack_path.stem  # e.g. M67_r

    result = {"file": stack_path.name, "psf_file": output_psf_path.name}

    # Skip if already done
    if not force and output_psf_path.exists():
        result["status"] = "SKIPPED"
        # Read existing PSF for stats
        hdr = fits.getheader(str(output_psf_path))
        result["fwhm_px"] = hdr.get("FWHM_PX")
        result["fwhm_as"] = hdr.get("FWHM_AS")
        result["ellipticity"] = hdr.get("ELLIPT")
        result["n_stars"] = hdr.get("NSTARS")
        return result

    # Get FWHM estimate from stack header
    hdr = fits.getheader(str(stack_path))
    fwhm_as = hdr.get("FWHM_AS", 3.5)

    with tempfile.TemporaryDirectory(prefix="stack_psf_") as tmpdir:
        tmpdir = Path(tmpdir)

        # SExtractor needs a single-extension FITS (stacks already are)
        # Run SExtractor → LDAC
        ldac_path = str(tmpdir / "ldac.cat")

        # Write custom SExtractor config for stacks
        # Use weight map if available
        param_path = tmpdir / "default.param"
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
            f"VIGNET({config.PSF_VIGNET_SIZE},{config.PSF_VIGNET_SIZE})\n"
        )

        aper_px = fwhm_as / PIXEL_SCALE * 3.0

        weight_opts = ""
        if weight_path and Path(weight_path).exists():
            weight_opts = (
                f"WEIGHT_TYPE      MAP_WEIGHT\n"
                f"WEIGHT_IMAGE     {weight_path}\n"
            )
        else:
            weight_opts = "WEIGHT_TYPE      NONE\n"

        sex_path = tmpdir / "default.sex"
        sex_path.write_text(
            f"CATALOG_NAME     {ldac_path}\n"
            f"CATALOG_TYPE     FITS_LDAC\n"
            f"PARAMETERS_NAME  {param_path}\n"
            f"\n"
            f"DETECT_TYPE      CCD\n"
            f"DETECT_THRESH    5.0\n"
            f"ANALYSIS_THRESH  5.0\n"
            f"\n"
            f"FILTER           N\n"
            f"\n"
            f"DEBLEND_NTHRESH  32\n"
            f"DEBLEND_MINCONT  0.005\n"
            f"\n"
            f"CLEAN            Y\n"
            f"CLEAN_PARAM      1.0\n"
            f"\n"
            f"PHOT_APERTURES   {aper_px:.1f}\n"
            f"\n"
            f"SATUR_LEVEL      1e10\n"
            f"PIXEL_SCALE      {PIXEL_SCALE}\n"
            f"SEEING_FWHM      {fwhm_as:.1f}\n"
            f"\n"
            f"BACK_SIZE        128\n"
            f"BACK_FILTERSIZE  3\n"
            f"\n"
            f"{weight_opts}"
            f"CHECKIMAGE_TYPE  NONE\n"
            f"VERBOSE_TYPE     QUIET\n"
        )

        import subprocess
        proc = subprocess.run(
            [config.SEXTRACTOR_CMD, str(stack_path), "-c", str(sex_path)],
            capture_output=True, text=True, timeout=120,
        )
        if proc.returncode != 0:
            result["status"] = "SEX_FAILED"
            result["error"] = proc.stderr.strip()[-200:]
            return result

        if not Path(ldac_path).exists():
            result["status"] = "SEX_NO_CAT"
            return result

        # Run PSFex
        psf_outdir = tmpdir / "psf_out"
        psf_outdir.mkdir()

        psfex_config = tmpdir / "default.psfex"
        psfex_config.write_text(
            f"BASIS_TYPE       PIXEL_AUTO\n"
            f"PSF_SIZE         25,25\n"
            f"PSF_SAMPLING     0.0\n"
            f"\n"
            f"PSFVAR_KEYS      X_IMAGE,Y_IMAGE\n"
            f"PSFVAR_GROUPS    1,1\n"
            f"PSFVAR_DEGREES   2\n"
            f"\n"
            f"SAMPLE_AUTOSELECT Y\n"
            f"SAMPLE_FWHMRANGE  2.0,15.0\n"
            f"SAMPLE_VARIABILITY 0.3\n"
            f"SAMPLE_MINSN     10\n"
            f"SAMPLE_MAXELLIP  0.3\n"
            f"\n"
            f"PSF_DIR          {psf_outdir}\n"
            f"\n"
            f"WRITE_XML        Y\n"
            f"XML_NAME         {tmpdir / 'psfex.xml'}\n"
            f"\n"
            f"CHECKPLOT_TYPE   NONE\n"
            f"CHECKIMAGE_TYPE  NONE\n"
            f"VERBOSE_TYPE     QUIET\n"
        )

        proc = subprocess.run(
            [config.PSFEX_CMD, ldac_path, "-c", str(psfex_config)],
            capture_output=True, text=True, timeout=120,
        )
        if proc.returncode != 0:
            result["status"] = "PSFEX_FAILED"
            result["error"] = proc.stderr.strip()[-200:]
            return result

        # Parse XML stats
        stats = _parse_psfex_xml(tmpdir / "psfex.xml", TELESCOPE)

        # Find the .psf file PSFex produced
        psf_files = list(psf_outdir.glob("*.psf"))
        if not psf_files:
            result["status"] = "NO_PSF_MODEL"
            return result

        # Extract centre-of-field PSF
        psf_2d, n_basis = extract_psf_image(psf_files[0])
        if psf_2d is None:
            result["status"] = "PSF_READ_FAILED"
            return result

        # Compute encircled energy radii
        cy, cx = psf_2d.shape[0] / 2, psf_2d.shape[1] / 2
        Y, X = np.ogrid[:psf_2d.shape[0], :psf_2d.shape[1]]
        r = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
        r_flat = r.ravel()
        psf_flat = psf_2d.ravel()
        order = np.argsort(r_flat)
        r_sorted = r_flat[order]
        ee = np.cumsum(psf_flat[order])
        r50 = float(np.interp(0.5, ee, r_sorted))
        r80 = float(np.interp(0.8, ee, r_sorted))

        # Get stats from PSFex XML
        fwhm_px = stats.get("fwhm_px", 2 * r50) if stats else 2 * r50
        fwhm_as = round(fwhm_px * PIXEL_SCALE, 2)
        ellipticity = stats.get("ellipticity", 0) if stats else 0
        chi2 = stats.get("chi2", 0) if stats else 0
        n_stars_loaded = stats.get("n_stars_loaded", 0) if stats else 0
        n_stars_accepted = stats.get("n_stars_accepted", 0) if stats else 0

        # Write output PSF FITS
        out_hdr = fits.Header()
        out_hdr["OBJECT"] = hdr.get("OBJECT", "")
        out_hdr["FILTER"] = hdr.get("FILTER", "")
        out_hdr["PIXSCALE"] = (PIXEL_SCALE, "arcsec/px")
        out_hdr["FWHM_PX"] = (round(fwhm_px, 2), "PSF FWHM in pixels")
        out_hdr["FWHM_AS"] = (fwhm_as, "PSF FWHM in arcsec")
        out_hdr["EE50_AS"] = (round(r50 * PIXEL_SCALE, 2),
                              "50% encircled energy radius (arcsec)")
        out_hdr["EE80_AS"] = (round(r80 * PIXEL_SCALE, 2),
                              "80% encircled energy radius (arcsec)")
        out_hdr["ELLIPT"] = (round(ellipticity, 3), "Mean PSF ellipticity")
        out_hdr["CHI2"] = (round(chi2, 2), "PSFex chi2")
        out_hdr["NSTARS"] = (n_stars_accepted, "Stars used for PSF model")
        out_hdr["NLOADED"] = (n_stars_loaded, "Stars loaded by PSFex")
        out_hdr["NBASIS"] = (n_basis, "Number of PSFex basis images")
        out_hdr["PSFVAR"] = ("degree-2 polynomial",
                             "Spatial variation model")
        out_hdr.add_comment("PSF extracted from stack via SExtractor + PSFex")
        out_hdr.add_comment("Centre-of-field component, normalised to sum=1")

        fits.writeto(str(output_psf_path), psf_2d.astype(np.float32),
                     out_hdr, overwrite=True)

        result["status"] = "OK"
        result["fwhm_px"] = round(fwhm_px, 2)
        result["fwhm_as"] = fwhm_as
        result["ee50_as"] = round(r50 * PIXEL_SCALE, 2)
        result["ee80_as"] = round(r80 * PIXEL_SCALE, 2)
        result["ellipticity"] = round(ellipticity, 3)
        result["chi2"] = round(chi2, 2)
        result["n_stars"] = n_stars_accepted

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Extract PSF models from stacked images")
    parser.add_argument("--target", type=str, default=None,
                        help="Process only this target")
    parser.add_argument("--force", action="store_true",
                        help="Reprocess even if PSF exists")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show inventory only")
    args = parser.parse_args()

    stacks_root = config.DATA_ROOT / "stacks"

    # Find all stacks
    stacks = []
    for target_dir in sorted(stacks_root.iterdir()):
        if not target_dir.is_dir():
            continue
        target = target_dir.name
        if args.target and target.lower() != args.target.lower():
            continue

        t120_dir = target_dir / "T120"
        if not t120_dir.is_dir():
            continue

        for fits_file in sorted(t120_dir.glob("*.fits")):
            if "weight" in fits_file.name or "psf" in fits_file.name:
                continue
            weight_file = fits_file.with_name(
                fits_file.stem + ".weight.fits")
            psf_output = fits_file.with_name(
                fits_file.stem + "_psf.fits")
            stacks.append({
                "stack": fits_file,
                "weight": weight_file if weight_file.exists() else None,
                "psf_output": psf_output,
                "target": target,
            })

    print(f"\n{'=' * 60}")
    print(f"  Stack PSF extraction — {len(stacks)} stacks")
    print(f"{'=' * 60}\n")

    if args.dry_run:
        for s in stacks:
            exists = "EXISTS" if s["psf_output"].exists() else "missing"
            print(f"  {s['stack'].name:30s}  → {s['psf_output'].name:30s}  [{exists}]")
        return

    results = []
    for i, s in enumerate(stacks, 1):
        print(f"  [{i:3d}/{len(stacks)}] {s['stack'].name:30s}", end="",
              flush=True)

        result = process_stack(
            s["stack"],
            s["weight"],
            s["psf_output"],
            force=args.force,
        )
        results.append(result)

        status = result["status"]
        if status == "OK":
            print(f"  FWHM={result['fwhm_as']:.2f}\"  "
                  f"e={result['ellipticity']:.3f}  "
                  f"n={result['n_stars']}")
        elif status == "SKIPPED":
            print(f"  skipped (exists)")
        else:
            print(f"  {status}")

    # Summary
    n_ok = sum(1 for r in results if r["status"] == "OK")
    n_skip = sum(1 for r in results if r["status"] == "SKIPPED")
    n_fail = sum(1 for r in results
                 if r["status"] not in ("OK", "SKIPPED"))

    print(f"\n{'=' * 60}")
    print(f"  {n_ok} new PSFs extracted, {n_skip} skipped, {n_fail} failed")

    ok_results = [r for r in results if r["status"] in ("OK", "SKIPPED")
                  and r.get("fwhm_as")]
    if ok_results:
        fwhms = [r["fwhm_as"] for r in ok_results]
        print(f"  FWHM: median={np.median(fwhms):.2f}\"  "
              f"range=[{min(fwhms):.2f}\", {max(fwhms):.2f}\"]")
    print(f"{'=' * 60}")

    # Save report
    report_path = stacks_root / "stack_psf_report.json"
    with open(report_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Report: {report_path}")


if __name__ == "__main__":
    main()
