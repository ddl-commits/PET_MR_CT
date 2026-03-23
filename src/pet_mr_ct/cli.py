"""CLI entry point for the PET_MR_CT pipeline."""

from __future__ import annotations

from pathlib import Path

import click


@click.group()
@click.version_option()
def main() -> None:
    """DICOM-SEG pipeline: AI segmentation to MIM/3D Slicer round-trip."""


@main.group()
def convert() -> None:
    """Convert between DICOM, NIfTI, DICOM-SEG, and RT-STRUCT formats."""


@convert.command()
@click.option(
    "--source-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Directory containing source DICOM series.",
)
@click.option(
    "--mask-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Directory containing NIfTI segmentation masks (.nii.gz or .nii).",
)
@click.option(
    "--tool",
    required=True,
    type=str,
    help="Segmentation tool name -- must match a configs/segment_metadata_{tool}.json file (e.g. totalseg, moose).",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output DICOM-SEG file path.",
)
def dcmseg(source_dir: Path, mask_dir: Path, tool: str, output: Path) -> None:
    """Create DICOM-SEG from NIfTI masks + source DICOM series."""
    from pet_mr_ct.convert.nifti_to_dcmseg import create_dcmseg

    # Discover NIfTI masks
    masks: dict[str, Path] = {}
    for ext in ("*.nii.gz", "*.nii"):
        for f in mask_dir.glob(ext):
            stem = f.name.replace(".nii.gz", "").replace(".nii", "")
            if stem not in masks:  # .nii.gz takes priority
                masks[stem] = f

    if not masks:
        raise click.ClickException(f"No .nii.gz or .nii files found in {mask_dir}")

    click.echo(f"Found {len(masks)} mask(s) in {mask_dir}")
    result = create_dcmseg(source_dir, masks, tool, output)
    click.echo(f"DICOM-SEG written to {result}")


@convert.command()
@click.option("--source-dir", required=True, type=click.Path(path_type=Path))
@click.option("--mask-dir", required=True, type=click.Path(path_type=Path))
@click.option("--output", required=True, type=click.Path(path_type=Path))
def rtstruct(source_dir: Path, mask_dir: Path, output: Path) -> None:
    """Create RT-STRUCT from NIfTI masks (not yet implemented)."""
    click.echo("RT-STRUCT creation is not yet implemented.")
    click.echo("Use 'pet-mr-ct convert dcmseg' for DICOM-SEG output.")
    raise SystemExit(1)


if __name__ == "__main__":
    main()
