"""CLI entry point for the PET_MR_CT pipeline."""

import click


@click.group()
@click.version_option()
def main() -> None:
    """DICOM-SEG pipeline: AI segmentation to MIM/3D Slicer round-trip."""


if __name__ == "__main__":
    main()
