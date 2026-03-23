"""Convert DICOM series to NIfTI format."""

from __future__ import annotations

from pathlib import Path


def convert_dcm_to_nifti(dcm_dir: Path, output_path: Path) -> Path:
    """Convert a DICOM series to NIfTI format.

    Args:
        dcm_dir: Directory containing DICOM files.
        output_path: Where to write the NIfTI file.

    Returns:
        Path to the written NIfTI file.

    Raises:
        NotImplementedError: Not yet implemented.
    """
    raise NotImplementedError("DICOM-to-NIfTI conversion not yet implemented.")
