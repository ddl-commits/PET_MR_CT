"""Convert NIfTI segmentation masks to RT-STRUCT via rt-utils."""

from __future__ import annotations

from pathlib import Path


def create_rtstruct(
    source_dcm_dir: Path,
    nifti_masks: dict[str, Path],
    output_path: Path,
) -> Path:
    """Create an RT-STRUCT from NIfTI segmentation masks via rt-utils.

    Args:
        source_dcm_dir: Directory containing the source DICOM series.
        nifti_masks: Mapping of label name to NIfTI file path.
        output_path: Where to write the RT-STRUCT file.

    Returns:
        Path to the written RT-STRUCT file.

    Raises:
        NotImplementedError: Not yet implemented.
    """
    raise NotImplementedError("RT-STRUCT creation not yet implemented.")
