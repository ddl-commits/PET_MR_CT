"""Convert NIfTI segmentation masks to DICOM-SEG via highdicom."""

from __future__ import annotations

import logging
from pathlib import Path

import highdicom as hd
import nibabel as nib
import numpy as np
import pydicom

from pet_mr_ct.convert.terminology import get_segment_description

logger = logging.getLogger(__name__)


def _load_dicom_series(dcm_dir: Path) -> list[pydicom.Dataset]:
    """Load and z-sort a DICOM series from a directory.

    Tries .dcm files first; falls back to all non-hidden files.
    Catches InvalidDicomError and skips non-DICOM files.
    """
    dcm_files = sorted(dcm_dir.glob("*.dcm"))
    if not dcm_files:
        dcm_files = [
            f for f in sorted(dcm_dir.iterdir())
            if f.is_file() and not f.name.startswith(".")
        ]

    if not dcm_files:
        raise FileNotFoundError(f"No files found in {dcm_dir}")

    datasets: list[pydicom.Dataset] = []
    for f in dcm_files:
        try:
            ds = pydicom.dcmread(f)
            datasets.append(ds)
        except pydicom.errors.InvalidDicomError:
            logger.warning("Skipping non-DICOM file: %s", f.name)

    if not datasets:
        raise FileNotFoundError(f"No valid DICOM files in {dcm_dir}")

    datasets.sort(key=lambda ds: float(ds.ImagePositionPatient[2]))
    return datasets


def _load_nifti_mask(
    mask_path: Path,
    num_slices: int,
    rows: int,
    cols: int,
) -> np.ndarray | None:
    """Load a NIfTI mask and convert from RAS+ to DICOM (slice, row, col).

    Returns a boolean array of shape (num_slices, rows, cols), or None if
    the mask cannot be loaded or dimensions cannot be reconciled.
    """
    nii = nib.load(mask_path)

    # Reorient to closest canonical (RAS+), then convert to LPS
    canonical = nib.as_closest_canonical(nii)
    mask_ras = np.asarray(canonical.dataobj, dtype=np.float32)

    # RAS+ -> LPS+: flip R->L and A->P
    mask_lps = mask_ras[::-1, ::-1, :]

    # (L, P, S) -> (S, P, L) = (slice, row, col)
    mask_dicom = np.transpose(mask_lps, (2, 1, 0))

    # Verify dimensions
    if mask_dicom.shape != (num_slices, rows, cols):
        if mask_dicom.shape == (num_slices, cols, rows):
            logger.warning(
                "Swapping rows/cols for %s: got %s, expected (%d, %d, %d)",
                mask_path.name, mask_dicom.shape, num_slices, rows, cols,
            )
            mask_dicom = np.swapaxes(mask_dicom, 1, 2)
        else:
            logger.error(
                "Dimension mismatch for %s: got %s, expected (%d, %d, %d) -- skipping",
                mask_path.name, mask_dicom.shape, num_slices, rows, cols,
            )
            return None

    # Binarize
    mask_binary = mask_dicom > 0.5

    # Handle z-axis inversion
    affine = canonical.affine
    if affine[2, 2] < 0:
        mask_binary = mask_binary[::-1, :, :]

    return mask_binary


def create_dcmseg(
    source_dcm_dir: Path,
    nifti_masks: dict[str, Path],
    tool: str,
    output_path: Path,
) -> Path:
    """Create a DICOM-SEG file from NIfTI segmentation masks.

    Args:
        source_dcm_dir: Directory containing the source DICOM series.
        nifti_masks: Mapping of label name to NIfTI file path.
        tool: Segmentation tool name for terminology lookup.
        output_path: Where to write the DICOM-SEG file.

    Returns:
        Path to the written DICOM-SEG file.

    Raises:
        ValueError: If nifti_masks is empty or no valid masks found.
        FileNotFoundError: If source_dcm_dir has no valid DICOM files.
    """
    if not nifti_masks:
        raise ValueError("No masks provided")

    # Load source DICOM series
    datasets = _load_dicom_series(source_dcm_dir)
    num_slices = len(datasets)
    rows = datasets[0].Rows
    cols = datasets[0].Columns

    # Load and convert each NIfTI mask
    valid_labels: list[str] = []
    mask_arrays: list[np.ndarray] = []

    for label, mask_path in nifti_masks.items():
        if not mask_path.exists():
            logger.warning("Mask file not found, skipping: %s", mask_path)
            continue

        mask = _load_nifti_mask(mask_path, num_slices, rows, cols)
        if mask is None:
            continue

        valid_labels.append(label)
        mask_arrays.append(mask)

    if not valid_labels:
        raise ValueError("No valid masks found")

    # Combine into 4D array: (slices, rows, cols, num_segments)
    pixel_array = np.stack(mask_arrays, axis=-1)

    # Build segment descriptions
    segment_descriptions = [
        get_segment_description(label, tool, seg_num)
        for seg_num, label in enumerate(valid_labels, start=1)
    ]

    # Create DICOM-SEG
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seg = hd.seg.Segmentation(
        source_images=datasets,
        pixel_array=pixel_array,
        segmentation_type=hd.seg.SegmentationTypeValues.BINARY,
        segment_descriptions=segment_descriptions,
        series_instance_uid=hd.UID(),
        series_number=100,
        sop_instance_uid=hd.UID(),
        instance_number=1,
        manufacturer=getattr(datasets[0], "Manufacturer", "Unknown"),
        manufacturer_model_name=getattr(
            datasets[0], "ManufacturerModelName", "Unknown"
        ),
        software_versions=f"pet-mr-ct 0.1.0 ({tool})",
        device_serial_number=getattr(datasets[0], "DeviceSerialNumber", "0000"),
        content_description=f"{tool} segmentation",
        content_label=tool.upper(),
        omit_empty_frames=True,
    )

    seg.save_as(output_path)
    logger.info(
        "Created DICOM-SEG with %d segments: %s",
        len(valid_labels), output_path,
    )
    return output_path
