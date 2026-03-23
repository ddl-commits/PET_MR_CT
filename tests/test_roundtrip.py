"""Tests for round-trip validation (export -> reimport -> compare)."""

from __future__ import annotations

from pathlib import Path

import highdicom as hd
import nibabel as nib
import numpy as np
import pydicom

from pet_mr_ct.convert.nifti_to_dcmseg import create_dcmseg


def test_voxel_agreement(
    tmp_dicom_series: Path,
    tmp_nifti_masks: Path,
    tmp_path: Path,
) -> None:
    """DICOM-SEG round-trip preserves >99% voxel agreement per segment."""
    output = tmp_path / "output" / "roundtrip.dcm"
    mask_names = ["liver", "spleen", "kidney_right"]
    masks = {name: tmp_nifti_masks / f"{name}.nii.gz" for name in mask_names}

    create_dcmseg(tmp_dicom_series, masks, "totalseg", output)

    # Read back the DICOM-SEG
    seg = hd.seg.segread(str(output))

    # Collect source SOP Instance UIDs in z-sorted order (matching create_dcmseg)
    dcm_files = sorted(tmp_dicom_series.glob("*.dcm"))
    source_datasets = [pydicom.dcmread(f) for f in dcm_files]
    source_datasets.sort(key=lambda ds: float(ds.ImagePositionPatient[2]))
    source_uids = [ds.SOPInstanceUID for ds in source_datasets]

    for seg_num, name in enumerate(mask_names, start=1):
        # Load original NIfTI and convert to DICOM orientation (same logic as
        # nifti_to_dcmseg._load_nifti_mask)
        nii = nib.load(masks[name])
        canonical = nib.as_closest_canonical(nii)
        original_ras = np.asarray(canonical.dataobj, dtype=np.float32)
        original_lps = original_ras[::-1, ::-1, :]
        original_dicom = np.transpose(original_lps, (2, 1, 0))

        if canonical.affine[2, 2] < 0:
            original_dicom = original_dicom[::-1, :, :]

        original_mask = original_dicom > 0.5

        # Extract segment from DICOM-SEG using source instance UIDs
        seg_array = seg.get_pixels_by_source_instance(
            source_sop_instance_uids=source_uids,
            segment_numbers=[seg_num],
            assert_missing_frames_are_empty=True,
        )
        # Shape: (num_instances, rows, cols, 1) -> squeeze last dim
        seg_mask = seg_array.squeeze(axis=-1).astype(bool)

        assert seg_mask.shape == original_mask.shape, (
            f"Shape mismatch for {name}: "
            f"seg={seg_mask.shape}, original={original_mask.shape}"
        )

        total_voxels = original_mask.size
        matching = int(np.sum(seg_mask == original_mask))
        agreement = matching / total_voxels

        assert agreement > 0.99, (
            f"Voxel agreement for {name}: {agreement:.4f} "
            f"({matching}/{total_voxels})"
        )
