"""Tests for the convert module."""

from __future__ import annotations

from pathlib import Path

import highdicom as hd
import pytest

from pet_mr_ct.convert.nifti_to_dcmseg import create_dcmseg


def test_create_dcmseg_produces_file(
    tmp_dicom_series: Path,
    tmp_nifti_masks: Path,
    tmp_path: Path,
) -> None:
    """create_dcmseg writes a non-empty DICOM-SEG file."""
    output = tmp_path / "output" / "test.dcm"
    masks = {
        "liver": tmp_nifti_masks / "liver.nii.gz",
        "spleen": tmp_nifti_masks / "spleen.nii.gz",
        "kidney_right": tmp_nifti_masks / "kidney_right.nii.gz",
    }
    result = create_dcmseg(tmp_dicom_series, masks, "totalseg", output)
    assert result == output
    assert output.exists()
    assert output.stat().st_size > 0


def test_dcmseg_has_correct_segment_count(
    tmp_dicom_series: Path,
    tmp_nifti_masks: Path,
    tmp_path: Path,
) -> None:
    """DICOM-SEG contains exactly 3 segments."""
    output = tmp_path / "output" / "test.dcm"
    masks = {
        "liver": tmp_nifti_masks / "liver.nii.gz",
        "spleen": tmp_nifti_masks / "spleen.nii.gz",
        "kidney_right": tmp_nifti_masks / "kidney_right.nii.gz",
    }
    create_dcmseg(tmp_dicom_series, masks, "totalseg", output)
    seg = hd.seg.segread(str(output))
    assert seg.number_of_segments == 3


def test_dcmseg_segment_labels(
    tmp_dicom_series: Path,
    tmp_nifti_masks: Path,
    tmp_path: Path,
) -> None:
    """DICOM-SEG segment labels match input label display names."""
    output = tmp_path / "output" / "test.dcm"
    masks = {
        "liver": tmp_nifti_masks / "liver.nii.gz",
        "spleen": tmp_nifti_masks / "spleen.nii.gz",
        "kidney_right": tmp_nifti_masks / "kidney_right.nii.gz",
    }
    create_dcmseg(tmp_dicom_series, masks, "totalseg", output)
    seg = hd.seg.segread(str(output))
    labels = {
        seg.get_segment_description(i + 1).segment_label
        for i in range(seg.number_of_segments)
    }
    assert "Liver" in labels
    assert "Spleen" in labels
    assert "Right Kidney" in labels


def test_create_dcmseg_empty_masks_raises(
    tmp_dicom_series: Path,
    tmp_path: Path,
) -> None:
    """Empty masks dict raises ValueError."""
    output = tmp_path / "output" / "test.dcm"
    with pytest.raises(ValueError, match="No masks provided"):
        create_dcmseg(tmp_dicom_series, {}, "totalseg", output)


def test_create_dcmseg_missing_mask_skipped(
    tmp_dicom_series: Path,
    tmp_nifti_masks: Path,
    tmp_path: Path,
) -> None:
    """Missing mask files are skipped with warning, valid ones still processed."""
    output = tmp_path / "output" / "test.dcm"
    masks = {
        "liver": tmp_nifti_masks / "liver.nii.gz",
        "nonexistent": tmp_nifti_masks / "nonexistent.nii.gz",
    }
    create_dcmseg(tmp_dicom_series, masks, "totalseg", output)
    seg = hd.seg.segread(str(output))
    assert seg.number_of_segments == 1
