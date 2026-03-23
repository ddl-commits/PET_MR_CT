"""Shared pytest fixtures for PET_MR_CT tests."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pydicom
import pytest
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

# Small fixture geometry for fast tests
ROWS = 64
COLS = 64
NUM_SLICES = 8
PIXEL_SPACING = [1.0, 1.0]
SLICE_SPACING = 2.5

# Shared UIDs for consistent geometry
STUDY_UID = generate_uid()
SERIES_UID = generate_uid()
FRAME_OF_REF_UID = generate_uid()

# Test organs: (name, center_zyx_frac, radii_zyx_frac)
TEST_ORGANS = [
    ("liver", (0.5, 0.5, 0.6), (0.3, 0.3, 0.2)),
    ("spleen", (0.5, 0.5, 0.3), (0.2, 0.2, 0.15)),
    ("kidney_right", (0.5, 0.6, 0.7), (0.15, 0.15, 0.1)),
]


def _make_ellipsoid(
    shape: tuple[int, int, int],
    center_frac: tuple[float, float, float],
    radii_frac: tuple[float, float, float],
) -> np.ndarray:
    """Create a binary ellipsoid mask."""
    zz, yy, xx = np.ogrid[0 : shape[0], 0 : shape[1], 0 : shape[2]]
    cz, cy, cx = [c * s for c, s in zip(center_frac, shape)]
    rz, ry, rx = [r * s for r, s in zip(radii_frac, shape)]
    dist = ((zz - cz) / rz) ** 2 + ((yy - cy) / ry) ** 2 + ((xx - cx) / rx) ** 2
    return (dist <= 1.0).astype(np.uint8)


@pytest.fixture()
def tmp_dicom_series(tmp_path: Path) -> Path:
    """Generate an 8-slice synthetic CT DICOM series.

    Returns the directory containing the DICOM files.
    """
    dcm_dir = tmp_path / "dicom_ct"
    dcm_dir.mkdir()

    for i in range(NUM_SLICES):
        sop_uid = generate_uid()
        z_pos = i * SLICE_SPACING

        img = np.zeros((ROWS, COLS), dtype=np.int16)

        file_meta = pydicom.Dataset()
        file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        file_meta.MediaStorageSOPInstanceUID = sop_uid
        file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

        filename = dcm_dir / f"CT_{i:04d}.dcm"
        ds = pydicom.dataset.FileDataset(
            str(filename), {}, file_meta=file_meta, preamble=b"\x00" * 128,
        )

        ds.PatientName = "Test^Phantom"
        ds.PatientID = "TEST001"
        ds.PatientBirthDate = "19700101"
        ds.PatientSex = "O"

        ds.StudyInstanceUID = STUDY_UID
        ds.StudyDate = "20260322"
        ds.StudyTime = "120000"
        ds.StudyDescription = "Test"
        ds.StudyID = "1"
        ds.AccessionNumber = ""

        ds.SeriesInstanceUID = SERIES_UID
        ds.SeriesNumber = 1
        ds.SeriesDescription = "Test CT"
        ds.Modality = "CT"

        ds.FrameOfReferenceUID = FRAME_OF_REF_UID
        ds.PositionReferenceIndicator = ""

        ds.Manufacturer = "TestManufacturer"
        ds.InstitutionName = "TestLab"

        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        ds.SOPInstanceUID = sop_uid
        ds.InstanceNumber = i + 1
        ds.ImagePositionPatient = [0.0, 0.0, z_pos]
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        ds.PixelSpacing = PIXEL_SPACING
        ds.SliceThickness = SLICE_SPACING
        ds.SpacingBetweenSlices = SLICE_SPACING
        ds.SliceLocation = z_pos

        ds.Rows = ROWS
        ds.Columns = COLS
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.RescaleIntercept = 0
        ds.RescaleSlope = 1

        ds.PixelData = img.tobytes()
        ds.save_as(filename)

    return dcm_dir


@pytest.fixture()
def tmp_nifti_masks(tmp_path: Path) -> Path:
    """Generate 3 binary NIfTI masks matching tmp_dicom_series geometry.

    Returns the directory containing the .nii.gz files.
    """
    mask_dir = tmp_path / "segmentations"
    mask_dir.mkdir()

    volume_shape = (NUM_SLICES, ROWS, COLS)  # z, y, x

    for name, center, radii in TEST_ORGANS:
        mask_zyx = _make_ellipsoid(volume_shape, center, radii)

        # Transpose to (x, y, z) for NIfTI RAS+ convention
        mask_xyz = np.transpose(mask_zyx, (2, 1, 0))

        affine = np.diag([PIXEL_SPACING[1], PIXEL_SPACING[0], SLICE_SPACING, 1.0])
        affine[:3, 3] = [0.0, 0.0, 0.0]

        nii = nib.Nifti1Image(mask_xyz, affine)
        nib.save(nii, mask_dir / f"{name}.nii.gz")

    return mask_dir
