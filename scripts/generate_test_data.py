#!/usr/bin/env python3
"""Generate synthetic CT DICOM series + NIfTI segmentation masks for testing.

Creates:
  data/test/dicom_ct/     — 64 synthetic CT DICOM slices (256×256)
  data/test/segmentations/ — 5 binary NIfTI masks (.nii.gz)

The phantom has simple ellipsoidal "organs" at known positions so we can
verify spatial alignment after DICOM-SEG round-trip.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import nibabel as nib
import numpy as np
import pydicom
from pydicom.dataset import FileDataset
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CT_DIR = PROJECT_ROOT / "data" / "test" / "dicom_ct"
SEG_DIR = PROJECT_ROOT / "data" / "test" / "segmentations"

# Volume geometry
NUM_SLICES = 64
ROWS = 256
COLS = 256
PIXEL_SPACING = [1.0, 1.0]  # mm
SLICE_THICKNESS = 2.5  # mm
SLICE_SPACING = 2.5  # mm

# UIDs (shared across the series)
STUDY_UID = generate_uid()
SERIES_UID = generate_uid()
FRAME_OF_REF_UID = generate_uid()

# Synthetic organ definitions: (name, center_zyx_fraction, radii_zyx_fraction)
ORGANS = [
    ("liver",        (0.45, 0.55, 0.60), (0.15, 0.18, 0.15)),
    ("spleen",       (0.45, 0.55, 0.30), (0.08, 0.08, 0.06)),
    ("kidney_right", (0.50, 0.60, 0.70), (0.07, 0.05, 0.05)),
    ("kidney_left",  (0.50, 0.60, 0.30), (0.07, 0.05, 0.05)),
    ("pancreas",     (0.45, 0.50, 0.50), (0.04, 0.03, 0.12)),
]


def make_ellipsoid_mask(
    shape: tuple[int, int, int],
    center_frac: tuple[float, float, float],
    radii_frac: tuple[float, float, float],
) -> np.ndarray:
    """Create a binary ellipsoid mask in a volume of given shape."""
    zz, yy, xx = np.ogrid[
        0:shape[0], 0:shape[1], 0:shape[2]
    ]
    cz = center_frac[0] * shape[0]
    cy = center_frac[1] * shape[1]
    cx = center_frac[2] * shape[2]
    rz = radii_frac[0] * shape[0]
    ry = radii_frac[1] * shape[1]
    rx = radii_frac[2] * shape[2]

    dist = ((zz - cz) / rz) ** 2 + ((yy - cy) / ry) ** 2 + ((xx - cx) / rx) ** 2
    return (dist <= 1.0).astype(np.uint8)


def create_ct_dicom_series() -> None:
    """Generate a synthetic CT DICOM series."""
    CT_DIR.mkdir(parents=True, exist_ok=True)

    now = datetime.now()
    date_str = now.strftime("%Y%m%d")
    time_str = now.strftime("%H%M%S")

    for i in range(NUM_SLICES):
        sop_uid = generate_uid()
        z_pos = i * SLICE_SPACING  # mm, ascending

        # Create a simple phantom image: body ellipse + some noise
        img = np.random.randint(-1000, -900, (ROWS, COLS), dtype=np.int16)
        # Body ellipse
        yy, xx = np.ogrid[0:ROWS, 0:COLS]
        body = ((yy - ROWS // 2) / (ROWS * 0.4)) ** 2 + ((xx - COLS // 2) / (COLS * 0.35)) ** 2
        img[body <= 1.0] = np.random.randint(-100, 100, img[body <= 1.0].shape, dtype=np.int16)

        # File meta
        file_meta = pydicom.Dataset()
        file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"  # CT
        file_meta.MediaStorageSOPInstanceUID = sop_uid
        file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

        filename = CT_DIR / f"CT_{i:04d}.dcm"
        ds = FileDataset(str(filename), {}, file_meta=file_meta, preamble=b"\x00" * 128)

        # Patient
        ds.PatientName = "Phantom^Test"
        ds.PatientID = "PHANTOM001"
        ds.PatientBirthDate = "19700101"
        ds.PatientSex = "O"

        # Study
        ds.StudyInstanceUID = STUDY_UID
        ds.StudyDate = date_str
        ds.StudyTime = time_str
        ds.StudyDescription = "Synthetic CT for DICOM-SEG trial"
        ds.StudyID = "1"
        ds.AccessionNumber = ""

        # Series
        ds.SeriesInstanceUID = SERIES_UID
        ds.SeriesNumber = 1
        ds.SeriesDescription = "Synthetic CT Axial"
        ds.Modality = "CT"

        # Frame of Reference
        ds.FrameOfReferenceUID = FRAME_OF_REF_UID
        ds.PositionReferenceIndicator = ""

        # Equipment
        ds.Manufacturer = "Synthetic"
        ds.InstitutionName = "Test Lab"

        # Image
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        ds.SOPInstanceUID = sop_uid
        ds.InstanceNumber = i + 1
        ds.ImagePositionPatient = [0.0, 0.0, z_pos]
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]  # standard axial
        ds.PixelSpacing = PIXEL_SPACING
        ds.SliceThickness = SLICE_THICKNESS
        ds.SpacingBetweenSlices = SLICE_SPACING
        ds.SliceLocation = z_pos

        ds.Rows = ROWS
        ds.Columns = COLS
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1  # signed
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.RescaleIntercept = 0
        ds.RescaleSlope = 1
        ds.WindowCenter = 0
        ds.WindowWidth = 400

        ds.PixelData = img.tobytes()

        ds.save_as(filename)

    print(f"Created {NUM_SLICES} CT DICOM slices in {CT_DIR}/")


def create_nifti_masks() -> None:
    """Generate NIfTI segmentation masks matching the CT geometry."""
    SEG_DIR.mkdir(parents=True, exist_ok=True)

    volume_shape = (NUM_SLICES, ROWS, COLS)  # (z, y, x) in image space

    # Build the affine that maps voxel indices to mm coordinates.
    # Our DICOM has:
    #   ImageOrientationPatient = [1,0,0, 0,1,0] (standard axial)
    #   ImagePositionPatient[0] = [0, 0, 0]
    #   PixelSpacing = [1, 1], SliceSpacing = 2.5
    #
    # NIfTI convention is RAS+. Our DICOM ImageOrientationPatient [1,0,0,0,1,0]
    # means: col direction = (1,0,0)=R, row direction = (0,1,0)=A.
    # Slice direction = cross product = (0,0,1)=S.
    # So DICOM voxel (col, row, slice) → (R, A, S) which IS RAS+.
    #
    # nibabel NIfTI shape is (i, j, k) and affine maps (i,j,k) → (x,y,z) in RAS.
    # We store the volume as (x=col, y=row, z=slice) = (R, A, S).
    # So we need to transpose from our (z, y, x) to (x, y, z).

    for name, center, radii in ORGANS:
        mask_zyx = make_ellipsoid_mask(volume_shape, center, radii)

        # Transpose to (x, y, z) = (col, row, slice) for NIfTI (RAS+)
        mask_xyz = np.transpose(mask_zyx, (2, 1, 0))

        # Affine: maps (x_vox, y_vox, z_vox) → (R_mm, A_mm, S_mm)
        affine = np.diag([
            PIXEL_SPACING[1],   # col spacing → R
            PIXEL_SPACING[0],   # row spacing → A
            SLICE_SPACING,      # slice spacing → S
            1.0,
        ])
        # Origin = ImagePositionPatient = (0, 0, 0)
        affine[:3, 3] = [0.0, 0.0, 0.0]

        nii = nib.Nifti1Image(mask_xyz, affine)
        out_path = SEG_DIR / f"{name}.nii.gz"
        nib.save(nii, out_path)

        voxels = int(mask_zyx.sum())
        print(f"  {name}: {voxels:,} voxels → {out_path.name}")

    print(f"Created {len(ORGANS)} NIfTI masks in {SEG_DIR}/")


def main() -> None:
    print("Generating synthetic test data...")
    print()
    create_ct_dicom_series()
    print()
    create_nifti_masks()
    print()
    print("Done! Run the trial script:")
    print("  python scripts/trial_dcmseg.py")


if __name__ == "__main__":
    main()
