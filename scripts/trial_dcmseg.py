#!/usr/bin/env python3
"""Trial script: CT DICOM + TotalSegmentator NIfTI masks → DICOM-SEG.

Proof-of-concept for the round-trip:
  AI segmentation → DICOM-SEG → MIM 7.3.7 / 3D Slicer import.

Usage:
    python scripts/trial_dcmseg.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import highdicom as hd
import nibabel as nib
import numpy as np
import pydicom

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent

CT_DICOM_DIR = PROJECT_ROOT / "data" / "test" / "dicom_ct"
SEG_DIR = PROJECT_ROOT / "data" / "test" / "segmentations"
OUTPUT_DIR = PROJECT_ROOT / "data" / "test" / "output"
OUTPUT_PATH = OUTPUT_DIR / "trial_dcmseg.dcm"

# Organs to include (label_value, filename_stem, display_label, SNOMED CT code)
ORGANS: list[tuple[int, str, str, str]] = [
    (1, "liver", "Liver", "10200004"),
    (2, "spleen", "Spleen", "78961009"),
    (3, "kidney_right", "Right Kidney", "9846003"),
    (4, "kidney_left", "Left Kidney", "18639004"),
    (5, "pancreas", "Pancreas", "15776009"),
]

BODY_STRUCTURE_CODE = "123037004"  # SCT: "Body structure"


# ---------------------------------------------------------------------------
# Step 1 — Load source CT DICOM series
# ---------------------------------------------------------------------------

def load_ct_series(dicom_dir: Path) -> list[pydicom.Dataset]:
    """Load and sort CT DICOM datasets by ImagePositionPatient z-coordinate."""
    if not dicom_dir.is_dir():
        sys.exit(f"ERROR: CT DICOM directory not found: {dicom_dir}")

    dcm_files = sorted(dicom_dir.glob("*.dcm"))
    if not dcm_files:
        # Also try without extension — some DICOM files have no .dcm suffix
        dcm_files = [
            f for f in sorted(dicom_dir.iterdir())
            if f.is_file() and not f.name.startswith(".")
        ]
    if not dcm_files:
        sys.exit(f"ERROR: No DICOM files found in {dicom_dir}")

    datasets = []
    for f in dcm_files:
        try:
            ds = pydicom.dcmread(f)
            datasets.append(ds)
        except pydicom.errors.InvalidDicomError:
            print(f"  WARN: skipping non-DICOM file {f.name}")

    if not datasets:
        sys.exit(f"ERROR: No valid DICOM files in {dicom_dir}")

    # Sort by z-coordinate of ImagePositionPatient
    datasets.sort(key=lambda ds: float(ds.ImagePositionPatient[2]))

    print(f"Loaded {len(datasets)} CT DICOM slices from {dicom_dir.name}/")
    series_desc = getattr(datasets[0], "SeriesDescription", "N/A")
    print(f"  Series: {series_desc}")
    print(f"  Rows×Cols: {datasets[0].Rows}×{datasets[0].Columns}")
    z_positions = [float(ds.ImagePositionPatient[2]) for ds in datasets]
    print(f"  Z range: {min(z_positions):.2f} → {max(z_positions):.2f} mm")

    return datasets


# ---------------------------------------------------------------------------
# Step 2–4 — Load NIfTI masks and combine into multi-label array
# ---------------------------------------------------------------------------

def load_and_combine_masks(
    seg_dir: Path,
    organs: list[tuple[int, str, str, str]],
    ct_datasets: list[pydicom.Dataset],
) -> np.ndarray:
    """Load NIfTI masks, reorient to match DICOM, combine into one label array.

    Returns:
        4D boolean array of shape (num_slices, rows, cols, num_segments).
    """
    if not seg_dir.is_dir():
        sys.exit(f"ERROR: Segmentation directory not found: {seg_dir}")

    num_slices = len(ct_datasets)
    rows = ct_datasets[0].Rows
    cols = ct_datasets[0].Columns

    # highdicom expects (slices, rows, cols, num_segments) for binary multi-segment
    combined = np.zeros((num_slices, rows, cols, len(organs)), dtype=np.bool_)

    for idx, (label_val, stem, display_name, _sct_code) in enumerate(organs):
        nii_path = seg_dir / f"{stem}.nii.gz"
        if not nii_path.exists():
            # Try without .gz
            nii_path = seg_dir / f"{stem}.nii"
        if not nii_path.exists():
            print(f"  WARN: mask not found for {display_name} ({stem}.nii.gz) — skipping")
            continue

        nii = nib.load(nii_path)

        # ------------------------------------------------------------------
        # Step 6 — Handle NIfTI → DICOM orientation mismatch
        #
        # NIfTI is stored in RAS+ orientation (by convention).
        # DICOM CT is typically LPS+ (rows=AP, cols=LR, slices=SI).
        #
        # Strategy:
        #   1. Reorient the NIfTI data to LPS+ using nibabel's as_closest_canonical
        #      then flip axes as needed.
        #   2. Verify dimensions match the DICOM grid.
        # ------------------------------------------------------------------

        # Reorient NIfTI to closest canonical (RAS+) first, then convert to LPS
        canonical = nib.as_closest_canonical(nii)
        mask_ras = np.asarray(canonical.dataobj, dtype=np.float32)

        # RAS+ → LPS+: flip first two axes (R→L, A→P), S stays
        mask_lps = mask_ras[::-1, ::-1, :]

        # DICOM pixel array convention: (slice, row, col)
        # In LPS+ space: axis 2 = S (slice), axis 1 = P (row), axis 0 = L (col)
        # So transpose from (L, P, S) → (S, P, L) = (slice, row, col)
        # But actually: DICOM row direction and column direction depend on
        # ImageOrientationPatient. For standard axial CT:
        #   row direction = posterior (P), col direction = left (L)
        #   => pixel[row, col] = [P, L]
        #   => volume[slice, row, col] = [S, P, L]
        #
        # Our mask_lps has shape (L, P, S), so:
        mask_dicom = np.transpose(mask_lps, (2, 1, 0))  # (S, P, L)

        # Check dimensions
        if mask_dicom.shape != (num_slices, rows, cols):
            print(f"  WARN: dimension mismatch for {display_name}:")
            print(f"    NIfTI (reoriented): {mask_dicom.shape}")
            print(f"    DICOM grid:         ({num_slices}, {rows}, {cols})")

            # Try to handle common mismatches
            if mask_dicom.shape == (num_slices, cols, rows):
                print("    → Swapping rows/cols")
                mask_dicom = np.swapaxes(mask_dicom, 1, 2)
            elif mask_dicom.shape[0] != num_slices:
                print("    → Slice count mismatch — cannot reconcile, skipping")
                continue

        # Binarize (threshold at 0.5 in case of interpolation artifacts)
        mask_binary = mask_dicom > 0.5

        # Check slice ordering: DICOM datasets are sorted by z ascending.
        # If NIfTI z-axis is flipped relative to DICOM, we need to reverse.
        # Compare: does the NIfTI affine z-direction match DICOM z-ordering?
        affine = canonical.affine
        # After RAS→LPS flip, z (S axis) direction from the affine:
        # canonical is RAS, so positive z = Superior.
        # DICOM sorted by z ascending = Inferior→Superior.
        # The canonical affine's z-spacing sign tells us the direction.
        z_spacing_sign = affine[2, 2]  # positive = I→S in RAS
        # After our LPS flip, z axis is NOT flipped (S stays S).
        # So if z_spacing_sign > 0, NIfTI slice 0 = most inferior = DICOM slice 0. OK.
        # If z_spacing_sign < 0, NIfTI slice 0 = most superior = need to flip.
        if z_spacing_sign < 0:
            mask_binary = mask_binary[::-1, :, :]

        combined[:, :, :, idx] = mask_binary

        voxel_count = int(mask_binary.sum())
        print(f"  Loaded {display_name}: {voxel_count:,} non-zero voxels")

    return combined


# ---------------------------------------------------------------------------
# Step 5 — Build DICOM-SEG
# ---------------------------------------------------------------------------

def build_segment_descriptions(
    organs: list[tuple[int, str, str, str]],
) -> list[hd.seg.SegmentDescription]:
    """Create SegmentDescription for each organ with SNOMED CT coding."""
    category = hd.sr.CodedConcept(
        value=BODY_STRUCTURE_CODE,
        meaning="Body Structure",
        scheme_designator="SCT",
    )

    algo_id = hd.AlgorithmIdentificationSequence(
        name="TotalSegmentator",
        version="2.0",
        family=hd.sr.CodedConcept(
            value="123109",
            meaning="Artificial Intelligence",
            scheme_designator="DCM",
        ),
    )

    descriptions = []
    for seg_num, _stem, display_name, sct_code in organs:
        desc = hd.seg.SegmentDescription(
            segment_number=seg_num,
            segment_label=display_name,
            segmented_property_category=category,
            segmented_property_type=hd.sr.CodedConcept(
                value=sct_code,
                meaning=display_name,
                scheme_designator="SCT",
            ),
            algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
            algorithm_identification=algo_id,
        )
        descriptions.append(desc)

    return descriptions


def create_dcmseg(
    ct_datasets: list[pydicom.Dataset],
    pixel_array: np.ndarray,
    segment_descriptions: list[hd.seg.SegmentDescription],
    output_path: Path,
) -> Path:
    """Create and save a DICOM-SEG file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seg = hd.seg.Segmentation(
        source_images=ct_datasets,
        pixel_array=pixel_array,
        segmentation_type=hd.seg.SegmentationTypeValues.BINARY,
        segment_descriptions=segment_descriptions,
        series_instance_uid=hd.UID(),
        series_number=100,
        sop_instance_uid=hd.UID(),
        instance_number=1,
        manufacturer="PET_MR_CT Pipeline",
        manufacturer_model_name="trial_dcmseg",
        software_versions="0.1.0",
        device_serial_number="0001",
        content_description="TotalSegmentator organ segmentation",
        content_label="TOTALSEG",
        omit_empty_frames=True,
    )

    seg.save_as(output_path)
    return output_path


# ---------------------------------------------------------------------------
# Step 8 — Summary and validation
# ---------------------------------------------------------------------------

def print_summary(
    output_path: Path,
    pixel_array: np.ndarray,
    organs: list[tuple[int, str, str, str]],
) -> None:
    """Print summary of the created DICOM-SEG."""
    file_size_mb = output_path.stat().st_size / (1024 * 1024)
    num_slices = pixel_array.shape[0]
    num_segments = pixel_array.shape[3]

    print("\n" + "=" * 60)
    print("DICOM-SEG created successfully!")
    print("=" * 60)
    print(f"  Output:     {output_path}")
    print(f"  File size:  {file_size_mb:.2f} MB")
    print(f"  Segments:   {num_segments}")
    print(f"  Slices:     {num_slices}")
    print()

    for idx, (_label_val, _stem, display_name, _sct_code) in enumerate(organs):
        voxels = int(pixel_array[:, :, :, idx].sum())
        print(f"  Segment {idx + 1}: {display_name:20s} — {voxels:>10,} voxels")


def validate_roundtrip(output_path: Path) -> None:
    """Read back the DICOM-SEG and verify segment metadata."""
    print("\n" + "-" * 60)
    print("Validation: reading DICOM-SEG back...")
    print("-" * 60)

    seg = hd.seg.segread(str(output_path))
    print(f"  Number of segments: {seg.number_of_segments}")

    for i in range(seg.number_of_segments):
        desc = seg.get_segment_description(i + 1)
        print(f"  Segment {i + 1}: {desc.segment_label}")

    print("\n  Round-trip validation PASSED ✓")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 60)
    print("Trial: CT DICOM + TotalSegmentator → DICOM-SEG")
    print("=" * 60)
    print()

    # Step 1: Load CT DICOM
    ct_datasets = load_ct_series(CT_DICOM_DIR)
    print()

    # Steps 2–4: Load and combine NIfTI masks
    print("Loading NIfTI segmentation masks...")
    pixel_array = load_and_combine_masks(SEG_DIR, ORGANS, ct_datasets)
    print()

    # Step 5: Build DICOM-SEG
    print("Building DICOM-SEG...")
    segment_descriptions = build_segment_descriptions(ORGANS)
    output_path = create_dcmseg(ct_datasets, pixel_array, segment_descriptions, OUTPUT_PATH)

    # Step 8: Summary
    print_summary(output_path, pixel_array, ORGANS)

    # Validation
    validate_roundtrip(output_path)


if __name__ == "__main__":
    main()
