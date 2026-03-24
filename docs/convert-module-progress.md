# Convert Module — Implementation Progress

**Date:** 2026-03-23
**Status:** Complete (first pass)
**Branch:** main

---

## What Was Built

The convert module refactors the proof-of-concept script (`scripts/trial_dcmseg.py`) into a production Python package that creates DICOM-SEG files from AI segmentation NIfTI masks + source DICOM series.

### Pipeline Flow

```
NIfTI masks (.nii.gz)  ──┐
                          ├──> create_dcmseg() ──> DICOM-SEG (.dcm)
Source DICOM series    ──┘        │
                                  ├── RAS+ → LPS+ coordinate conversion
                                  ├── SNOMED CT coded terminology
                                  └── highdicom Segmentation object
```

### Files Implemented

| File | Description | Status |
|------|-------------|--------|
| `configs/segment_metadata_totalseg.json` | 117 TotalSegmentator labels → SNOMED CT codes | Done |
| `src/pet_mr_ct/convert/terminology.py` | `load_terminology()`, `get_segment_description()` | Done |
| `src/pet_mr_ct/convert/nifti_to_dcmseg.py` | `create_dcmseg()` — core conversion | Done |
| `src/pet_mr_ct/convert/dcm_to_nifti.py` | Typed stub (`NotImplementedError`) | Stub |
| `src/pet_mr_ct/convert/nifti_to_rtstruct.py` | Typed stub (`NotImplementedError`) | Stub |
| `src/pet_mr_ct/convert/__init__.py` | Public API re-exports | Done |
| `src/pet_mr_ct/cli.py` | `pet-mr-ct convert dcmseg` command | Done |

### Test Coverage

| Test File | Tests | Description |
|-----------|-------|-------------|
| `tests/test_terminology.py` | 6 | Config loading, SNOMED codes, error handling |
| `tests/test_convert.py` | 5 | DICOM-SEG creation, segment count, labels, error cases |
| `tests/test_roundtrip.py` | 1 | Voxel-level round-trip agreement (>99%) |
| **Total** | **12** | All passing (0.32s) |

---

## Verification Results

### End-to-End CLI Test

```
$ pet-mr-ct convert dcmseg \
    --source-dir data/test/dicom_ct \
    --mask-dir data/test/segmentations \
    --tool totalseg \
    --output data/test/output/cli_test.dcm

Found 117 mask(s)
DICOM-SEG written to data/test/output/cli_test.dcm
```

| Metric | Value |
|--------|-------|
| Output file size | 154 MB |
| Total segments | 117 |
| Segments with voxels | 115 |
| Empty segments | 2 |
| Total segmented voxels | 3,025,808 |
| CT slices | 345 |
| CT resolution | 512 x 512 |

### Top 5 Segments by Voxel Count

| Segment | Voxels |
|---------|--------|
| Liver | 350,427 |
| Brain | 180,984 |
| Small Bowel | 163,948 |
| Right Gluteus Maximus | 133,972 |
| Left Gluteus Maximus | 126,814 |

### Visual Verification

Overlay figure (`data/test/output/verify_overlays.png`) shows segments correctly registered to CT anatomy across 4 representative axial slices (head, chest, abdomen, pelvis). Voxel count distribution (`data/test/output/verify_voxel_counts.png`) shows expected anatomical size ordering.

---

## Key Technical Decisions

1. **Coordinate handling:** NIfTI RAS+ → DICOM LPS+ via nibabel canonical reorientation + axis flips + transpose. Z-axis inversion detected from affine sign.

2. **DICOM loading:** Tries `.dcm` extension first, falls back to all non-hidden files with `InvalidDicomError` catching. Handles real-world DICOM directories.

3. **Terminology caching:** `load_terminology()` uses `@functools.lru_cache` to avoid re-reading JSON on every segment.

4. **omit_empty_frames:** DICOM-SEG created with `omit_empty_frames=True` for smaller file size. Round-trip test handles sparse frames via `assert_missing_frames_are_empty=True`.

5. **Tool-agnostic CLI:** `--tool` accepts any string matching `configs/segment_metadata_{tool}.json`, not a hardcoded list.

---

## What's Next (Not Yet Implemented)

| Module | Purpose | Priority |
|--------|---------|----------|
| RT-STRUCT fallback | MIM may struggle with large DICOM-SEG; RT-STRUCT is native | Medium |
| DICOM-to-NIfTI | Standalone conversion (segmentation tools do this already) | Low |
| Segment module | TotalSegmentator/MOOSE CLI wrappers | Medium |
| Network module | pynetdicom C-STORE to MIM PACS | High (needed for clinical use) |
| Register module | ANTsPy cross-timepoint registration | Low |

### Open Questions for Discussion

1. **MIM compatibility:** The 154 MB / 117-segment DICOM-SEG needs testing in MIM 7.3.7. If it fails, we have two options: split into per-region files or fall back to RT-STRUCT.

2. **SNOMED CT code validation:** 5 codes are verified from the POC. The remaining ~112 were populated from training data knowledge and should be cross-checked against the SNOMED CT browser before production use.

3. **PET DICOM support:** The test data includes 2,070 real PET DICOM files. The convert module was tested with CT — PET may have different orientation/geometry conventions worth testing.

4. **Next priority:** Network module (C-STORE to MIM) or segment module (TotalSegmentator wrapper)?
