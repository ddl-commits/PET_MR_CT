# Convert Module Design Spec

**Date:** 2026-03-22
**Status:** Approved
**Scope:** First implementation pass of the convert module — refactor `scripts/trial_dcmseg.py` into reusable modules with CLI and tests.

## Goal

Replace the proof-of-concept script with a production convert module that creates DICOM-SEG files from NIfTI segmentation masks + source DICOM series, preserving spatial registration and SNOMED CT coded terminology.

## What's In Scope

- `terminology.py` — full implementation with JSON config loading
- `nifti_to_dcmseg.py` — full implementation extracted from POC
- `dcm_to_nifti.py` — typed interface only (raises `NotImplementedError`)
- `nifti_to_rtstruct.py` — typed interface only (raises `NotImplementedError`)
- CLI subcommands for convert
- Populated `segment_metadata_totalseg.json` (~117 TotalSegmentator labels)
- pytest suite: terminology, conversion, round-trip voxel agreement
- Verification: pytest + highdicom round-trip (automated), 3D Slicer (manual spot-check)

## What's Out of Scope

- RT-STRUCT implementation (MIM fallback — deferred)
- DICOM-to-NIfTI conversion (segmentation tools handle this)
- DICOM-SEG splitting for large segment counts (not needed until MIM testing)
- Network module (C-STORE to PACS)
- Segment/register modules

## Module Design

### `convert/terminology.py`

Two public functions:

```python
def load_terminology(tool: str) -> dict[str, SegmentInfo]:
    """Load SNOMED CT terminology mapping from configs/segment_metadata_{tool}.json.

    Args:
        tool: Segmentation tool name ("totalseg", "moose").

    Returns:
        Mapping of label name to SegmentInfo (snomed_code, snomed_meaning,
        category_code, category_meaning).

    Raises:
        FileNotFoundError: If config file does not exist.
        KeyError: If config file has no "segments" key.
    """

def get_segment_description(
    label: str,
    tool: str,
    segment_number: int,
) -> hd.seg.SegmentDescription:
    """Build a highdicom SegmentDescription with SNOMED CT coding.

    Uses CID 7151 (Segmentation Property Category) and CID 7152/7153
    (Segmentation Property Type). Algorithm identification set to the
    tool name.

    Args:
        label: Segment label string (e.g. "liver").
        tool: Segmentation tool name for terminology lookup.
        segment_number: 1-based segment index.

    Returns:
        highdicom SegmentDescription ready for Segmentation constructor.

    Raises:
        KeyError: If label not found in terminology config.
    """
```

Supporting dataclass:

```python
@dataclasses.dataclass(frozen=True)
class SegmentInfo:
    snomed_code: str
    snomed_meaning: str
    category_code: str
    category_meaning: str
    display_label: str  # User-friendly name shown in MIM/3D Slicer (e.g. "Liver")
```

### JSON Config Schema

`configs/segment_metadata_totalseg.json`:

```json
{
  "_comment": "TotalSegmentator label index -> SNOMED CT mapping (CID 7151/7152/7153)",
  "segments": {
    "liver": {
      "snomed_code": "10200004",
      "snomed_meaning": "Liver structure",
      "category_code": "123037004",
      "category_meaning": "Body structure",
      "display_label": "Liver"
    },
    "spleen": { ... },
    ...
  }
}
```

Populate all ~117 TotalSegmentator v2 labels. Leave `segment_metadata_moose.json` as an empty stub for now.

### `convert/nifti_to_dcmseg.py`

One public function:

```python
def create_dcmseg(
    source_dcm_dir: Path,
    nifti_masks: dict[str, Path],
    tool: str,
    output_path: Path,
) -> Path:
    """Create a DICOM-SEG file from NIfTI segmentation masks.

    Loads the source DICOM series for geometry and Frame of Reference,
    combines NIfTI masks into a multi-label pixel array, and writes
    a DICOM-SEG file via highdicom.

    Args:
        source_dcm_dir: Directory containing the source DICOM series.
        nifti_masks: Mapping of label name to NIfTI file path.
        tool: Segmentation tool name for terminology lookup.
        output_path: Where to write the DICOM-SEG file.

    Returns:
        Path to the written DICOM-SEG file.
    """
```

Internal responsibilities (extracted from `scripts/trial_dcmseg.py`):

1. **Load source DICOMs** — try `.dcm` files first; if none found, iterate all non-hidden files and attempt `pydicom.dcmread()`, catching `InvalidDicomError` (logged, skipped). Sort by `ImagePositionPatient[2]`.
2. **Load NIfTI masks** — for each mask:
   - Skip masks whose path does not exist, logging a warning
   - `nibabel.as_closest_canonical()` for consistent orientation
   - Flip axes 0 and 1 for RAS+ to LPS+ conversion
   - Transpose from (L, P, S) to (slice, row, col) for DICOM pixel array
   - Handle z-axis inversion when affine has negative z-spacing
   - Verify dimensions match DICOM grid (rows, cols, slices); try rows/cols swap if transposed; raise `ValueError` if no match
3. **Build segment descriptions** — via `terminology.get_segment_description()`
4. **Create highdicom Segmentation** — binary type, `omit_empty_frames=True`, `series_number=100` (default), manufacturer/software metadata sourced from source DICOM series
5. **Write to disk** — `seg.save_as(output_path)`

**Error handling:**
- Empty `nifti_masks` dict → raise `ValueError("No masks provided")`
- All masks skipped (none exist on disk) → raise `ValueError("No valid masks found")`

### `convert/dcm_to_nifti.py` (stub)

```python
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
```

### `convert/nifti_to_rtstruct.py` (stub)

```python
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
```

### `convert/__init__.py`

Re-export public API:

```python
from pet_mr_ct.convert.nifti_to_dcmseg import create_dcmseg
from pet_mr_ct.convert.nifti_to_rtstruct import create_rtstruct
from pet_mr_ct.convert.dcm_to_nifti import convert_dcm_to_nifti
from pet_mr_ct.convert.terminology import load_terminology, get_segment_description
```

### CLI (`cli.py`)

```
pet-mr-ct convert dcmseg --source-dir PATH --mask-dir PATH --tool TEXT --output PATH
pet-mr-ct convert rtstruct --source-dir PATH --mask-dir PATH --output PATH  (-> "not yet implemented")
```

- `convert` is a Click group under the main `pet-mr-ct` group
- `dcmseg` command: discovers all `.nii.gz` and `.nii` files in `--mask-dir`, builds label→path dict from filenames (stem without extension), calls `create_dcmseg()`
- `rtstruct` command: prints message that RT-STRUCT is not yet implemented, exits cleanly

## Test Design

### `tests/conftest.py`

Shared fixtures using logic from `scripts/generate_test_data.py`:

- `tmp_dicom_series(tmp_path)` — generates 8-slice synthetic CT DICOM series (small, fast)
- `tmp_nifti_masks(tmp_path)` — generates 3 binary NIfTI masks (liver, spleen, kidney_right) with known geometry matching the DICOM series
- Keep fixtures small: 64x64 pixels, 8 slices, 3 segments

### `tests/test_terminology.py`

- `test_load_totalseg_terminology` — loads config, verifies liver SNOMED code is "10200004"
- `test_unknown_label_raises` — requesting unknown label raises `KeyError`
- `test_get_segment_description_returns_valid_type` — returns `hd.seg.SegmentDescription`

### `tests/test_convert.py`

- `test_create_dcmseg_produces_file` — output file exists and is non-empty
- `test_dcmseg_has_correct_segment_count` — read back, verify 3 segments
- `test_dcmseg_segment_labels` — read back, verify labels match input names

### `tests/test_roundtrip.py`

- `test_voxel_agreement` — create DICOM-SEG, read back with highdicom, extract pixel array per segment, compare against input NIfTI masks, require >99% voxel agreement

## Verification Strategy

- **Automated:** pytest suite covers terminology loading, DICOM-SEG creation, metadata correctness, and voxel-level round-trip agreement
- **Manual:** Load output DICOM-SEG in 3D Slicer to visually confirm segments render at correct anatomical locations with correct labels
- **No dciodvfy** in this pass (external C binary; add later if needed)

## Dependencies

All already declared in `pyproject.toml`:
- pydicom >= 2.4.0
- highdicom >= 0.22.0
- nibabel >= 5.2
- click >= 8.1
- numpy
- pytest (dev dependency)

## File Changes Summary

| File | Action |
|------|--------|
| `src/pet_mr_ct/convert/__init__.py` | Rewrite — re-export public API |
| `src/pet_mr_ct/convert/terminology.py` | Rewrite — full implementation |
| `src/pet_mr_ct/convert/nifti_to_dcmseg.py` | Rewrite — full implementation |
| `src/pet_mr_ct/convert/dcm_to_nifti.py` | Rewrite — typed stub |
| `src/pet_mr_ct/convert/nifti_to_rtstruct.py` | Rewrite — typed stub |
| `src/pet_mr_ct/cli.py` | Rewrite — add convert group + subcommands |
| `configs/segment_metadata_totalseg.json` | Populate ~117 TotalSegmentator labels |
| `tests/conftest.py` | Rewrite — shared fixtures |
| `tests/test_terminology.py` | Rewrite — terminology tests |
| `tests/test_convert.py` | Rewrite — conversion tests |
| `tests/test_roundtrip.py` | Rewrite — voxel agreement tests |
