# Convert Module Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor the POC script (`scripts/trial_dcmseg.py`) into a production convert module that creates DICOM-SEG files from NIfTI masks + source DICOM series, with CLI and tests.

**Architecture:** The convert module has 4 files: `terminology.py` (SNOMED CT lookup from JSON configs), `nifti_to_dcmseg.py` (core conversion with RAS+→LPS+ coordinate handling via highdicom), and two stubs (`dcm_to_nifti.py`, `nifti_to_rtstruct.py`). CLI exposes `pet-mr-ct convert dcmseg`. Tests use synthetic DICOM+NIfTI fixtures for fast, deterministic validation.

**Tech Stack:** Python 3.11+, pydicom, highdicom, nibabel, numpy, click, pytest

**Spec:** `docs/superpowers/specs/2026-03-22-convert-module-design.md`

**Reference POC:** `scripts/trial_dcmseg.py` (working end-to-end proof-of-concept)

**Reference test data generator:** `scripts/generate_test_data.py` (synthetic DICOM + NIfTI creation logic)

---

## File Structure

| File | Responsibility |
|------|---------------|
| `configs/segment_metadata_totalseg.json` | SNOMED CT mappings for all ~117 TotalSegmentator v2 labels |
| `src/pet_mr_ct/convert/terminology.py` | `SegmentInfo` dataclass, `load_terminology()`, `get_segment_description()` |
| `src/pet_mr_ct/convert/nifti_to_dcmseg.py` | `create_dcmseg()` — DICOM loading, NIfTI→LPS+ conversion, highdicom Segmentation creation |
| `src/pet_mr_ct/convert/dcm_to_nifti.py` | Typed stub — `convert_dcm_to_nifti()` raises `NotImplementedError` |
| `src/pet_mr_ct/convert/nifti_to_rtstruct.py` | Typed stub — `create_rtstruct()` raises `NotImplementedError` |
| `src/pet_mr_ct/convert/__init__.py` | Re-exports public API |
| `src/pet_mr_ct/cli.py` | Click group `convert` with subcommands `dcmseg` and `rtstruct` |
| `tests/conftest.py` | Shared fixtures: `tmp_dicom_series`, `tmp_nifti_masks` (64×64, 8 slices, 3 segments) |
| `tests/test_terminology.py` | Tests for terminology loading, unknown labels, segment description building |
| `tests/test_convert.py` | Tests for DICOM-SEG creation, segment count, labels |
| `tests/test_roundtrip.py` | Voxel-level round-trip agreement test |

**Important ordering note:** Tasks 2-5 depend on `src/pet_mr_ct/convert/__init__.py` remaining minimal (just a docstring). Do NOT add re-exports to `__init__.py` until Task 6, after all modules exist. The existing stub `__init__.py` is safe.

---

## Task 1: Populate TotalSegmentator SNOMED CT Config

**Files:**
- Modify: `configs/segment_metadata_totalseg.json`

This is a data-only task — no code. The JSON must contain all ~117 TotalSegmentator v2 label names mapped to SNOMED CT codes, meanings, categories, and display labels.

- [ ] **Step 1: Write the populated JSON config**

Replace the empty stub with all TotalSegmentator v2 labels. Each entry follows this schema:

```json
{
  "label_name": {
    "snomed_code": "SNOMED_CT_CODE",
    "snomed_meaning": "SNOMED formal term",
    "category_code": "123037004",
    "category_meaning": "Body structure",
    "display_label": "Human-readable label"
  }
}
```

The full TotalSegmentator v2 label list (use these exact filename stems as keys):
`spleen`, `kidney_right`, `kidney_left`, `gallbladder`, `liver`, `stomach`, `pancreas`, `adrenal_gland_right`, `adrenal_gland_left`, `lung_upper_lobe_left`, `lung_lower_lobe_left`, `lung_upper_lobe_right`, `lung_middle_lobe_right`, `lung_lower_lobe_right`, `esophagus`, `trachea`, `thyroid_gland`, `small_bowel`, `duodenum`, `colon`, `urinary_bladder`, `prostate`, `kidney_cyst_left`, `kidney_cyst_right`, `sacrum`, `vertebrae_S1`, `vertebrae_L5`, `vertebrae_L4`, `vertebrae_L3`, `vertebrae_L2`, `vertebrae_L1`, `vertebrae_T12`, `vertebrae_T11`, `vertebrae_T10`, `vertebrae_T9`, `vertebrae_T8`, `vertebrae_T7`, `vertebrae_T6`, `vertebrae_T5`, `vertebrae_T4`, `vertebrae_T3`, `vertebrae_T2`, `vertebrae_T1`, `vertebrae_C7`, `vertebrae_C6`, `vertebrae_C5`, `vertebrae_C4`, `vertebrae_C3`, `vertebrae_C2`, `vertebrae_C1`, `heart`, `aorta`, `pulmonary_vein`, `brachiocephalic_trunk`, `subclavian_artery_right`, `subclavian_artery_left`, `common_carotid_artery_right`, `common_carotid_artery_left`, `brachiocephalic_vein_left`, `brachiocephalic_vein_right`, `atrial_appendage_left`, `superior_vena_cava`, `inferior_vena_cava`, `portal_vein_and_splenic_vein`, `iliac_artery_left`, `iliac_artery_right`, `iliac_vena_left`, `iliac_vena_right`, `humerus_left`, `humerus_right`, `scapula_left`, `scapula_right`, `clavicula_left`, `clavicula_right`, `femur_left`, `femur_right`, `hip_left`, `hip_right`, `spinal_cord`, `gluteus_maximus_left`, `gluteus_maximus_right`, `gluteus_medius_left`, `gluteus_medius_right`, `gluteus_minimus_left`, `gluteus_minimus_right`, `autochthon_left`, `autochthon_right`, `iliopsoas_left`, `iliopsoas_right`, `brain`, `skull`, `rib_left_1`, `rib_left_2`, `rib_left_3`, `rib_left_4`, `rib_left_5`, `rib_left_6`, `rib_left_7`, `rib_left_8`, `rib_left_9`, `rib_left_10`, `rib_left_11`, `rib_left_12`, `rib_right_1`, `rib_right_2`, `rib_right_3`, `rib_right_4`, `rib_right_5`, `rib_right_6`, `rib_right_7`, `rib_right_8`, `rib_right_9`, `rib_right_10`, `rib_right_11`, `rib_right_12`, `sternum`, `costal_cartilages`

For each label, look up the correct SNOMED CT code. Use category `"123037004"` / `"Body structure"` for anatomical structures. For specific categories use:
- Bones/skeletal: category `"123037004"` / `"Body structure"`
- Organs: category `"123037004"` / `"Body structure"`
- Vessels: category `"123037004"` / `"Body structure"`
- Muscles: category `"123037004"` / `"Body structure"`

Reference: The POC uses these known-good codes: liver=`10200004`, spleen=`78961009`, kidney_right=`9846003`, kidney_left=`18639004`, pancreas=`15776009`.

- [ ] **Step 2: Validate JSON is parseable**

Run: `python -c "import json; json.load(open('configs/segment_metadata_totalseg.json')); print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add configs/segment_metadata_totalseg.json
git commit -m "feat: populate TotalSegmentator SNOMED CT terminology config (~117 labels)"
```

---

## Task 2: Implement `terminology.py`

**Files:**
- Create: `src/pet_mr_ct/convert/terminology.py`
- Test: `tests/test_terminology.py`

**Reference:** POC `build_segment_descriptions()` at `scripts/trial_dcmseg.py:197-233`

- [ ] **Step 1: Write the failing tests**

Write `tests/test_terminology.py`:

```python
"""Tests for SNOMED CT terminology mappings."""

from __future__ import annotations

import pytest

import highdicom as hd

from pet_mr_ct.convert.terminology import (
    SegmentInfo,
    get_segment_description,
    load_terminology,
)


def test_load_totalseg_terminology():
    """Loading totalseg config returns dict with known liver entry."""
    terms = load_terminology("totalseg")
    assert "liver" in terms
    assert isinstance(terms["liver"], SegmentInfo)
    assert terms["liver"].snomed_code == "10200004"
    assert terms["liver"].display_label == "Liver"


def test_load_totalseg_has_all_labels():
    """Config should have at least 100 labels for TotalSegmentator."""
    terms = load_terminology("totalseg")
    assert len(terms) >= 100


def test_load_unknown_tool_raises():
    """Loading a nonexistent tool config raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        load_terminology("nonexistent_tool")


def test_unknown_label_raises():
    """Requesting an unknown label raises KeyError."""
    with pytest.raises(KeyError, match="not_a_real_organ"):
        get_segment_description("not_a_real_organ", "totalseg", 1)


def test_get_segment_description_returns_valid_type():
    """get_segment_description returns a highdicom SegmentDescription."""
    desc = get_segment_description("liver", "totalseg", 1)
    assert isinstance(desc, hd.seg.SegmentDescription)
    assert desc.segment_number == 1
    assert desc.segment_label == "Liver"


def test_get_segment_description_uses_snomed_codes():
    """Segment description has correct SNOMED CT property type."""
    desc = get_segment_description("spleen", "totalseg", 2)
    assert desc.segmented_property_type.value == "78961009"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_terminology.py -v`
Expected: FAIL (ImportError — module not implemented yet)

- [ ] **Step 3: Write the implementation**

Write `src/pet_mr_ct/convert/terminology.py`:

```python
"""SNOMED CT coded terminology for organ segment labels."""

from __future__ import annotations

import dataclasses
import json
import functools
import logging
from importlib import resources
from pathlib import Path

import highdicom as hd

logger = logging.getLogger(__name__)


def _get_configs_dir() -> Path:
    """Resolve the configs directory, handling both dev and installed layouts."""
    # Development layout: src/pet_mr_ct/convert/terminology.py → project_root/configs/
    dev_path = Path(__file__).resolve().parent.parent.parent.parent / "configs"
    if dev_path.is_dir():
        return dev_path
    # Fallback: package data via importlib.resources
    try:
        return Path(str(resources.files("pet_mr_ct") / ".." / ".." / "configs"))
    except (TypeError, FileNotFoundError):
        raise FileNotFoundError(
            "Cannot locate configs directory. Ensure the package is installed "
            "with configs/ accessible or set PET_MR_CT_CONFIGS_DIR env var."
        )


@dataclasses.dataclass(frozen=True)
class SegmentInfo:
    """SNOMED CT terminology for a single segment label."""

    snomed_code: str
    snomed_meaning: str
    category_code: str
    category_meaning: str
    display_label: str


@functools.lru_cache(maxsize=8)
def load_terminology(tool: str) -> dict[str, SegmentInfo]:
    """Load SNOMED CT terminology mapping from configs/segment_metadata_{tool}.json.

    Results are cached — repeated calls with the same tool name are free.

    Args:
        tool: Segmentation tool name ("totalseg", "moose").

    Returns:
        Mapping of label name to SegmentInfo.

    Raises:
        FileNotFoundError: If config file does not exist.
        KeyError: If config file has no "segments" key.
    """
    config_path = _get_configs_dir() / f"segment_metadata_{tool}.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Terminology config not found: {config_path}")

    with open(config_path) as f:
        data = json.load(f)

    if "segments" not in data:
        raise KeyError(f"Config file missing 'segments' key: {config_path}")

    result: dict[str, SegmentInfo] = {}
    for label, info in data["segments"].items():
        result[label] = SegmentInfo(
            snomed_code=info["snomed_code"],
            snomed_meaning=info["snomed_meaning"],
            category_code=info["category_code"],
            category_meaning=info["category_meaning"],
            display_label=info["display_label"],
        )
    return result


def get_segment_description(
    label: str,
    tool: str,
    segment_number: int,
) -> hd.seg.SegmentDescription:
    """Build a highdicom SegmentDescription with SNOMED CT coding.

    Args:
        label: Segment label string (e.g. "liver").
        tool: Segmentation tool name for terminology lookup.
        segment_number: 1-based segment index.

    Returns:
        highdicom SegmentDescription ready for Segmentation constructor.

    Raises:
        KeyError: If label not found in terminology config.
    """
    terms = load_terminology(tool)
    if label not in terms:
        raise KeyError(
            f"Label '{label}' not found in {tool} terminology config. "
            f"Available labels: {', '.join(sorted(terms.keys())[:10])}..."
        )

    info = terms[label]

    category = hd.sr.CodedConcept(
        value=info.category_code,
        meaning=info.category_meaning,
        scheme_designator="SCT",
    )

    property_type = hd.sr.CodedConcept(
        value=info.snomed_code,
        meaning=info.display_label,
        scheme_designator="SCT",
    )

    algo_id = hd.AlgorithmIdentificationSequence(
        name=tool,
        version="1.0",
        family=hd.sr.CodedConcept(
            value="123109",
            meaning="Artificial Intelligence",
            scheme_designator="DCM",
        ),
    )

    return hd.seg.SegmentDescription(
        segment_number=segment_number,
        segment_label=info.display_label,
        segmented_property_category=category,
        segmented_property_type=property_type,
        algorithm_type=hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC,
        algorithm_identification=algo_id,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_terminology.py -v`
Expected: All 6 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/pet_mr_ct/convert/terminology.py tests/test_terminology.py
git commit -m "feat: implement terminology module with SNOMED CT lookup and tests"
```

---

## Task 3: Create Test Fixtures

**Files:**
- Create: `tests/conftest.py`

**Reference:** `scripts/generate_test_data.py` for synthetic data creation logic. Use the same geometry approach but smaller (64×64, 8 slices) for speed.

- [ ] **Step 1: Write conftest.py with shared fixtures**

```python
"""Shared pytest fixtures for PET_MR_CT tests."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pydicom
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

import pytest

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
```

- [ ] **Step 2: Verify fixtures work**

Run: `python -c "import tests.conftest; print('Fixtures importable')"`
Expected: `Fixtures importable` (or may fail if not on sys.path — that's OK, pytest handles this)

Run: `pytest --collect-only tests/`
Expected: Shows collected tests from test_terminology.py, fixtures listed

- [ ] **Step 3: Commit**

```bash
git add tests/conftest.py
git commit -m "feat: add shared pytest fixtures for synthetic DICOM and NIfTI test data"
```

---

## Task 4: Implement `nifti_to_dcmseg.py`

**Files:**
- Create: `src/pet_mr_ct/convert/nifti_to_dcmseg.py`
- Test: `tests/test_convert.py`

**Reference:** POC functions `load_ct_series()` at `scripts/trial_dcmseg.py:48-84`, `load_and_combine_masks()` at `scripts/trial_dcmseg.py:91-190`, and `create_dcmseg()` at `scripts/trial_dcmseg.py:236-264`.

- [ ] **Step 1: Write the failing tests**

Write `tests/test_convert.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pytest tests/test_convert.py -v`
Expected: FAIL (ImportError — `create_dcmseg` not implemented)

- [ ] **Step 3: Write the implementation**

Write `src/pet_mr_ct/convert/nifti_to_dcmseg.py`. This extracts and refactors logic from the POC:

```python
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

    # RAS+ → LPS+: flip R→L and A→P
    mask_lps = mask_ras[::-1, ::-1, :]

    # (L, P, S) → (S, P, L) = (slice, row, col)
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
                "Dimension mismatch for %s: got %s, expected (%d, %d, %d) — skipping",
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pytest tests/test_convert.py -v`
Expected: All 5 tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/pet_mr_ct/convert/nifti_to_dcmseg.py tests/test_convert.py
git commit -m "feat: implement nifti_to_dcmseg with DICOM loading, RAS→LPS conversion, and tests"
```

---

## Task 5: Round-Trip Voxel Agreement Test

**Files:**
- Create: `tests/test_roundtrip.py`

**Reference:** POC `validate_roundtrip()` at `scripts/trial_dcmseg.py:295-308`, but this test goes further — it compares actual voxel data.

- [ ] **Step 1: Write the round-trip test**

```python
"""Tests for round-trip validation (export → reimport → compare)."""

from __future__ import annotations

from pathlib import Path

import highdicom as hd
import nibabel as nib
import numpy as np

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

    # Read back
    seg = hd.seg.segread(str(output))

    # Read back DICOM series to know total slice count
    import pydicom
    dcm_files = sorted(tmp_dicom_series.glob("*.dcm"))
    num_slices = len(dcm_files)
    ds0 = pydicom.dcmread(dcm_files[0])
    rows, cols = ds0.Rows, ds0.Columns

    for seg_num, name in enumerate(mask_names, start=1):
        # Extract segment from DICOM-SEG — reconstruct full volume
        # omit_empty_frames=True means we need to request reconstruction
        seg_array = seg.get_pixels_by_segment(
            segment_numbers=[seg_num],
            assert_spatial_locations_preserved=False,
        )
        # Shape: (frames, rows, cols, 1) — squeeze last dim
        seg_mask = seg_array.squeeze(axis=-1).astype(bool)

        # If omit_empty_frames removed slices, pad back to full volume
        if seg_mask.shape[0] < num_slices:
            full_mask = np.zeros((num_slices, rows, cols), dtype=bool)
            # Map frames back using segment frame info
            # For simple comparison, re-create from source NIfTI shape
            # and compare only non-empty region
            pass  # fall through to NIfTI-based comparison below

        # Load original NIfTI for comparison
        nii = nib.load(masks[name])
        canonical = nib.as_closest_canonical(nii)
        original_ras = np.asarray(canonical.dataobj, dtype=np.float32)
        original_lps = original_ras[::-1, ::-1, :]
        original_dicom = np.transpose(original_lps, (2, 1, 0))

        if canonical.affine[2, 2] < 0:
            original_dicom = original_dicom[::-1, :, :]

        original_mask = original_dicom > 0.5

        # Extract only the slices that have data in either mask
        # This handles omit_empty_frames correctly
        nonempty_slices = np.where(
            original_mask.any(axis=(1, 2)) | (
                seg_mask.any(axis=(1, 2)) if seg_mask.shape[0] == num_slices
                else np.zeros(num_slices, dtype=bool)
            )
        )[0]

        if seg_mask.shape[0] == num_slices:
            # Full volume returned — direct comparison
            total_voxels = original_mask.size
            matching = np.sum(seg_mask == original_mask)
        else:
            # Sparse frames — compare only non-empty slices from original
            original_nonempty = original_mask[nonempty_slices]
            total_voxels = original_nonempty.size
            matching = np.sum(seg_mask[:len(nonempty_slices)] == original_nonempty)

        agreement = matching / total_voxels

        assert agreement > 0.99, (
            f"Voxel agreement for {name}: {agreement:.4f} "
            f"({matching}/{total_voxels})"
        )
```

- [ ] **Step 2: Run to verify it passes**

Run: `pytest tests/test_roundtrip.py -v`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add tests/test_roundtrip.py
git commit -m "test: add voxel-level round-trip agreement test for DICOM-SEG"
```

---

## Task 6: Stubs and Module Init

**Files:**
- Modify: `src/pet_mr_ct/convert/dcm_to_nifti.py`
- Modify: `src/pet_mr_ct/convert/nifti_to_rtstruct.py`
- Modify: `src/pet_mr_ct/convert/__init__.py`

- [ ] **Step 1: Write the stubs**

Write `src/pet_mr_ct/convert/dcm_to_nifti.py`:

```python
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
```

Write `src/pet_mr_ct/convert/nifti_to_rtstruct.py`:

```python
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
```

- [ ] **Step 2: Write the module __init__.py**

Write `src/pet_mr_ct/convert/__init__.py`:

```python
"""Convert module: DICOM ↔ NIfTI ↔ DICOM-SEG / RT-STRUCT."""

from pet_mr_ct.convert.dcm_to_nifti import convert_dcm_to_nifti
from pet_mr_ct.convert.nifti_to_dcmseg import create_dcmseg
from pet_mr_ct.convert.nifti_to_rtstruct import create_rtstruct
from pet_mr_ct.convert.terminology import (
    SegmentInfo,
    get_segment_description,
    load_terminology,
)

__all__ = [
    "convert_dcm_to_nifti",
    "create_dcmseg",
    "create_rtstruct",
    "load_terminology",
    "get_segment_description",
    "SegmentInfo",
]
```

- [ ] **Step 3: Verify imports work**

Run: `python -c "from pet_mr_ct.convert import create_dcmseg, create_rtstruct, convert_dcm_to_nifti, load_terminology; print('OK')"`
Expected: `OK`

- [ ] **Step 4: Run full test suite**

Run: `pytest tests/ -v`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/pet_mr_ct/convert/
git commit -m "feat: add convert module stubs and public API re-exports"
```

---

## Task 7: CLI Commands

**Files:**
- Modify: `src/pet_mr_ct/cli.py`

- [ ] **Step 1: Write the CLI**

Write `src/pet_mr_ct/cli.py`:

```python
"""CLI entry point for the PET_MR_CT pipeline."""

from __future__ import annotations

from pathlib import Path

import click


@click.group()
@click.version_option()
def main() -> None:
    """DICOM-SEG pipeline: AI segmentation to MIM/3D Slicer round-trip."""


@main.group()
def convert() -> None:
    """Convert between DICOM, NIfTI, DICOM-SEG, and RT-STRUCT formats."""


@convert.command()
@click.option(
    "--source-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Directory containing source DICOM series.",
)
@click.option(
    "--mask-dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Directory containing NIfTI segmentation masks (.nii.gz or .nii).",
)
@click.option(
    "--tool",
    required=True,
    type=str,
    help="Segmentation tool name — must match a configs/segment_metadata_{tool}.json file (e.g. totalseg, moose).",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output DICOM-SEG file path.",
)
def dcmseg(source_dir: Path, mask_dir: Path, tool: str, output: Path) -> None:
    """Create DICOM-SEG from NIfTI masks + source DICOM series."""
    from pet_mr_ct.convert.nifti_to_dcmseg import create_dcmseg

    # Discover NIfTI masks
    masks: dict[str, Path] = {}
    for ext in ("*.nii.gz", "*.nii"):
        for f in mask_dir.glob(ext):
            stem = f.name.replace(".nii.gz", "").replace(".nii", "")
            if stem not in masks:  # .nii.gz takes priority
                masks[stem] = f

    if not masks:
        raise click.ClickException(f"No .nii.gz or .nii files found in {mask_dir}")

    click.echo(f"Found {len(masks)} mask(s) in {mask_dir}")
    result = create_dcmseg(source_dir, masks, tool, output)
    click.echo(f"DICOM-SEG written to {result}")


@convert.command()
@click.option("--source-dir", required=True, type=click.Path(path_type=Path))
@click.option("--mask-dir", required=True, type=click.Path(path_type=Path))
@click.option("--output", required=True, type=click.Path(path_type=Path))
def rtstruct(source_dir: Path, mask_dir: Path, output: Path) -> None:
    """Create RT-STRUCT from NIfTI masks (not yet implemented)."""
    click.echo("RT-STRUCT creation is not yet implemented.")
    click.echo("Use 'pet-mr-ct convert dcmseg' for DICOM-SEG output.")
    raise SystemExit(1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Verify CLI help works**

Run: `python -m pet_mr_ct.cli --help`
Expected: Shows main help with `convert` group

Run: `python -m pet_mr_ct.cli convert --help`
Expected: Shows `dcmseg` and `rtstruct` subcommands

Run: `python -m pet_mr_ct.cli convert dcmseg --help`
Expected: Shows `--source-dir`, `--mask-dir`, `--tool`, `--output` options

- [ ] **Step 3: Run full test suite**

Run: `pytest tests/ -v`
Expected: All tests still PASS

- [ ] **Step 4: Commit**

```bash
git add src/pet_mr_ct/cli.py
git commit -m "feat: add CLI convert group with dcmseg and rtstruct subcommands"
```

---

## Task 8: Final Integration Verification

- [ ] **Step 1: Run full test suite with coverage**

Run: `pytest tests/ -v --tb=short`
Expected: All tests PASS (test_terminology: 6, test_convert: 5, test_roundtrip: 1 = 12 total)

- [ ] **Step 2: Run ruff linter**

Run: `ruff check src/pet_mr_ct/convert/ tests/`
Expected: No errors (or fix any that appear)

- [ ] **Step 3: Test CLI end-to-end with real test data**

Run using the existing test data from `scripts/generate_test_data.py`:

```bash
python -m pet_mr_ct.cli convert dcmseg \
  --source-dir data/test/dicom_ct \
  --mask-dir data/test/segmentations \
  --tool totalseg \
  --output data/test/output/cli_test.dcm
```

Expected: DICOM-SEG file created at `data/test/output/cli_test.dcm`

- [ ] **Step 4: Commit any final fixes**

If ruff or integration test required fixes, commit them:

```bash
git add -u
git commit -m "fix: address linting and integration test issues"
```
