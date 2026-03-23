# PET_MR_CT — DICOM-SEG Pipeline

## What this project does

Python pipeline that takes AI segmentation output (TotalSegmentator, MOOSE, Merlin) and converts it to DICOM-SEG or RT-STRUCT for import into MIM 7.3.7 and 3D Slicer. Replaces an existing PowerShell workflow (nii2dcm.ps1) that used 7-Zip + XMedCon + DCMTK storescu and produced plain DICOM images (losing all segment labels).

## Architecture

The pipeline has 4 stages:
1. **convert**: DICOM ↔ NIfTI ↔ DICOM-SEG / RT-STRUCT. Core module. Uses highdicom for DICOM-SEG creation, rt-utils for RT-STRUCT fallback.
2. **segment**: Thin wrappers around TotalSegmentator and MOOSE CLIs. Tool-agnostic — accepts any NIfTI label map.
3. **register**: Cross-timepoint registration via ANTsPy/SimpleITK. Two strategies: register-then-segment (rigid organs) vs. segment-then-match (bowel/mesentery).
4. **network**: pynetdicom C-STORE to push DICOM-SEG/RT-STRUCT to MIM PACS.

## Key constraints

- **DICOM-SEG must preserve spatial registration** to the source DICOM series: same Frame of Reference UID, same image geometry. Use highdicom's `Segmentation` class with the source `Dataset` objects.
- **Every segment needs SNOMED CT coded terminology.** Map files are in `configs/segment_metadata_*.json`. Use CID 7151 (Segmentation Property Category) and CID 7152/7153 (Segmentation Property Type).
- **MIM 7.3.7 compatibility is uncertain for large DICOM-SEG** (>50 segments). Build with the option to split into per-region DICOM-SEG files or fall back to RT-STRUCT.
- **RT-STRUCT is the MIM fallback path.** MIM handles RT-STRUCT natively. Use rt-utils to convert NIfTI masks → RT-STRUCT contours.
- Pipeline must be tool-agnostic for the segmentation source. Don't hardcode assumptions about TotalSegmentator vs. MOOSE label conventions.

## Code conventions

- Python 3.11+, type hints everywhere
- Ruff for linting (line length 100)
- pytest for tests
- Click for CLI
- Use pathlib.Path, not os.path
- Docstrings: Google style
- All DICOM operations go through pydicom; never shell out to dcmtk or other binaries
- NIfTI I/O through nibabel only (handles .nii.gz natively, no 7-zip)

## Testing

- `pytest tests/` to run all tests
- Round-trip test: create DICOM-SEG from known NIfTI mask + source DICOM, read it back, compare voxel-level agreement
- Test with small fixtures (3-5 segments), not full TotalSegmentator output

## Target platforms

- Development: macOS / Linux / WSL2
- Deployment: WashU RIS HPC (Docker/Singularity)
- Consumers: MIM 7.3.7 (Windows), 3D Slicer (cross-platform)

## MIM network config (for DICOM send)

- AET: CCIR-e7RL2
- Called AET: CCIR_MIMV_2
- Host: 10.25.48.20
- Port: 4008
