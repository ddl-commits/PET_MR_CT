# PET_MR_CT

DICOM-SEG pipeline for PET/MR/CT: AI segmentation (TotalSegmentator, MOOSE, Merlin) → DICOM-SEG/RT-STRUCT → MIM 7.3.7 / 3D Slicer.

## Install

```bash
pip install -e ".[dev]"
```

For segmentation tools (requires GPU):
```bash
pip install -e ".[segment]"
```

## Usage

```bash
pet-mr-ct --help
```
