"""Convert module: DICOM <-> NIfTI <-> DICOM-SEG / RT-STRUCT."""

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
