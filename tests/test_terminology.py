"""Tests for SNOMED CT terminology mappings."""

from __future__ import annotations

import highdicom as hd
import pytest

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
