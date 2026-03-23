"""SNOMED CT coded terminology for organ segment labels."""

from __future__ import annotations

import dataclasses
import functools
import json
import logging
from importlib import resources
from pathlib import Path

import highdicom as hd

logger = logging.getLogger(__name__)


def _get_configs_dir() -> Path:
    """Resolve the configs directory, handling both dev and installed layouts."""
    dev_path = Path(__file__).resolve().parent.parent.parent.parent / "configs"
    if dev_path.is_dir():
        return dev_path
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

    Results are cached -- repeated calls with the same tool name are free.

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
