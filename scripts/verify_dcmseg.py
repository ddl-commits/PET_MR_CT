#!/usr/bin/env python3
"""Quick visual verification of DICOM-SEG output.

Generates a multi-panel figure showing:
  - Source CT slices with segment overlays
  - Per-segment voxel counts
  - Round-trip agreement stats

Usage:
    python scripts/verify_dcmseg.py [--dcmseg PATH] [--source-dir PATH]
"""

from __future__ import annotations

from pathlib import Path

import highdicom as hd
import matplotlib.pyplot as plt
import numpy as np
import pydicom

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DCMSEG = PROJECT_ROOT / "data" / "test" / "output" / "cli_test.dcm"
DEFAULT_SOURCE = PROJECT_ROOT / "data" / "test" / "dicom_ct"
OUTPUT_DIR = PROJECT_ROOT / "data" / "test" / "output"


def load_ct_volume(dcm_dir: Path) -> tuple[np.ndarray, list[pydicom.Dataset]]:
    """Load CT DICOM series as a 3D volume."""
    dcm_files = sorted(dcm_dir.glob("*.dcm"))
    datasets = [pydicom.dcmread(f) for f in dcm_files]
    datasets.sort(key=lambda ds: float(ds.ImagePositionPatient[2]))

    volume = np.stack([ds.pixel_array for ds in datasets], axis=0)
    return volume, datasets


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Verify DICOM-SEG visually")
    parser.add_argument("--dcmseg", type=Path, default=DEFAULT_DCMSEG)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args()

    print(f"Loading DICOM-SEG: {args.dcmseg}")
    seg = hd.seg.segread(str(args.dcmseg))
    n_segments = seg.number_of_segments

    print(f"Loading source CT: {args.source_dir}")
    ct_vol, ct_datasets = load_ct_volume(args.source_dir)
    n_slices = ct_vol.shape[0]

    # Get segment labels and colors
    cmap = plt.cm.get_cmap("tab20", n_segments)
    labels = []
    for i in range(n_segments):
        desc = seg.get_segment_description(i + 1)
        labels.append(desc.segment_label)

    # Extract all segment masks
    print(f"Extracting {n_segments} segments...")
    all_masks = {}
    voxel_counts = {}
    for seg_num in range(1, n_segments + 1):
        mask = seg.get_pixels_by_source_instance(
            source_sop_instance_uids=[ds.SOPInstanceUID for ds in ct_datasets],
            segment_numbers=[seg_num],
            assert_missing_frames_are_empty=True,
        )
        mask = mask.squeeze(axis=-1).astype(bool)
        all_masks[seg_num] = mask
        voxel_counts[seg_num] = int(mask.sum())

    # --- Figure 1: Segment overlay on CT slices ---
    # Pick 4 evenly spaced slices that have segments
    has_seg = np.zeros(n_slices, dtype=bool)
    for mask in all_masks.values():
        has_seg |= mask.any(axis=(1, 2))
    seg_slices = np.where(has_seg)[0]

    if len(seg_slices) >= 4:
        pick_idx = np.linspace(0, len(seg_slices) - 1, 4, dtype=int)
        show_slices = seg_slices[pick_idx]
    else:
        show_slices = seg_slices[:4] if len(seg_slices) > 0 else [n_slices // 2]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("DICOM-SEG Verification: CT + Segment Overlays", fontsize=14)

    for ax_idx, sl in enumerate(show_slices):
        ax = axes[ax_idx // 2, ax_idx % 2]
        ax.imshow(ct_vol[sl], cmap="gray", vmin=-200, vmax=400)

        # Overlay each segment with a different color
        overlay = np.zeros((*ct_vol[sl].shape, 4))
        for seg_num, mask in all_masks.items():
            if mask[sl].any():
                color = cmap(seg_num - 1)
                region = mask[sl]
                overlay[region, :3] = color[:3]
                overlay[region, 3] = 0.4

        ax.imshow(overlay)
        ax.set_title(f"Slice {sl} (z={float(ct_datasets[sl].ImagePositionPatient[2]):.1f}mm)")
        ax.axis("off")

    plt.tight_layout()
    out1 = OUTPUT_DIR / "verify_overlays.png"
    fig.savefig(out1, dpi=150, bbox_inches="tight")
    print(f"Saved: {out1}")

    # --- Figure 2: Voxel count bar chart ---
    fig2, ax2 = plt.subplots(figsize=(14, 6))

    # Sort by voxel count, show top 20
    sorted_segs = sorted(voxel_counts.items(), key=lambda x: x[1], reverse=True)
    top_n = min(20, len(sorted_segs))
    top_segs = sorted_segs[:top_n]

    seg_labels = [labels[s - 1] for s, _ in top_segs]
    counts = [c for _, c in top_segs]
    colors = [cmap(s - 1) for s, _ in top_segs]

    bars = ax2.barh(range(top_n), counts, color=colors)
    ax2.set_yticks(range(top_n))
    ax2.set_yticklabels(seg_labels, fontsize=9)
    ax2.set_xlabel("Voxel Count")
    ax2.set_title(f"Top {top_n} Segments by Voxel Count (of {n_segments} total)")
    ax2.invert_yaxis()

    for bar, count in zip(bars, counts):
        ax2.text(bar.get_width() + max(counts) * 0.01, bar.get_y() + bar.get_height() / 2,
                 f"{count:,}", va="center", fontsize=8)

    plt.tight_layout()
    out2 = OUTPUT_DIR / "verify_voxel_counts.png"
    fig2.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"Saved: {out2}")

    # --- Summary stats ---
    print("\n" + "=" * 60)
    print("DICOM-SEG Verification Summary")
    print("=" * 60)
    print(f"  File: {args.dcmseg}")
    print(f"  File size: {args.dcmseg.stat().st_size / (1024*1024):.2f} MB")
    print(f"  Total segments: {n_segments}")
    print(f"  CT slices: {n_slices}")
    print(f"  CT dimensions: {ct_vol.shape[1]}x{ct_vol.shape[2]}")
    print(f"  Segments with voxels: {sum(1 for c in voxel_counts.values() if c > 0)}")
    print(f"  Empty segments: {sum(1 for c in voxel_counts.values() if c == 0)}")
    total_voxels = sum(voxel_counts.values())
    print(f"  Total segmented voxels: {total_voxels:,}")
    print(f"\n  Top 5 segments:")
    for seg_num, count in sorted_segs[:5]:
        print(f"    {labels[seg_num-1]:30s} {count:>10,} voxels")
    print(f"\n  Figures saved to: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
