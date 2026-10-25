"""Residue boundary parsing, multi-segment chopping overlap (IoU/coverage), and positional analysis."""

from __future__ import annotations

import re
from typing import Iterable

from ted_namer.models import OverlapCategory, SegmentRange

_SEGMENT_RE = re.compile(r"(\d+)\s*-\s*(\d+)")


def parse_chopping(chopping: str) -> list[SegmentRange]:
    """Parse a TED or InterPro chopping string (e.g. '109-287' or '15-65_140-210') into SegmentRanges."""
    if not chopping:
        return []
    segments: list[SegmentRange] = []
    for match in _SEGMENT_RE.finditer(chopping):
        start = int(match.group(1))
        end = int(match.group(2))
        if end >= start > 0:
            segments.append(SegmentRange(start=start, end=end))
    return sorted(segments, key=lambda s: (s.start, s.end))


def format_segments(segments: Iterable[SegmentRange]) -> str:
    """Format a sequence of SegmentRanges into canonical TED chopping notation ('15-65_140-210')."""
    sorted_segs = sorted(segments, key=lambda s: (s.start, s.end))
    return "_".join(f"{s.start}-{s.end}" for s in sorted_segs if s.end >= s.start)


def segments_to_residue_set(segments: Iterable[SegmentRange]) -> set[int]:
    """Convert one or more segments into a 1-based set of residue indices."""
    residues: set[int] = set()
    for seg in segments:
        if seg.end >= seg.start:
            residues.update(range(seg.start, seg.end + 1))
    return residues


def compute_overlap_metrics(
    ted_chopping: str | Iterable[SegmentRange],
    annot_chopping: str | Iterable[SegmentRange],
) -> tuple[float, float, float, OverlapCategory]:
    """Compute (iou, ted_coverage, annot_coverage, overlap_category) between two domain choppings."""
    ted_segs = parse_chopping(ted_chopping) if isinstance(ted_chopping, str) else list(ted_chopping)
    annot_segs = parse_chopping(annot_chopping) if isinstance(annot_chopping, str) else list(annot_chopping)

    ted_res = segments_to_residue_set(ted_segs)
    annot_res = segments_to_residue_set(annot_segs)

    if not ted_res or not annot_res:
        return 0.0, 0.0, 0.0, OverlapCategory.DISJOINT

    inter = len(ted_res & annot_res)
    if inter == 0:
        return 0.0, 0.0, 0.0, OverlapCategory.DISJOINT

    union = len(ted_res | annot_res)
    iou = inter / union
    ted_cov = inter / len(ted_res)
    annot_cov = inter / len(annot_res)

    category = classify_overlap(iou=iou, ted_coverage=ted_cov, annot_coverage=annot_cov)
    return round(iou, 4), round(ted_cov, 4), round(annot_cov, 4), category


def classify_overlap(iou: float, ted_coverage: float, annot_coverage: float) -> OverlapCategory:
    """Classify the spatial relationship between a sequence annotation and a TED domain."""
    if iou >= 0.70 or (ted_coverage >= 0.80 and annot_coverage >= 0.80):
        return OverlapCategory.EXACT_MATCH
    if annot_coverage >= 0.80 and ted_coverage < 0.65:
        return OverlapCategory.SUB_DOMAIN
    if ted_coverage >= 0.80 and annot_coverage < 0.65:
        return OverlapCategory.SUPRA_DOMAIN
    if iou >= 0.10 or ted_coverage >= 0.15 or annot_coverage >= 0.25:
        return OverlapCategory.PARTIAL_OVERLAP
    return OverlapCategory.DISJOINT


def compute_relative_position(chopping: str | Iterable[SegmentRange], protein_length: int) -> str:
    """Determine whether a TED domain is N-terminal, central, C-terminal, or full-length within a protein."""
    segs = parse_chopping(chopping) if isinstance(chopping, str) else list(chopping)
    if not segs or protein_length <= 0:
        return "domain"

    res = segments_to_residue_set(segs)
    if not res:
        return "domain"

    coverage = len(res) / protein_length
    if coverage >= 0.75:
        return "full-length"

    midpoint = (min(res) + max(res)) / 2.0
    rel_mid = midpoint / protein_length

    if rel_mid <= 0.36:
        return "N-terminal"
    if rel_mid >= 0.64:
        return "C-terminal"
    return "central"
