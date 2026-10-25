"""Unit tests for continuous and discontinuous TED chopping parsing, residue IoU/coverage, and relative position."""

from ted_namer.models import OverlapCategory
from ted_namer.overlap import (
    compute_overlap_metrics,
    compute_relative_position,
    format_segments,
    parse_chopping,
)


def test_parse_and_format_continuous_and_discontinuous_chopping() -> None:
    segs = parse_chopping("109-287")
    assert len(segs) == 1
    assert segs[0].start == 109
    assert segs[0].end == 287
    assert segs[0].length == 179

    multi = parse_chopping("140-210_15-65")
    assert len(multi) == 2
    assert format_segments(multi) == "15-65_140-210"


def test_exact_match_overlap() -> None:
    iou, ted_cov, annot_cov, cat = compute_overlap_metrics("109-287", "109-288")
    assert iou > 0.99
    assert ted_cov == 1.0
    assert annot_cov > 0.99
    assert cat == OverlapCategory.EXACT_MATCH


def test_discontinuous_overlap_metrics() -> None:
    # TED domain has two segments: 10-50 (41 res) and 100-149 (50 res) = 91 residues
    # Annotation covers 10-50 (41 res) and 100-140 (41 res) = 82 residues
    iou, ted_cov, annot_cov, cat = compute_overlap_metrics("10-50_100-149", "10-50_100-140")
    assert round(iou, 4) == round(82 / 91, 4)
    assert round(ted_cov, 4) == round(82 / 91, 4)
    assert annot_cov == 1.0
    assert cat == OverlapCategory.EXACT_MATCH


def test_sub_domain_and_supra_domain_classification() -> None:
    # Small motif (120-140, 21 res) inside a large TED domain (100-300, 201 res)
    iou_sub, ted_cov_sub, annot_cov_sub, cat_sub = compute_overlap_metrics("100-300", "120-140")
    assert annot_cov_sub == 1.0
    assert ted_cov_sub < 0.15
    assert cat_sub == OverlapCategory.SUB_DOMAIN

    # Whole-protein family (1-393) spanning a single TED domain (109-287, 179 res)
    iou_sup, ted_cov_sup, annot_cov_sup, cat_sup = compute_overlap_metrics("109-287", "1-393")
    assert ted_cov_sup == 1.0
    assert annot_cov_sup < 0.50
    assert cat_sup == OverlapCategory.SUPRA_DOMAIN


def test_disjoint_domains() -> None:
    iou, ted_cov, annot_cov, cat = compute_overlap_metrics("10-80", "120-200")
    assert iou == 0.0
    assert cat == OverlapCategory.DISJOINT


def test_compute_relative_position() -> None:
    assert compute_relative_position("10-120", 500) == "N-terminal"
    assert compute_relative_position("200-300", 500) == "central"
    assert compute_relative_position("380-490", 500) == "C-terminal"
    assert compute_relative_position("5-195", 200) == "full-length"
