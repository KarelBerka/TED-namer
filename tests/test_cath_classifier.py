"""Unit tests for rule-based and 3D-geometric CATH classification and strict hierarchical '?' wildcards."""

from ted_namer.cath_classifier import (
    classify_cath_from_sse,
    classify_cath_with_rules,
    format_hierarchical_cath,
    format_t_prefixed_cath,
    infer_cath_classification,
    is_estimated_cath_code,
    is_novel_fold_source,
    strip_t_prefixes,
)
from ted_namer.models import OverlapCategory, OverlappingAnnotation, TedDomainMember


def test_format_and_strip_wildcards_and_legacy_t_prefixes() -> None:
    code = format_hierarchical_cath(["3", "40", "50", "100"], known_depth=1)
    assert code == "3.?.?.?"
    assert is_estimated_cath_code(code) is True
    assert strip_t_prefixes(code) == "3"

    arch_direct = format_t_prefixed_cath(["3", "40", "50", "10"], known_depth=2)
    assert arch_direct == "3.40.?.?"
    assert strip_t_prefixes(arch_direct) == "3.40"

    topo_direct = format_t_prefixed_cath(["2", "60", "40", "10"], known_depth=3)
    assert topo_direct == "2.60.40.?"
    assert is_estimated_cath_code(topo_direct) is True
    assert strip_t_prefixes(topo_direct) == "2.60.40"

    official = format_t_prefixed_cath(["2", "60", "40", "720"], known_depth=4)
    assert official == "2.60.40.720"
    assert is_estimated_cath_code(official) is False
    assert strip_t_prefixes(official) == "2.60.40.720"

    # Legacy 'T'-prefixed strings remain backwards compatible
    assert is_estimated_cath_code("3.T40.T50.T10") is True
    assert strip_t_prefixes("3.T40.T50.T10") == "3.40.50.10"


def test_novel_fold_source_helper() -> None:
    assert is_novel_fold_source("rule_sse_3layer_rossmann_3.40") is True
    assert is_novel_fold_source("rule_3d_tim_barrel_3.20.20") is True
    assert is_novel_fold_source("ted_official") is False
    assert is_novel_fold_source("g3dsa_overlap") is False
    assert is_novel_fold_source("profile_cath_map") is False
    assert is_novel_fold_source("ted_t_plus_g3dsa") is False


def test_infer_cath_from_g3dsa_overlap_when_ted_unassigned() -> None:
    member = TedDomainMember(
        ted_id="AF-B7Z3M3-F1-model_v4_TED01",
        uniprot_acc="B7Z3M3",
        chopping="10-250",
        nres_domain=241,
        num_helix=8,
        num_strand=7,
        cath_label="-",
        cath_assignment_level="-",
    )
    # High-IoU G3DSA match (>=0.65 IoU, >=0.75 TED cov) -> direct 4-level code
    annots_high = [
        OverlappingAnnotation(
            source_db="cathgene3d",
            accession="G3DSA:3.40.640.10",
            name="Pyridoxal phosphate-dependent transferase",
            annotation_type="homologous_superfamily",
            chopping="15-245",
            iou=0.92,
            ted_coverage=0.95,
            annot_coverage=0.96,
            overlap_category=OverlapCategory.EXACT_MATCH,
        )
    ]
    formatted, numeric, source = infer_cath_classification(member, annots_high)
    assert formatted == "3.40.640.10"
    assert numeric == "3.40.640.10"
    assert source == "g3dsa_overlap"
    assert is_novel_fold_source(source) is False

    # Moderate-IoU G3DSA match establishes C.A.T and marks H as '?' wildcard
    annots_mod = [
        OverlappingAnnotation(
            source_db="cathgene3d",
            accession="G3DSA:3.40.640.10",
            name="Pyridoxal phosphate-dependent transferase",
            annotation_type="homologous_superfamily",
            chopping="15-160",
            iou=0.50,
            ted_coverage=0.60,
            annot_coverage=0.90,
            overlap_category=OverlapCategory.PARTIAL_OVERLAP,
        )
    ]
    formatted_m, numeric_m, _ = infer_cath_classification(member, annots_mod)
    assert formatted_m == "3.40.640.?"
    assert numeric_m == "3.40.640.10"

    # Low-IoU G3DSA match keeps only C direct and uses '?' for A.T.H
    annots_partial = [
        OverlappingAnnotation(
            source_db="cathgene3d",
            accession="G3DSA:3.40.640.10",
            name="Pyridoxal phosphate-dependent transferase",
            annotation_type="homologous_superfamily",
            chopping="15-120",
            iou=0.35,
            ted_coverage=0.50,
            annot_coverage=0.90,
            overlap_category=OverlapCategory.PARTIAL_OVERLAP,
        )
    ]
    formatted_p, numeric_p, _ = infer_cath_classification(member, annots_partial)
    assert formatted_p == "3.?.?.?"
    assert numeric_p == "3.40.640.10"


def test_infer_cath_when_ted_is_t_level() -> None:
    member = TedDomainMember(
        ted_id="AF-P0DP43-F1-model_v4_TED01",
        uniprot_acc="P0DP43",
        chopping="20-130",
        nres_domain=111,
        num_helix=0,
        num_strand=8,
        cath_label="2.60.40",
        cath_assignment_level="T",
    )
    formatted, numeric, source = infer_cath_classification(member, [])
    assert formatted == "2.60.40.?"
    assert numeric == "2.60.40.10"
    assert source == "ted_t_plus_sse"
    assert is_novel_fold_source(source) is False


def test_direct_rules_for_architecture_and_topology() -> None:
    # 1. Single alpha-helix coiled-coil rule -> Class 1 direct, '1.?.?.?' with fold hint '1.20.5.10'
    coiled_coil = TedDomainMember(
        ted_id="TED01",
        uniprot_acc="P00001",
        chopping="1-80",
        nres_domain=80,
        num_helix=1,
        num_strand=0,
        norm_rg=0.48,
    )
    f_cc, num_cc, rule_cc = infer_cath_classification(coiled_coil, [])
    assert f_cc == "1.?.?.?"
    assert num_cc == "1.20.5.10"
    assert rule_cc == "rule_alpha_coiled_coil_1.20.5"
    assert is_novel_fold_source(rule_cc) is True

    # 2. Class 3 Rossmann-like SSE rule without 3D PDB -> '3.?.?.?' with fold hint '3.40.50.10'
    ab_dom = TedDomainMember(
        ted_id="TED01",
        uniprot_acc="P00003",
        chopping="1-180",
        nres_domain=180,
        num_helix=5,
        num_strand=6,
    )
    f_ab, num_ab, src_ab = infer_cath_classification(ab_dom, [])
    assert f_ab == "3.?.?.?"
    assert num_ab == "3.40.50.10"
    assert src_ab == "rule_sse_3layer_rossmann_3.40"

    # 3. Class 2 Beta sandwich SSE rule -> '2.?.?.?' with fold hint '2.60.40.10'
    beta_dom = TedDomainMember(
        ted_id="TED01",
        uniprot_acc="P00004",
        chopping="1-110",
        nres_domain=110,
        num_helix=0,
        num_strand=8,
        packing_density=11.5,
    )
    c2, a2, t2, h2 = classify_cath_from_sse(beta_dom)
    assert (c2, a2, t2, h2) == ("2", "60", "40", "10")
    f_beta, num_beta, _ = infer_cath_classification(beta_dom, [])
    assert f_beta == "2.?.?.?"
    assert num_beta == "2.60.40.10"

    # 4. Class 4 Few Secondary Structures -> '4.10.?.?' (known_depth = 2)
    few_sse = TedDomainMember(
        ted_id="TED02",
        uniprot_acc="P00005",
        chopping="1-30",
        nres_domain=30,
        num_helix=0,
        num_strand=0,
        num_turn=3,
    )
    f_few, num_few, src_few = infer_cath_classification(few_sse, [])
    assert f_few == "4.10.?.?"
    assert num_few == "4.10.10.10"
    assert is_novel_fold_source(src_few) is True
