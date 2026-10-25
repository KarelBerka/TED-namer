"""Unit tests for rule-based and 3D-geometric CATH classification and 'T'-prefixed estimated levels."""

from ted_namer.cath_classifier import (
    classify_cath_from_sse,
    classify_cath_with_rules,
    format_t_prefixed_cath,
    infer_cath_classification,
    is_estimated_cath_code,
    strip_t_prefixes,
)
from ted_namer.models import OverlapCategory, OverlappingAnnotation, TedDomainMember


def test_format_and_strip_t_prefixes() -> None:
    code = format_t_prefixed_cath(["3", "40", "50", "100"], known_depth=1)
    assert code == "3.T40.T50.T100"
    assert is_estimated_cath_code(code) is True
    assert strip_t_prefixes(code) == "3.40.50.100"

    arch_direct = format_t_prefixed_cath(["3", "40", "50", "10"], known_depth=2)
    assert arch_direct == "3.40.T50.T10"

    topo_direct = format_t_prefixed_cath(["2", "60", "40", "10"], known_depth=3)
    assert topo_direct == "2.60.40.T10"
    assert is_estimated_cath_code(topo_direct) is True
    assert strip_t_prefixes(topo_direct) == "2.60.40.10"

    official = format_t_prefixed_cath(["2", "60", "40", "720"], known_depth=4)
    assert official == "2.60.40.720"
    assert is_estimated_cath_code(official) is False


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
    # High-IoU G3DSA match directly establishes C.A.T (3.40.640) and prefixes H (.T10)
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
    assert formatted == "3.40.640.T10"
    assert numeric == "3.40.640.10"
    assert source == "g3dsa_overlap"

    # Partial-IoU G3DSA match keeps only C direct and prefixes A.T.H
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
    assert formatted_p == "3.T40.T640.T10"
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
    assert formatted == "2.60.40.T10"
    assert numeric == "2.60.40.10"
    assert source == "ted_t_plus_sse"


def test_direct_rules_for_architecture_and_topology() -> None:
    # 1. Single alpha-helix coiled-coil rule -> direct 1.20.5 (depth=3 -> 1.20.5.T10)
    coiled_coil = TedDomainMember(
        ted_id="TED01",
        uniprot_acc="P00001",
        chopping="1-80",
        nres_domain=80,
        num_helix=1,
        num_strand=0,
        norm_rg=0.48,
    )
    f_cc, _, rule_cc = infer_cath_classification(coiled_coil, [])
    assert f_cc == "1.20.5.T10"
    assert rule_cc == "rule_alpha_coiled_coil_1.20.5"

    # 2. Classic (beta/alpha)8 TIM barrel rule -> direct 3.20.20 (depth=3 -> 3.20.20.T10)
    tim_dom = TedDomainMember(
        ted_id="TED01",
        uniprot_acc="P00002",
        chopping="1-250",
        nres_domain=250,
        num_helix=8,
        num_strand=8,
    )
    f_tim, _, rule_tim = infer_cath_classification(tim_dom, [])
    assert f_tim == "3.20.20.T10"
    assert rule_tim == "rule_sse_tim_barrel_3.20.20"

    # 3. Class 3 Rossmann-like SSE rule without 3D PDB -> direct C.A (3.40) and T-prefixed Topology/Superfamily (3.40.T50.T10)
    ab_dom = TedDomainMember(
        ted_id="TED01",
        uniprot_acc="P00003",
        chopping="1-180",
        nres_domain=180,
        num_helix=5,
        num_strand=6,
    )
    f_ab, num_ab, src_ab = infer_cath_classification(ab_dom, [])
    assert f_ab == "3.40.T50.T10"
    assert num_ab == "3.40.50.10"
    assert src_ab == "rule_sse_3layer_rossmann_3.40"

    # 4. Class 2 Beta sandwich SSE rule -> direct C.A (2.60) and T-prefixed T/H (2.60.T40.T10)
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
    f_beta, _, _ = infer_cath_classification(beta_dom, [])
    assert f_beta == "2.60.T40.T10"
