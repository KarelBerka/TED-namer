"""Unit tests for cluster evidence aggregation, synonym unification, TrEMBL noise filtering, and novel fold naming."""

from ted_namer.aggregator import aggregate_cluster_evidence
from ted_namer.models import (
    CathResolution,
    EnrichedTedDomain,
    OverlapCategory,
    OverlappingAnnotation,
    ProteinContext,
    SpecificityTier,
    TedDomainMember,
)


def test_synonym_unification_and_tier_ranking() -> None:
    """CDD, Pfam, and InterPro domain hits sharing an integrated IPR ID unite and outrank broad family hits."""
    m1 = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-P04637-F1-model_v4_TED01",
            uniprot_acc="P04637",
            chopping="109-287",
            nres_domain=179,
            plddt=95.1,
            cath_label="2.60.40.720",
            cath_assignment_level="H",
        ),
        cath_info=CathResolution(
            cath_label="2.60.40.720",
            assignment_level="H",
            superfamily_name="p53/RUNT-type transcription factor, DNA-binding domain superfamily",
            integrated_ipr="IPR012346",
            topology_code="2.60.40",
            topology_name="Immunoglobulin-like",
        ),
        protein_context=ProteinContext(
            uniprot_acc="P04637",
            recommended_name="Cellular tumor antigen p53",
            is_reviewed=True,
            protein_length=393,
            similarity_families=["p53 family"],
        ),
        overlapping_annotations=[
            OverlappingAnnotation(
                source_db="cdd",
                accession="cd08367",
                name="P53 DNA-binding domain",
                integrated_ipr="IPR011615",
                annotation_type="domain",
                chopping="109-288",
                iou=0.994,
                ted_coverage=1.0,
                annot_coverage=0.994,
                overlap_category=OverlapCategory.EXACT_MATCH,
            ),
            OverlappingAnnotation(
                source_db="pfam",
                accession="PF00870",
                name="P53 DNA-binding domain",
                integrated_ipr="IPR011615",
                annotation_type="domain",
                chopping="99-289",
                iou=0.937,
                ted_coverage=1.0,
                annot_coverage=0.937,
                overlap_category=OverlapCategory.EXACT_MATCH,
            ),
            OverlappingAnnotation(
                source_db="panther",
                accession="PTHR11447",
                name="CELLULAR TUMOR ANTIGEN P53",
                integrated_ipr="IPR002117",
                annotation_type="family",
                chopping="1-393",
                iou=0.455,
                ted_coverage=1.0,
                annot_coverage=0.455,
                overlap_category=OverlapCategory.SUPRA_DOMAIN,
            ),
        ],
        relative_position="central",
    )

    res = aggregate_cluster_evidence("p53_dbd_cluster", [m1])
    assert res.specificity_tier == SpecificityTier.SPECIFIC_DOMAIN
    assert "dna-binding domain" in res.consensus_name.lower()
    assert res.confidence_score >= 0.90
    # Verify CDD and Pfam were unified under the top candidate
    top = res.top_candidates[0]
    assert "cdd" in top.source_databases
    assert "pfam" in top.source_databases


def test_informative_swissprot_outranks_majority_uncharacterized_trembl() -> None:
    """Even when 3 out of 4 cluster members are uncharacterized/DUF TrEMBL proteins, the 1 annotated member wins."""
    unchar_members = [
        EnrichedTedDomain(
            domain=TedDomainMember(
                ted_id=f"AF-A0A00{i}-F1-model_v4_TED01",
                uniprot_acc=f"A0A00{i}",
                chopping="10-150",
                nres_domain=141,
                cath_label="-",
            ),
            protein_context=ProteinContext(
                uniprot_acc=f"A0A00{i}",
                recommended_name="Uncharacterized protein",
                is_reviewed=False,
                protein_length=300,
            ),
            overlapping_annotations=[
                OverlappingAnnotation(
                    source_db="pfam",
                    accession="PF09999",
                    name="Domain of unknown function (DUF9999)",
                    annotation_type="domain",
                    chopping="10-150",
                    iou=1.0,
                    ted_coverage=1.0,
                    annot_coverage=1.0,
                    overlap_category=OverlapCategory.EXACT_MATCH,
                )
            ],
            relative_position="N-terminal",
        )
        for i in range(1, 4)
    ]

    annotated_member = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-P12345-F1-model_v4_TED01",
            uniprot_acc="P12345",
            chopping="12-148",
            nres_domain=137,
            cath_label="-",
        ),
        protein_context=ProteinContext(
            uniprot_acc="P12345",
            recommended_name="CRISPR-associated endonuclease Cas9",
            is_reviewed=True,
            protein_length=300,
        ),
        overlapping_annotations=[
            OverlappingAnnotation(
                source_db="cdd",
                accession="cd99999",
                name="Cas9 REC lobe recognition domain",
                annotation_type="domain",
                chopping="10-150",
                iou=0.97,
                ted_coverage=1.0,
                annot_coverage=0.97,
                overlap_category=OverlapCategory.EXACT_MATCH,
            )
        ],
        relative_position="N-terminal",
    )

    res = aggregate_cluster_evidence("mixed_cluster", unchar_members + [annotated_member])
    assert "cas9 rec lobe recognition domain" in res.consensus_name.lower()
    assert res.specificity_tier == SpecificityTier.SPECIFIC_DOMAIN


def test_novel_fold_positional_fallback() -> None:
    """When a novel fold (cath_label='-') has no sequence domain annotations, it uses protein family + position."""
    m = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-Q99999-F1-model_v4_TED02",
            uniprot_acc="Q99999",
            chopping="400-520",
            nres_domain=121,
            cath_label="-",
        ),
        protein_context=ProteinContext(
            uniprot_acc="Q99999",
            recommended_name="ATP-dependent helicase HelX",
            is_reviewed=True,
            protein_length=530,
            similarity_families=["DEAD-box helicase family"],
        ),
        overlapping_annotations=[],
        relative_position="C-terminal",
    )

    res = aggregate_cluster_evidence("novel_fold_cluster", [m])
    assert res.specificity_tier == SpecificityTier.PROTEIN_FAMILY_CONTEXT
    assert "c-terminal domain" in res.consensus_name.lower()
    assert "helicase" in res.consensus_name.lower()


def test_curated_subdomain_lobe_and_integrated_ipr_override() -> None:
    """Verify multi-lobe sequence domains resolve to positional subdomains and integrated IPR overrides raw G3DSA PDB titles."""
    m = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-O13958-F1-model_v4_TED01",
            uniprot_acc="O13958",
            chopping="25-110",
            nres_domain=86,
            cath_label="3.30.200.20",
            cath_assignment_level="H",
        ),
        protein_context=ProteinContext(
            uniprot_acc="O13958",
            recommended_name="Serine/threonine-protein kinase",
            is_reviewed=True,
            protein_length=450,
        ),
        overlapping_annotations=[
            OverlappingAnnotation(
                source_db="pfam",
                accession="PF00069",
                name="Protein kinase domain",
                integrated_ipr="IPR000719",
                integrated_ipr_name="Protein kinase domain",
                annotation_type="domain",
                chopping="25-285",
                iou=0.33,
                ted_coverage=1.0,
                annot_coverage=0.33,
                overlap_category=OverlapCategory.SUPRA_DOMAIN,
            ),
            OverlappingAnnotation(
                source_db="cathgene3d",
                accession="G3DSA:3.30.200.20",
                name="Phosphorylase Kinase; domain 1",
                integrated_ipr="IPR011009",
                integrated_ipr_name="Protein kinase-like domain superfamily",
                annotation_type="homologous_superfamily",
                chopping="20-112",
                iou=0.91,
                ted_coverage=1.0,
                annot_coverage=0.91,
                overlap_category=OverlapCategory.EXACT_MATCH,
            ),
        ],
        relative_position="N-terminal",
    )

    res = aggregate_cluster_evidence("kinase_nlobe", [m])
    assert res.specificity_tier == SpecificityTier.SPECIFIC_DOMAIN
    assert res.consensus_name == "Protein kinase N-terminal subdomain"
    # Also verify that the superfamily candidate used the integrated InterPro name rather than 'Phosphorylase Kinase; domain 1'
    sf_cands = [c for c in res.top_candidates if c.specificity_tier == SpecificityTier.SUPERFAMILY]
    assert sf_cands
    assert sf_cands[0].canonical_name == "Protein kinase-like domain superfamily"
