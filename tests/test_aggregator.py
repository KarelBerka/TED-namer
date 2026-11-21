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


def test_unintegrated_legacy_g3dsa_overrides_and_kinase_lobe_competition() -> None:
    """Unintegrated legacy G3DSA titles map to curated superfamily names and never outscore clean Tier 1 kinase domains."""
    m = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-Q6DBX4-F1-model_v4_TED01",
            uniprot_acc="Q6DBX4",
            chopping="40-300",
            nres_domain=261,
            cath_label="3.30.200.20",
            cath_assignment_level="H",
        ),
        cath_info=CathResolution(
            cath_label="3.30.200.20",
            assignment_level="H",
            superfamily_name="Phosphorylase Kinase; domain 1",
            integrated_ipr=None,
        ),
        protein_context=ProteinContext(
            uniprot_acc="Q6DBX4",
            recommended_name="Serine/threonine-protein kinase",
            is_reviewed=True,
            protein_length=600,
        ),
        overlapping_annotations=[
            OverlappingAnnotation(
                source_db="pfam",
                accession="PF00069",
                name="Protein kinase domain",
                integrated_ipr="IPR000719",
                annotation_type="domain",
                chopping="35-295",
                iou=0.90,
                ted_coverage=0.92,
                annot_coverage=0.97,
                overlap_category=OverlapCategory.EXACT_MATCH,
            ),
            OverlappingAnnotation(
                source_db="cathgene3d",
                accession="G3DSA:3.30.200.20",
                name="Phosphorylase Kinase; domain 1",
                integrated_ipr=None,
                annotation_type="homologous_superfamily",
                chopping="38-150",
                iou=0.43,
                ted_coverage=0.43,
                annot_coverage=1.0,
                overlap_category=OverlapCategory.SUB_DOMAIN,
            ),
        ],
        relative_position="N-terminal",
    )
    res = aggregate_cluster_evidence("kinase_test", [m])
    assert res.consensus_name == "Protein kinase domain"
    assert res.specificity_tier == SpecificityTier.SPECIFIC_DOMAIN
    sf_cands = [c for c in res.top_candidates if c.specificity_tier == SpecificityTier.SUPERFAMILY]
    assert sf_cands
    assert sf_cands[0].canonical_name == "Protein kinase-like N-terminal lobe superfamily"

    # Also test standalone C-lobe 1.10.510.10 unintegrated override when only CATH H / G3DSA exists
    m_clobe = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-Q8IWB3-F1-model_v4_TED01",
            uniprot_acc="Q8IWB3",
            chopping="150-350",
            nres_domain=201,
            cath_label="1.10.510.10",
            cath_assignment_level="H",
        ),
        cath_info=CathResolution(
            cath_label="1.10.510.10",
            assignment_level="H",
            superfamily_name="Transferase(Phosphotransferase) domain 1",
            integrated_ipr=None,
        ),
        overlapping_annotations=[
            OverlappingAnnotation(
                source_db="cathgene3d",
                accession="G3DSA:1.10.510.10",
                name="Transferase(Phosphotransferase) domain 1",
                annotation_type="homologous_superfamily",
                chopping="150-345",
                iou=0.97,
                ted_coverage=0.97,
                annot_coverage=1.0,
                overlap_category=OverlapCategory.EXACT_MATCH,
            )
        ],
        relative_position="central",
    )
    res_c = aggregate_cluster_evidence("clobe_only", [m_clobe])
    assert res_c.consensus_name == "Protein kinase-like C-terminal catalytic lobe superfamily"
    assert res_c.specificity_tier == SpecificityTier.SUPERFAMILY


def test_annotation_type_strictness_hamap_ncbifam_pfam_family() -> None:
    """Whole-protein HAMAP/NCBIfam families go to Tier 3, while strict 1:1 Pfam families qualify for Tier 1."""
    m = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-A0A1G9Z131-F1-model_v4_TED02",
            uniprot_acc="A0A1G9Z131",
            chopping="182-316",
            nres_domain=135,
            cath_label="3.40.50.720",
        ),
        protein_context=ProteinContext(
            uniprot_acc="A0A1G9Z131",
            recommended_name="Bifunctional protein BioU",
            is_reviewed=False,
            protein_length=450,
        ),
        overlapping_annotations=[
            OverlappingAnnotation(
                source_db="hamap",
                accession="MF_01920",
                name="8-amino-7-oxononanoate synthase BioU [bioU]",
                annotation_type="family",
                chopping="1-450",
                iou=0.30,
                ted_coverage=1.0,
                annot_coverage=0.30,
                overlap_category=OverlapCategory.SUPRA_DOMAIN,
            ),
            OverlappingAnnotation(
                source_db="ncbifam",
                accession="TIGR01234",
                name="Biotin biosynthesis enzyme BioU",
                annotation_type="family",
                chopping="1-450",
                iou=0.30,
                ted_coverage=1.0,
                annot_coverage=0.30,
                overlap_category=OverlapCategory.SUPRA_DOMAIN,
            ),
            OverlappingAnnotation(
                source_db="cathgene3d",
                accession="G3DSA:3.40.50.720",
                name="NAD(P)-binding Rossmann-like Domain",
                integrated_ipr="IPR036291",
                integrated_ipr_name="NAD(P)-binding domain superfamily",
                annotation_type="homologous_superfamily",
                chopping="180-318",
                iou=0.96,
                ted_coverage=1.0,
                annot_coverage=0.96,
                overlap_category=OverlapCategory.EXACT_MATCH,
            ),
        ],
        relative_position="C-terminal",
    )
    res = aggregate_cluster_evidence("hamap_family_test", [m])
    # Superfamily (Tier 2) must outrank whole-protein HAMAP/NCBIfam families (now Tier 3 with '[bioU]' stripped)
    assert res.specificity_tier == SpecificityTier.SUPERFAMILY
    assert res.consensus_name == "NAD(P)-binding domain superfamily"
    fam_cands = [c for c in res.top_candidates if c.specificity_tier == SpecificityTier.PROTEIN_FAMILY_CONTEXT]
    assert fam_cands
    assert all("[" not in c.canonical_name for c in fam_cands)

    # Strict 1:1 Pfam family (IoU >= 0.65, annot_cov >= 0.70, ted_cov >= 0.65) still qualifies as Tier 1
    m_pf_fam = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-P11111-F1-model_v4_TED01",
            uniprot_acc="P11111",
            chopping="10-150",
            nres_domain=141,
            cath_label="1.20.58.10",
        ),
        overlapping_annotations=[
            OverlappingAnnotation(
                source_db="pfam",
                accession="PF12345",
                name="Small archaeal fold-X",
                annotation_type="family",
                chopping="8-152",
                iou=0.92,
                ted_coverage=0.98,
                annot_coverage=0.94,
                overlap_category=OverlapCategory.EXACT_MATCH,
            )
        ],
        relative_position="full-length",
    )
    res_pf = aggregate_cluster_evidence("pfam_1to1_family", [m_pf_fam])
    assert res_pf.specificity_tier == SpecificityTier.SPECIFIC_DOMAIN
    assert "small archaeal fold-x" in res_pf.consensus_name.lower()


def test_naming_hygiene_sentences_bare_positions_and_quality_flags() -> None:
    """Verify PROSITE conserved_site demotion, sentence penalization, bare position enrichment, and quality flags."""
    m = EnrichedTedDomain(
        domain=TedDomainMember(
            ted_id="AF-Q8U3D2-F1-model_v4_TED01",
            uniprot_acc="Q8U3D2",
            chopping="1-146_284-319",
            nres_domain=182,
            plddt=55.2,
            norm_rg=0.58,
            consensus_level="low",
            num_helix=8,
            num_strand=9,
            cath_label="-",
        ),
        protein_context=ProteinContext(
            uniprot_acc="Q8U3D2",
            recommended_name="Argonaute protein",
            is_reviewed=True,
            protein_length=770,
            similarity_families=["Argonaute family"],
        ),
        overlapping_annotations=[
            OverlappingAnnotation(
                source_db="prosite",
                accession="PS00107",
                name="Protein kinases active-site signature",
                annotation_type="conserved_site",
                chopping="100-115",
                iou=0.09,
                ted_coverage=0.09,
                annot_coverage=1.0,
                overlap_category=OverlapCategory.SUB_DOMAIN,
            ),
            OverlappingAnnotation(
                source_db="uniprot_feature",
                accession="UP_FEAT:Q8U3D2:10-145",
                name="Domain I interacts with DnaA modulators",
                annotation_type="region",
                chopping="10-145",
                iou=0.74,
                ted_coverage=0.74,
                annot_coverage=1.0,
                overlap_category=OverlapCategory.EXACT_MATCH,
            ),
            OverlappingAnnotation(
                source_db="uniprot_feature",
                accession="UP_FEAT:Q8U3D2:1-145_285-318",
                name="N-terminal domain",
                annotation_type="domain",
                chopping="1-145_285-318",
                iou=0.98,
                ted_coverage=0.98,
                annot_coverage=1.0,
                overlap_category=OverlapCategory.EXACT_MATCH,
            ),
        ],
        relative_position="N-terminal",
    )

    res = aggregate_cluster_evidence("hygiene_test", [m])
    assert res.specificity_tier == SpecificityTier.SPECIFIC_DOMAIN
    # Bare 'N-terminal domain' enriched with parent family -> 'Argonaute N-terminal domain'
    assert res.consensus_name == "Argonaute N-terminal domain"
    # Structure quality flags tracked without artificially ruining annotation confidence
    assert set(res.structure_quality_flags) == {"low_plddt", "low_consensus", "non_globular"}
    assert res.confidence_score >= 0.90
    assert res.novel_fold_candidate is True
    assert res.dominant_cath_label == "3.?.?.?"
    assert res.cath_fold_hint is not None

