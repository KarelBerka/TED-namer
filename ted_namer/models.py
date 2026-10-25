"""Data models for TED domains, residue boundaries, annotations, and cluster naming results."""

from __future__ import annotations

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field


class OverlapCategory(str, Enum):
    """Spatial relationship between a sequence annotation and a TED structural domain."""

    EXACT_MATCH = "exact_match"
    SUB_DOMAIN = "sub_domain"
    SUPRA_DOMAIN = "supra_domain"
    PARTIAL_OVERLAP = "partial_overlap"
    DISJOINT = "disjoint"


class SpecificityTier(str, Enum):
    """Specificity hierarchy for domain naming candidates."""

    SPECIFIC_DOMAIN = "specific_domain"  # Tier 1: CDD/Pfam/SMART/PROSITE/UniProt Domain
    SUPERFAMILY = "superfamily"  # Tier 2: CATH H-level / Gene3D / SUPERFAMILY / InterPro Homologous Superfamily
    PROTEIN_FAMILY_CONTEXT = "protein_family_context"  # Tier 3: InterPro Family / PANTHER / UniProt family + position
    TOPOLOGY_FOLD = "topology_fold"  # Tier 4: CATH T-level fold / topology
    NOVEL_UNCHARACTERIZED = "novel_uncharacterized"  # Fallback when only uncharacterized/DUF evidence exists


class SegmentRange(BaseModel):
    """A 1-based inclusive contiguous residue segment [start, end]."""

    start: int
    end: int

    @property
    def length(self) -> int:
        return max(0, self.end - self.start + 1)

    def to_set(self) -> set[int]:
        if self.end < self.start:
            return set()
        return set(range(self.start, self.end + 1))

    def __str__(self) -> str:
        return f"{self.start}-{self.end}"


class TedDomainMember(BaseModel):
    """A single TED structural domain instance on a UniProt protein."""

    ted_id: str
    uniprot_acc: str
    chopping: str
    nres_domain: int = 0
    num_segments: int = 1
    plddt: float = 0.0
    num_helix_strand_turn: int = 0
    num_helix: int = 0
    num_strand: int = 0
    num_helix_strand: int = 0
    num_turn: int = 0
    cath_label: str = "-"
    cath_assignment_level: str = "-"
    cath_assignment_method: str = "-"
    inferred_cath_label: Optional[str] = Field(
        default=None,
        description="Full 4-level CATH label where estimated levels are prefixed with 'T' (e.g. '3.T40.T50.T10')",
    )
    cath_estimation_source: Optional[str] = Field(
        default=None,
        description="How the CATH label was obtained: 'ted_official', 'g3dsa_overlap', 'sse_heuristic'",
    )
    packing_density: Optional[float] = None
    norm_rg: Optional[float] = None
    tax_common_name: Optional[str] = None
    tax_scientific_name: Optional[str] = None
    tax_lineage: Optional[str] = None

    @property
    def short_ted_num(self) -> str:
        """Extracts e.g. 'TED01' from 'AF-P04637-F1-model_v4_TED01'."""
        if "_TED" in self.ted_id:
            return "TED" + self.ted_id.rsplit("_TED", 1)[1]
        return self.ted_id


class OverlappingAnnotation(BaseModel):
    """A sequence or structural annotation overlapping a specific TED domain."""

    source_db: str = Field(description="e.g. cdd, pfam, smart, prosite, panther, cathgene3d, ssf, uniprot_feature")
    accession: str = Field(description="Database accession e.g. PF00870, cd08367, G3DSA:2.60.40.720")
    name: str = Field(description="Human-readable annotation name")
    short_name: Optional[str] = None
    integrated_ipr: Optional[str] = Field(default=None, description="Parent InterPro accession e.g. IPR011615")
    integrated_ipr_name: Optional[str] = None
    annotation_type: str = Field(
        default="domain",
        description="domain, repeat, homologous_superfamily, family, region, site, etc.",
    )
    chopping: str = Field(description="Residue boundaries of this annotation hit, e.g. '100-290'")
    iou: float = Field(description="Jaccard intersection-over-union with the TED domain residues")
    ted_coverage: float = Field(description="Fraction of TED domain residues covered by this annotation")
    annot_coverage: float = Field(description="Fraction of annotation residues covered by the TED domain")
    overlap_category: OverlapCategory = OverlapCategory.PARTIAL_OVERLAP


class CathResolution(BaseModel):
    """Resolved human-readable names for a CATH code (H-level and/or T-level), supporting 'T'-prefixed estimates."""

    cath_label: str = Field(description="Formatted CATH code (e.g. '2.60.40.720' or '3.T40.T50.T10')")
    numeric_cath_code: Optional[str] = Field(
        default=None,
        description="Underlying numeric CATH code without 'T' prefixes (e.g. '3.40.50.10')",
    )
    raw_ted_cath_label: str = "-"
    is_estimated: bool = False
    estimation_source: str = "ted_official"
    assignment_level: str = "-"  # "H", "T", or "inferred"
    superfamily_name: Optional[str] = None  # H-level name from Gene3D / InterPro / CATH
    integrated_ipr: Optional[str] = None
    topology_code: Optional[str] = None  # 3-part C.A.T code (with T prefixes if estimated)
    topology_name: Optional[str] = None  # T-level fold name from CATH cathtree
    architecture_name: Optional[str] = None  # A-level name
    class_name: Optional[str] = None  # C-level name
    funfam_names: list[str] = Field(default_factory=list)
    representative_pdb_title: Optional[str] = None


class ProteinContext(BaseModel):
    """Protein-level metadata from UniProtKB for a TED domain's parent protein."""

    uniprot_acc: str
    entry_name: Optional[str] = None
    recommended_name: Optional[str] = None
    gene_name: Optional[str] = None
    is_reviewed: bool = False
    protein_length: int = 0
    organism: Optional[str] = None
    similarity_families: list[str] = Field(default_factory=list, description="Extracted from SIMILARITY comments")
    go_molecular_functions: list[str] = Field(default_factory=list)
    go_biological_processes: list[str] = Field(default_factory=list)


class EnrichedTedDomain(BaseModel):
    """Full evidence collected for a single TED domain member."""

    domain: TedDomainMember
    cath_info: Optional[CathResolution] = None
    protein_context: Optional[ProteinContext] = None
    overlapping_annotations: list[OverlappingAnnotation] = Field(default_factory=list)
    relative_position: str = Field(
        default="full-length",
        description="N-terminal, central, C-terminal, or full-length relative to parent protein",
    )


class ClusterCandidateTerm(BaseModel):
    """An aggregated candidate name scored across a cluster of TED domains."""

    canonical_name: str
    short_name: Optional[str] = None
    specificity_tier: SpecificityTier
    source_databases: list[str] = Field(default_factory=list)
    supporting_accessions: list[str] = Field(default_factory=list)
    member_hits: int = Field(description="Number of cluster members carrying this annotation")
    cluster_frequency: float = Field(description="Fraction of cluster members [0.0 - 1.0]")
    mean_iou: float = 0.0
    mean_ted_coverage: float = 0.0
    mean_annot_coverage: float = 0.0
    reviewed_hit_count: int = 0
    is_low_information: bool = False
    composite_score: float = 0.0


class ClusterNamingResult(BaseModel):
    """Final naming output for a TED domain cluster or family."""

    cluster_id: str
    member_count: int
    consensus_name: str
    short_name: str
    confidence_score: float = Field(ge=0.0, le=1.0)
    specificity_tier: SpecificityTier
    naming_method: str = Field(description="'rule_based' or 'llm_consensus'")
    rationale: str
    dominant_cath_label: str = "-"
    cath_summary: Optional[CathResolution] = None
    consensus_relative_position: Optional[str] = None
    taxonomic_summary: list[str] = Field(default_factory=list)
    top_candidates: list[ClusterCandidateTerm] = Field(default_factory=list)
    members: list[EnrichedTedDomain] = Field(default_factory=list)
    extra_metadata: dict[str, Any] = Field(default_factory=dict)
