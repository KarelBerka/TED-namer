"""Cluster-level evidence aggregation, synonym unification, naming hygiene, and 4-tier candidate scoring."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional

from ted_namer.cath_classifier import (
    CURATED_UNINTEGRATED_CATH_H_NAMES,
    _G3DSA_RE,
    infer_cath_classification,
    is_estimated_cath_code,
    is_novel_fold_source,
    strip_t_prefixes,
)
from ted_namer.models import (
    CathResolution,
    ClusterCandidateTerm,
    ClusterNamingResult,
    EnrichedTedDomain,
    OverlapCategory,
    OverlappingAnnotation,
    SpecificityTier,
)
from ted_namer.overlap import parse_chopping

_LOW_INFO_PATTERNS = re.compile(
    r"(\b(?:uncharacterized|uncharacterised|hypothetical|putative|duf\d+|domain of unknown function|"
    r"unnamed protein|predicted protein|whole genome shotgun|orf|upf\d+|unknown\s+domain|unknown\s+protein)\b"
    r"|^unknown$)",
    re.IGNORECASE,
)

_LEGACY_PDB_HEADER_RE = re.compile(
    r"(;\s*domain\s+\d+|\bdomain\s+[1-9]\b|\bsubunit\s+[A-Z0-9]\b|\bchain\s+[A-Z0-9]\b|"
    r"Transferase\(|Phosphorylase Kinase|Vaccinia Virus)",
    re.IGNORECASE,
)

_SENTENCE_LIKE_RE = re.compile(
    r"\b(is involved in|involved in|interacts with|interacts|mediates|comprises|"
    r"required for|essential for|responsible for|binds to|important for|found in .+ and related)\b",
    re.IGNORECASE,
)

_BARE_POSITIONAL_RE = re.compile(
    r"^(N-terminal|C-terminal|N-term|C-term|central|middle|NTD|CTD|N|C|linker|tail|core)"
    r"(\s+domain|\s+subdomain|\s+region)?$",
    re.IGNORECASE,
)

# Strictly domain/repeat curated databases (excluding whole-protein HAMAP/NCBIfam families and short PROSITE sites)
_SPECIFIC_DOMAIN_DBS = {"cdd", "pfam", "smart", "profile", "uniprot_feature"}
_SUBLOBE_DOMAIN_DBS = {"cdd", "pfam", "smart", "profile", "interpro", "ncbifam", "uniprot_feature"}
_SUPERFAMILY_DBS = {"cathgene3d", "ssf", "supfam"}

_TIER_WEIGHTS: dict[SpecificityTier, float] = {
    SpecificityTier.SPECIFIC_DOMAIN: 1.00,
    SpecificityTier.SUPERFAMILY: 0.85,
    SpecificityTier.PROTEIN_FAMILY_CONTEXT: 0.65,
    SpecificityTier.TOPOLOGY_FOLD: 0.45,
    SpecificityTier.NOVEL_UNCHARACTERIZED: 0.25,
}


def is_low_information_name(name: str) -> bool:
    """Return True if a candidate name is generic/uncharacterized (e.g. DUF, Uncharacterized protein, Unknown domain)."""
    if not name or len(name.strip()) < 3:
        return True
    return bool(_LOW_INFO_PATTERNS.search(name.strip()))


def is_sentence_like_description(name: str) -> bool:
    """Return True if a string reads like functional prose rather than a concise domain name."""
    txt = name.strip()
    if len(txt) > 90:
        return True
    return bool(_SENTENCE_LIKE_RE.search(txt))


def is_bare_positional_name(name: str) -> bool:
    """Return True if a name gives only position without protein/family context (e.g. 'N-terminal domain')."""
    return bool(_BARE_POSITIONAL_RE.match(name.strip()))


def strip_bracketed_tags(name: str) -> str:
    """Remove trailing or inline bracketed gene symbols and EC numbers such as '[bioU]'."""
    return re.sub(r"\s*\[[^\]]+\]", "", name).strip()


def normalize_domain_name(name: str) -> str:
    """Normalize a domain/superfamily name for deduplication across databases."""
    cleaned = strip_bracketed_tags(name)
    cleaned = re.sub(r"[,_/()-]+", " ", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned).strip().lower()
    # Strip trailing generic words for key comparison
    cleaned = re.sub(r"\s+(homologous superfamily|superfamily|domain|family|fold|like)$", "", cleaned)
    return cleaned.strip()


def clean_display_domain_name(name: str, tier: SpecificityTier) -> str:
    """Standardize capitalization, strip bracketed gene tags and redundant suffixes."""
    txt = strip_bracketed_tags(name)
    txt = re.sub(r"\s+", " ", txt.strip().rstrip(".").replace(", ", " "))
    if not txt:
        return "Unassigned TED structural domain"

    # Strip Chain suffixes from CATH/Gene3D representative titles while leaving numbered legacy PDB headers detectable
    txt = re.sub(r";\s*Chain\s+[A-Za-z0-9]+.*$", "", txt, flags=re.IGNORECASE)
    txt = re.sub(r"\.\s*Chain\s+[A-Za-z0-9]+.*$", "", txt, flags=re.IGNORECASE)
    txt = txt.strip().rstrip(";.,")

    # Preserve capitalization unless the entire string is ALL CAPS (like PANTHER 'CELLULAR TUMOR ANTIGEN P53')
    if txt.isupper() and len(txt) > 6:
        txt = txt.capitalize()

    # Deduplicate redundant 'domain' words
    txt = re.sub(r"\bdomain-like\s+domain\b", "domain-like", txt, flags=re.IGNORECASE)
    txt = re.sub(r"\bdomain\s+domain\b", "domain", txt, flags=re.IGNORECASE)

    low = txt.lower()
    if tier == SpecificityTier.SPECIFIC_DOMAIN:
        # Remove trailing 'superfamily' if present on a specific domain
        txt = re.sub(r"\s+superfamily$", "", txt, flags=re.IGNORECASE)
        has_structural_noun = any(
            w in txt.lower()
            for w in (
                "domain",
                "subdomain",
                "repeat",
                "finger",
                "motif",
                "region",
                "lobe",
                "barrel",
                "helix",
                "bundle",
                "fold",
            )
        )
        if not has_structural_noun:
            txt = f"{txt} domain"
        elif txt.lower().endswith("motif") and "domain" not in txt.lower():
            txt = f"{txt} domain"
    elif tier == SpecificityTier.TOPOLOGY_FOLD:
        if "fold" not in low and "topology" not in low and "domain" not in low:
            txt = f"{txt} fold"
    return re.sub(r"\s+", " ", txt).strip()


def make_short_label(name: str, max_len: int = 24) -> str:
    """Generate a compact short label (<= max_len chars) suitable for tables and UI badges."""
    replacements = [
        (r"\btranscription factor\b", "TF"),
        (r"\bDNA-binding domain\b", "DNA-bd"),
        (r"\bRNA-binding domain\b", "RNA-bd"),
        (r"\btetramerisation domain\b", "tetramer"),
        (r"\btetramerization domain\b", "tetramer"),
        (r"\boligomerisation domain\b", "oligomer"),
        (r"\boligomerization domain\b", "oligomer"),
        (r"\btransactivation domain\b", "TAD"),
        (r"\bimmunoglobulin-like\b", "Ig-like"),
        (r"\bimmunoglobulin\b", "Ig"),
        (r"\bhomologous superfamily\b", "sf"),
        (r"\bsuperfamily\b", "sf"),
        (r"\bN-terminal\b", "N-term"),
        (r"\bC-terminal\b", "C-term"),
        (r"\bcentral\b", "central"),
        (r"\bdomain\b", "dom"),
        (r"\bfamily\b", "fam"),
    ]
    short = strip_bracketed_tags(name)
    for pat, rep in replacements:
        short = re.sub(pat, rep, short, flags=re.IGNORECASE)
    short = re.sub(r"[,()]+", "", short)
    short = re.sub(r"\s+", "_", short.strip())
    if len(short) <= max_len:
        return short
    return short[:max_len].rstrip("_-")


def _is_curated_subdomain_lobe(annot: OverlappingAnnotation) -> bool:
    """Return True if a TED domain is a structural lobe/subdomain inside a curated functional domain."""
    atype = annot.annotation_type.lower()
    db = annot.source_db.lower()
    return (
        atype == "domain"
        and db in _SUBLOBE_DOMAIN_DBS
        and annot.ted_coverage >= 0.80
        and 0.20 <= annot.annot_coverage < 0.50
        and annot.iou < 0.35
    )


def _compute_subdomain_position(ted_chopping: str, annot_chopping: str) -> str:
    """Compute the relative position of a TED structural lobe within an enclosing sequence domain."""
    t_segs = parse_chopping(ted_chopping)
    a_segs = parse_chopping(annot_chopping)
    if not t_segs or not a_segs:
        return "subdomain"
    a_start = a_segs[0].start
    a_end = a_segs[-1].end
    a_span = max(1, a_end - a_start + 1)
    t_mid = 0.5 * (t_segs[0].start + t_segs[-1].end)
    rel = (t_mid - a_start) / a_span
    if rel < 0.38:
        return "N-terminal subdomain"
    if rel > 0.62:
        return "C-terminal subdomain"
    return "central subdomain"


def determine_annotation_tier(annot: OverlappingAnnotation) -> SpecificityTier:
    """Map an OverlappingAnnotation to its SpecificityTier strictly honoring InterPro annotation_type."""
    if is_low_information_name(annot.name):
        return SpecificityTier.NOVEL_UNCHARACTERIZED

    atype = annot.annotation_type.lower()
    db = annot.source_db.lower()

    # Conserved sites, disordered regions, coiled coils, and PTMs are never full specific domains
    if atype in {"conserved_site", "active_site", "binding_site", "ptm", "disordered", "coiled_coil"}:
        return SpecificityTier.PROTEIN_FAMILY_CONTEXT

    # Direct domain / repeat matches (including NCBIfam entries explicitly typed 'domain'/'repeat')
    if atype in {"domain", "repeat"}:
        if annot.iou >= 0.35 or (annot.ted_coverage >= 0.50 and annot.annot_coverage >= 0.50):
            return SpecificityTier.SPECIFIC_DOMAIN
        if _is_curated_subdomain_lobe(annot):
            return SpecificityTier.SPECIFIC_DOMAIN
        return SpecificityTier.PROTEIN_FAMILY_CONTEXT

    # Family annotations: allow only single-domain Pfam families with strict 1:1 structural correspondence
    if atype == "family":
        if (
            db == "pfam"
            and annot.iou >= 0.65
            and annot.annot_coverage >= 0.70
            and annot.ted_coverage >= 0.65
        ):
            return SpecificityTier.SPECIFIC_DOMAIN
        return SpecificityTier.PROTEIN_FAMILY_CONTEXT

    # UniProt regions/motifs: promote to Tier 1 only on strict 1:1 overlap with a concise non-sentence label
    if atype in {"motif", "region"}:
        if annot.iou >= 0.65 and annot.ted_coverage >= 0.70 and not is_sentence_like_description(annot.name):
            return SpecificityTier.SPECIFIC_DOMAIN
        return SpecificityTier.PROTEIN_FAMILY_CONTEXT

    if atype == "homologous_superfamily" or db in _SUPERFAMILY_DBS:
        if annot.iou >= 0.30 or annot.ted_coverage >= 0.50:
            return SpecificityTier.SUPERFAMILY

    return SpecificityTier.PROTEIN_FAMILY_CONTEXT


def _extract_parent_context_prefix(member: EnrichedTedDomain) -> Optional[str]:
    """Extract an informative protein/family prefix for enriching bare positional domain names."""
    ctx = member.protein_context
    if not ctx:
        return None
    for fam in ctx.similarity_families:
        base_fam = re.sub(r"\s+(subfamily|family|superfamily).*$", "", fam, flags=re.IGNORECASE).strip()
        base_fam = strip_bracketed_tags(base_fam)
        if base_fam and not is_low_information_name(base_fam) and not is_bare_positional_name(base_fam):
            return base_fam
    if ctx.recommended_name:
        rec = strip_bracketed_tags(ctx.recommended_name).strip()
        rec = re.sub(r"\s+(protein|homolog)$", "", rec, flags=re.IGNORECASE).strip()
        if rec and not is_low_information_name(rec) and not is_bare_positional_name(rec) and len(rec) <= 65:
            return rec
    return None


def _enrich_bare_positional_label(raw_disp: str, parent_ctx: Optional[str]) -> str:
    """Turn bare positional labels like 'N-terminal domain' or 'NTD domain' into '<Parent> N-terminal domain'."""
    m = _BARE_POSITIONAL_RE.match(raw_disp.strip())
    if not m:
        return raw_disp
    pos_token = m.group(1).strip()
    suffix_token = (m.group(2) or "domain").strip().lower()
    pos_map = {
        "ntd": "N-terminal",
        "n-term": "N-terminal",
        "n": "N-terminal",
        "ctd": "C-terminal",
        "c-term": "C-terminal",
        "c": "C-terminal",
        "middle": "central",
    }
    norm_pos = pos_map.get(pos_token.lower(), pos_token)
    if norm_pos.lower() == "n-terminal":
        norm_pos = "N-terminal"
    elif norm_pos.lower() == "c-terminal":
        norm_pos = "C-terminal"
    elif norm_pos.lower() == "central":
        norm_pos = "central"

    if parent_ctx:
        return f"{parent_ctx} {norm_pos} {suffix_token}"
    return raw_disp


def _extract_cath_h_code(text: Optional[str]) -> Optional[str]:
    if not text:
        return None
    m = _G3DSA_RE.match(text.strip())
    if m:
        return f"{m.group(1)}.{m.group(2)}.{m.group(3)}.{m.group(4)}"
    return None


@dataclass
class _CandidateAccumulator:
    """Internal accumulator for grouping synonymous annotations across cluster members."""

    group_key: str
    names_counter: Counter[str] = field(default_factory=Counter)
    short_names: list[str] = field(default_factory=list)
    tiers: list[SpecificityTier] = field(default_factory=list)
    source_dbs: set[str] = field(default_factory=set)
    accessions: set[str] = field(default_factory=set)
    member_Best_iou: dict[str, float] = field(default_factory=dict)
    member_best_ted_cov: dict[str, float] = field(default_factory=dict)
    member_best_annot_cov: dict[str, float] = field(default_factory=dict)
    reviewed_members: set[str] = field(default_factory=set)


def _compute_structure_quality_flags(members: list[EnrichedTedDomain]) -> list[str]:
    """Compute structural quality warning flags independently from annotation confidence."""
    flags: list[str] = []
    plddts = [m.domain.plddt for m in members if m.domain.plddt and m.domain.plddt > 0]
    if plddts and ((sum(plddts) / len(plddts)) < 70.0 or any(p < 60.0 for p in plddts)):
        flags.append("low_plddt")
    if any((m.domain.consensus_level or "").strip().lower() == "low" for m in members):
        flags.append("low_consensus")
    if any(m.domain.norm_rg is not None and m.domain.norm_rg > 0.50 for m in members):
        flags.append("non_globular")
    return flags if flags else ["ok"]


def aggregate_cluster_evidence(
    cluster_id: str,
    members: list[EnrichedTedDomain],
) -> ClusterNamingResult:
    """Aggregate all annotations across a cluster of TED domains and compute ranked candidate names."""
    total_members = len(members)
    if total_members == 0:
        return ClusterNamingResult(
            cluster_id=cluster_id,
            member_count=0,
            consensus_name="Empty TED cluster",
            short_name="empty_cluster",
            confidence_score=0.0,
            specificity_tier=SpecificityTier.NOVEL_UNCHARACTERIZED,
            naming_method="rule_based",
            rationale="No TED domain members were provided for this cluster.",
            structure_quality_flags=["ok"],
        )

    # Ensure every member has inferred_cath_label, cath_fold_hint, and novel_fold_candidate populated
    for m in members:
        if not m.domain.inferred_cath_label:
            inferred_code, numeric_hint, est_src = infer_cath_classification(
                member=m.domain,
                overlapping_annotations=m.overlapping_annotations,
            )
            m.domain.inferred_cath_label = inferred_code
            m.domain.cath_estimation_source = est_src
            m.domain.cath_fold_hint = numeric_hint if is_estimated_cath_code(inferred_code) else None
            m.domain.novel_fold_candidate = is_novel_fold_source(est_src)
        else:
            m.domain.novel_fold_candidate = is_novel_fold_source(m.domain.cath_estimation_source)

    # Determine dominant CATH label and consensus relative position
    effective_cath_counts = Counter(
        (m.domain.inferred_cath_label or m.domain.cath_label)
        for m in members
        if (m.domain.inferred_cath_label or m.domain.cath_label)
        and (m.domain.inferred_cath_label or m.domain.cath_label) != "-"
    )
    dominant_cath = effective_cath_counts.most_common(1)[0][0] if effective_cath_counts else "-"

    cath_summary: Optional[CathResolution] = None
    for m in members:
        eff_code = m.domain.inferred_cath_label or m.domain.cath_label
        if m.cath_info and (m.cath_info.cath_label == dominant_cath or eff_code == dominant_cath):
            cath_summary = m.cath_info
            break

    hint_counts = Counter(
        m.domain.cath_fold_hint or (m.cath_info.fold_hint if m.cath_info else None)
        for m in members
        if (m.domain.cath_fold_hint or (m.cath_info and m.cath_info.fold_hint))
    )
    dominant_fold_hint = (
        hint_counts.most_common(1)[0][0]
        if (hint_counts and is_estimated_cath_code(dominant_cath))
        else None
    )
    novel_fold_candidate = all(m.domain.novel_fold_candidate for m in members)
    structure_quality_flags = _compute_structure_quality_flags(members)

    pos_counts = Counter(m.relative_position for m in members if m.relative_position)
    consensus_pos = pos_counts.most_common(1)[0][0] if pos_counts else "domain"

    tax_counts = Counter(
        m.domain.tax_scientific_name for m in members if m.domain.tax_scientific_name
    )
    taxonomic_summary = [f"{tax} ({cnt})" for tax, cnt in tax_counts.most_common(5)]

    accumulators: dict[str, _CandidateAccumulator] = {}

    def _add_hit(
        group_key: str,
        display_name: str,
        short_name: Optional[str],
        tier: SpecificityTier,
        source_db: str,
        accession: str,
        member_id: str,
        iou: float,
        ted_cov: float,
        annot_cov: float,
        is_reviewed: bool,
    ) -> None:
        acc_obj = accumulators.get(group_key)
        if acc_obj is None:
            acc_obj = _CandidateAccumulator(group_key=group_key)
            accumulators[group_key] = acc_obj

        # Prefer curated InterPro names (3) and specific domain DBs (2) over raw superfamily/PDB titles (1)
        if source_db == "interpro":
            weight = 3
        elif source_db in {"cdd", "pfam", "smart", "prosite", "profile", "uniprot_feature", "ncbifam"}:
            weight = 2
        else:
            weight = 1
        acc_obj.names_counter[display_name] += weight
        if short_name and short_name not in acc_obj.short_names:
            acc_obj.short_names.append(short_name)
        acc_obj.tiers.append(tier)
        acc_obj.source_dbs.add(source_db)
        if accession:
            acc_obj.accessions.add(accession)

        prev_iou = acc_obj.member_Best_iou.get(member_id, -1.0)
        if iou > prev_iou:
            acc_obj.member_Best_iou[member_id] = iou
            acc_obj.member_best_ted_cov[member_id] = ted_cov
            acc_obj.member_best_annot_cov[member_id] = annot_cov
        if is_reviewed:
            acc_obj.reviewed_members.add(member_id)

    for m in members:
        member_id = m.domain.ted_id
        is_rev = bool(m.protein_context and m.protein_context.is_reviewed)
        pos_qualifier = (
            f"{m.relative_position} domain"
            if m.relative_position not in {"full-length", "domain"}
            else "domain"
        )
        parent_ctx = _extract_parent_context_prefix(m)

        # 1. Overlapping sequence/structure annotations
        for annot in m.overlapping_annotations:
            tier = determine_annotation_tier(annot)
            ipr_group_id = annot.integrated_ipr or (
                annot.accession if annot.source_db == "interpro" and annot.accession.startswith("IPR") else None
            )

            disp = strip_bracketed_tags(annot.name)

            # Replace raw G3DSA/SSF legacy PDB titles with curated integrated IPR name or curated unintegrated H-name
            g3d_h = _extract_cath_h_code(annot.accession) if annot.source_db == "cathgene3d" else None
            if annot.source_db in {"cathgene3d", "ssf"}:
                if annot.integrated_ipr_name:
                    disp = strip_bracketed_tags(annot.integrated_ipr_name)
                elif g3d_h and g3d_h in CURATED_UNINTEGRATED_CATH_H_NAMES:
                    disp = CURATED_UNINTEGRATED_CATH_H_NAMES[g3d_h]
            elif annot.integrated_ipr_name and _LEGACY_PDB_HEADER_RE.search(disp):
                disp = strip_bracketed_tags(annot.integrated_ipr_name)

            # Sentence-like description sanitization: prefer concise integrated IPR or informative short_name
            if is_sentence_like_description(disp):
                if annot.integrated_ipr_name and not is_sentence_like_description(annot.integrated_ipr_name):
                    disp = strip_bracketed_tags(annot.integrated_ipr_name)
                elif (
                    annot.short_name
                    and len(annot.short_name.strip()) >= 4
                    and not is_sentence_like_description(annot.short_name)
                    and not is_low_information_name(annot.short_name)
                ):
                    disp = strip_bracketed_tags(annot.short_name)

            # Enrich bare positional labels (e.g. 'N-terminal domain', 'NTD domain') with parent protein context
            was_bare_pos = is_bare_positional_name(disp)
            is_lobe = _is_curated_subdomain_lobe(annot)
            if was_bare_pos and not is_lobe:
                disp = _enrich_bare_positional_label(disp, parent_ctx)

            eff_iou = annot.iou
            if is_lobe:
                sub_pos = _compute_subdomain_position(m.domain.chopping, annot.chopping)
                base_dom = (
                    (parent_ctx or "")
                    if was_bare_pos
                    else re.sub(r"\s+domain$", "", disp, flags=re.IGNORECASE).strip()
                )
                if base_dom:
                    disp = f"{base_dom} {sub_pos}"
                    gkey = f"SUBLOBE:{ipr_group_id or normalize_domain_name(base_dom)}:{sub_pos}"
                else:
                    disp = sub_pos
                    gkey = f"SUBLOBE:{ipr_group_id or 'bare'}:{sub_pos}"
                eff_iou = round(0.55 * annot.ted_coverage + 0.45 * annot.annot_coverage, 4)
            else:
                if ipr_group_id and tier in {SpecificityTier.SPECIFIC_DOMAIN, SpecificityTier.SUPERFAMILY}:
                    gkey = f"IPR:{ipr_group_id}:{tier.value}"
                elif g3d_h and not annot.integrated_ipr and tier == SpecificityTier.SUPERFAMILY:
                    gkey = f"CATH_H:{g3d_h}:{tier.value}"
                elif was_bare_pos and parent_ctx:
                    gkey = f"CTX_POS:{normalize_domain_name(disp)}:{tier.value}"
                else:
                    gkey = f"NORM:{normalize_domain_name(disp)}:{tier.value}"

                if tier == SpecificityTier.PROTEIN_FAMILY_CONTEXT and "domain" not in disp.lower():
                    # Qualify supra-domain family hit by the TED domain's relative position
                    base_fam = re.sub(r"\s+family$", "", disp, flags=re.IGNORECASE).strip()
                    disp = f"{base_fam} {pos_qualifier}"

            _add_hit(
                group_key=gkey,
                display_name=disp,
                short_name=annot.short_name,
                tier=tier,
                source_db=annot.source_db,
                accession=annot.accession,
                member_id=member_id,
                iou=eff_iou,
                ted_cov=annot.ted_coverage,
                annot_cov=annot.annot_coverage,
                is_reviewed=is_rev,
            )

        # 2. CATH resolution evidence (H-level superfamily and confirmed T-level fold)
        if m.cath_info:
            if m.cath_info.superfamily_name:
                h_code_num = (
                    _extract_cath_h_code(m.cath_info.cath_label)
                    or _extract_cath_h_code(m.cath_info.fold_hint)
                    or strip_t_prefixes(m.cath_info.cath_label)
                )
                sf_disp = strip_bracketed_tags(m.cath_info.superfamily_name)
                if (
                    not m.cath_info.integrated_ipr
                    and h_code_num in CURATED_UNINTEGRATED_CATH_H_NAMES
                ):
                    sf_disp = CURATED_UNINTEGRATED_CATH_H_NAMES[h_code_num]

                # Use real G3DSA overlap metrics when present instead of inflated hardcoded 0.85/0.90
                matched_g3d = None
                for a in m.overlapping_annotations:
                    if a.source_db == "cathgene3d" and h_code_num and a.accession.upper() == f"G3DSA:{h_code_num}":
                        matched_g3d = a
                        break

                h_iou = matched_g3d.iou if matched_g3d else 0.65
                h_tcov = matched_g3d.ted_coverage if matched_g3d else 0.75
                h_acov = matched_g3d.annot_coverage if matched_g3d else 0.65

                if m.cath_info.integrated_ipr:
                    gkey = f"IPR:{m.cath_info.integrated_ipr}:{SpecificityTier.SUPERFAMILY.value}"
                else:
                    gkey = f"CATH_H:{h_code_num}:{SpecificityTier.SUPERFAMILY.value}"

                _add_hit(
                    group_key=gkey,
                    display_name=sf_disp,
                    short_name=None,
                    tier=SpecificityTier.SUPERFAMILY,
                    source_db="cath_h",
                    accession=f"CATH:{h_code_num}",
                    member_id=member_id,
                    iou=h_iou,
                    ted_cov=h_tcov,
                    annot_cov=h_acov,
                    is_reviewed=is_rev,
                )

            # Only emit Topology fallback when T-level is actually confirmed (never for '?' topology)
            if (
                m.cath_info.topology_name
                and m.cath_info.topology_code
                and "?" not in m.cath_info.topology_code
                and not is_estimated_cath_code(m.cath_info.topology_code)
            ):
                gkey = f"CATH_T:{m.cath_info.topology_code}"
                topo_disp = m.cath_info.topology_name.strip()
                if not any(w in topo_disp.lower() for w in ("fold", "topology", "domain")):
                    topo_disp = f"{topo_disp} fold"
                _add_hit(
                    group_key=gkey,
                    display_name=topo_disp,
                    short_name=None,
                    tier=SpecificityTier.TOPOLOGY_FOLD,
                    source_db="cath_t",
                    accession=f"CATH:{m.cath_info.topology_code}",
                    member_id=member_id,
                    iou=0.45,
                    ted_cov=0.60,
                    annot_cov=0.45,
                    is_reviewed=is_rev,
                )

        # 3. Parent UniProt protein family & recommended name context (using true domain/protein length ratio)
        if m.protein_context:
            prot_len = (
                m.protein_context.protein_length
                if (m.protein_context.protein_length and m.protein_context.protein_length > 0)
                else max(m.domain.nres_domain, 250)
            )
            dom_ratio = round(min(1.0, max(0.05, m.domain.nres_domain / max(1, prot_len))), 4)
            ctx_iou = round(0.30 + 0.60 * dom_ratio, 4)

            for fam in m.protein_context.similarity_families:
                base_fam = re.sub(r"\s+(subfamily|family|superfamily).*$", "", fam, flags=re.IGNORECASE).strip()
                base_fam = strip_bracketed_tags(base_fam)
                if base_fam:
                    fam_dom_name = f"{base_fam} {pos_qualifier}"
                    tier = (
                        SpecificityTier.NOVEL_UNCHARACTERIZED
                        if is_low_information_name(base_fam)
                        else SpecificityTier.PROTEIN_FAMILY_CONTEXT
                    )
                    _add_hit(
                        group_key=f"UP_FAM:{normalize_domain_name(base_fam)}:{pos_qualifier}",
                        display_name=fam_dom_name,
                        short_name=None,
                        tier=tier,
                        source_db="uniprot_family",
                        accession=m.protein_context.uniprot_acc,
                        member_id=member_id,
                        iou=ctx_iou,
                        ted_cov=1.0,
                        annot_cov=dom_ratio,
                        is_reviewed=is_rev,
                    )

            if m.protein_context.recommended_name:
                rec = strip_bracketed_tags(m.protein_context.recommended_name).strip()
                if rec:
                    tier = (
                        SpecificityTier.NOVEL_UNCHARACTERIZED
                        if is_low_information_name(rec)
                        else SpecificityTier.PROTEIN_FAMILY_CONTEXT
                    )
                    rec_dom_name = f"{rec} {pos_qualifier}" if "domain" not in rec.lower() else rec
                    _add_hit(
                        group_key=f"UP_REC:{normalize_domain_name(rec)}:{pos_qualifier}",
                        display_name=rec_dom_name,
                        short_name=None,
                        tier=tier,
                        source_db="uniprot_protein",
                        accession=m.protein_context.uniprot_acc,
                        member_id=member_id,
                        iou=ctx_iou,
                        ted_cov=1.0,
                        annot_cov=dom_ratio,
                        is_reviewed=is_rev,
                    )

    # Score all candidate terms
    candidates: list[ClusterCandidateTerm] = []
    tier_priority = [
        SpecificityTier.SPECIFIC_DOMAIN,
        SpecificityTier.SUPERFAMILY,
        SpecificityTier.PROTEIN_FAMILY_CONTEXT,
        SpecificityTier.TOPOLOGY_FOLD,
        SpecificityTier.NOVEL_UNCHARACTERIZED,
    ]

    for acc_obj in accumulators.values():
        member_hits = len(acc_obj.member_Best_iou)
        if member_hits == 0:
            continue
        best_raw_name = acc_obj.names_counter.most_common(1)[0][0]
        best_tier = min(acc_obj.tiers, key=lambda t: tier_priority.index(t))
        low_info = is_low_information_name(best_raw_name)
        if low_info:
            best_tier = SpecificityTier.NOVEL_UNCHARACTERIZED

        canonical = clean_display_domain_name(best_raw_name, best_tier)
        short = acc_obj.short_names[0] if acc_obj.short_names else make_short_label(canonical)

        freq = member_hits / total_members
        mean_iou = sum(acc_obj.member_Best_iou.values()) / member_hits
        mean_ted_cov = sum(acc_obj.member_best_ted_cov.values()) / member_hits
        mean_annot_cov = sum(acc_obj.member_best_annot_cov.values()) / member_hits

        w_tier = _TIER_WEIGHTS[best_tier]
        if low_info:
            p_quality = 0.25
        else:
            p_quality = 1.0
            if _LEGACY_PDB_HEADER_RE.search(best_raw_name) or _LEGACY_PDB_HEADER_RE.search(canonical):
                p_quality *= 0.65
            if is_sentence_like_description(best_raw_name):
                p_quality *= 0.40
            if is_bare_positional_name(best_raw_name):
                p_quality *= 0.50

        # Slight bonus when backed by multiple independent databases or reviewed Swiss-Prot entries
        multi_db_bonus = min(1.08, 1.0 + 0.02 * (len(acc_obj.source_dbs) - 1))
        reviewed_bonus = 1.04 if len(acc_obj.reviewed_members) > 0 else 1.0

        # Avoid artificial +0.45 floor when naming a single domain (total_members == 1 -> freq is trivially 1.0)
        if total_members == 1:
            overlap_term = 0.65 * mean_iou + 0.35 * mean_ted_cov
        else:
            overlap_term = 0.35 * freq + 0.45 * mean_iou + 0.20 * mean_ted_cov

        raw_score = (
            w_tier
            * overlap_term
            * p_quality
            * multi_db_bonus
            * reviewed_bonus
        )
        score = round(min(1.0, raw_score), 4)

        candidates.append(
            ClusterCandidateTerm(
                canonical_name=canonical,
                short_name=short,
                specificity_tier=best_tier,
                source_databases=sorted(acc_obj.source_dbs),
                supporting_accessions=sorted(acc_obj.accessions)[:10],
                member_hits=member_hits,
                cluster_frequency=round(freq, 4),
                mean_iou=round(mean_iou, 4),
                mean_ted_coverage=round(mean_ted_cov, 4),
                mean_annot_coverage=round(mean_annot_cov, 4),
                reviewed_hit_count=len(acc_obj.reviewed_members),
                is_low_information=low_info,
                composite_score=score,
            )
        )

    candidates.sort(
        key=lambda c: (
            not c.is_low_information,
            c.composite_score,
            c.cluster_frequency,
            c.mean_iou,
        ),
        reverse=True,
    )

    if candidates:
        winner = candidates[0]
        consensus_name = winner.canonical_name
        short_name = make_short_label(consensus_name)
        confidence = round(min(0.99, max(0.15, winner.composite_score)), 3)
        tier = winner.specificity_tier
        dbs_str = ", ".join(winner.source_databases)
        accs_str = ", ".join(winner.supporting_accessions[:5])
        rationale = (
            f"Selected '{consensus_name}' ({tier.value}) supported by {winner.member_hits}/{total_members} "
            f"cluster members ({winner.cluster_frequency:.0%}) across [{dbs_str}] "
            f"({accs_str}) with mean residue IoU={winner.mean_iou:.2f} and TED coverage={winner.mean_ted_coverage:.2f}."
        )
        if dominant_cath != "-":
            topo_str = f" ({cath_summary.topology_name})" if cath_summary and cath_summary.topology_name else ""
            if is_estimated_cath_code(dominant_cath):
                if novel_fold_candidate:
                    hint_msg = f"; 3D/SSE fold hint: {dominant_fold_hint}" if dominant_fold_hint else ""
                    rationale += f" Partial CATH assignment: {dominant_cath} (novel fold candidate{hint_msg})."
                else:
                    rationale += f" Partial CATH assignment: {dominant_cath}{topo_str}."
            elif tier == SpecificityTier.SPECIFIC_DOMAIN:
                rationale += f" Structural fold classification: CATH {dominant_cath}{topo_str}."
    else:
        # Complete novel / unannotated fallback
        if dominant_cath != "-":
            consensus_name = f"CATH {dominant_cath} unannotated structural domain"
        else:
            consensus_name = f"Uncharacterized {consensus_pos} structural domain"
        short_name = make_short_label(consensus_name)
        confidence = 0.20
        tier = SpecificityTier.NOVEL_UNCHARACTERIZED
        rationale = (
            f"No sequence or structural domain annotations overlapped the {total_members} TED member(s); "
            f"assigned positional fallback name ({consensus_pos}, CATH={dominant_cath})."
        )

    return ClusterNamingResult(
        cluster_id=cluster_id,
        member_count=total_members,
        consensus_name=consensus_name,
        short_name=short_name,
        confidence_score=confidence,
        specificity_tier=tier,
        naming_method="rule_based",
        rationale=rationale,
        dominant_cath_label=dominant_cath,
        cath_fold_hint=dominant_fold_hint,
        novel_fold_candidate=novel_fold_candidate,
        structure_quality_flags=structure_quality_flags,
        cath_summary=cath_summary,
        consensus_relative_position=consensus_pos,
        taxonomic_summary=taxonomic_summary,
        top_candidates=candidates[:10],
        members=members,
    )
