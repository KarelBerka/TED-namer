"""Cluster-level evidence aggregation, synonym unification, and 4-tier candidate scoring."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Optional

from ted_namer.cath_classifier import infer_cath_classification, is_estimated_cath_code
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
    r"\b(uncharacterized|uncharacterised|hypothetical|putative|duf\d+|domain of unknown function|"
    r"unnamed protein|predicted protein|whole genome shotgun|orf\b|upf\d+)\b",
    re.IGNORECASE,
)

_SPECIFIC_DOMAIN_DBS = {"cdd", "pfam", "smart", "prosite", "profile", "hamap", "ncbifam", "uniprot_feature"}
_SUPERFAMILY_DBS = {"cathgene3d", "ssf", "supfam"}

_TIER_WEIGHTS: dict[SpecificityTier, float] = {
    SpecificityTier.SPECIFIC_DOMAIN: 1.00,
    SpecificityTier.SUPERFAMILY: 0.85,
    SpecificityTier.PROTEIN_FAMILY_CONTEXT: 0.65,
    SpecificityTier.TOPOLOGY_FOLD: 0.45,
    SpecificityTier.NOVEL_UNCHARACTERIZED: 0.25,
}


def is_low_information_name(name: str) -> bool:
    """Return True if a candidate name is generic/uncharacterized (e.g. DUF, Uncharacterized protein)."""
    if not name or len(name.strip()) < 3:
        return True
    return bool(_LOW_INFO_PATTERNS.search(name))


def normalize_domain_name(name: str) -> str:
    """Normalize a domain/superfamily name for deduplication across databases."""
    cleaned = name.strip()
    cleaned = re.sub(r"[,_/()-]+", " ", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned).strip().lower()
    # Strip trailing generic words for key comparison
    cleaned = re.sub(r"\s+(homologous superfamily|superfamily|domain|family|fold|like)$", "", cleaned)
    return cleaned.strip()


def clean_display_domain_name(name: str, tier: SpecificityTier) -> str:
    """Standardize capitalization and suffix of a domain name for display, stripping PDB-header artifacts."""
    txt = re.sub(r"\s+", " ", name.strip().rstrip(".").replace(", ", " "))
    if not txt:
        return "Unassigned TED structural domain"

    # Strip legacy PDB-header artifacts from CATH/Gene3D representative titles
    txt = re.sub(r";\s*Chain\s+[A-Za-z0-9]+.*$", "", txt, flags=re.IGNORECASE)
    txt = re.sub(r"\.\s*Chain\s+[A-Za-z0-9]+.*$", "", txt, flags=re.IGNORECASE)
    txt = re.sub(r";\s*domain\s+(\d+)$", r" domain \1", txt, flags=re.IGNORECASE)
    txt = txt.strip().rstrip(";.,")

    # Preserve capitalization unless the entire string is ALL CAPS (like PANTHER 'CELLULAR TUMOR ANTIGEN P53')
    if txt.isupper() and len(txt) > 6:
        txt = txt.capitalize()

    low = txt.lower()
    if tier == SpecificityTier.SPECIFIC_DOMAIN:
        # Remove trailing 'superfamily' if present on a specific domain
        txt = re.sub(r"\s+superfamily$", "", txt, flags=re.IGNORECASE)
        if not any(
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
        ):
            txt = f"{txt} domain"
        elif txt.lower().endswith("motif") and "domain" not in txt.lower():
            txt = f"{txt} domain"
    elif tier == SpecificityTier.TOPOLOGY_FOLD:
        if "fold" not in low and "topology" not in low and "domain" not in low:
            txt = f"{txt} fold"
    return txt


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
    short = name.strip()
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
        and (db in _SPECIFIC_DOMAIN_DBS or db == "interpro")
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
    """Map an OverlappingAnnotation to its SpecificityTier based on source_db, type, and overlap."""
    if is_low_information_name(annot.name):
        return SpecificityTier.NOVEL_UNCHARACTERIZED

    atype = annot.annotation_type.lower()
    db = annot.source_db.lower()

    # Direct high-overlap domain or repeat match -> Tier 1 (Specific Domain)
    if atype in {"domain", "repeat"} or db in _SPECIFIC_DOMAIN_DBS:
        if annot.iou >= 0.35 or (annot.ted_coverage >= 0.50 and annot.annot_coverage >= 0.50):
            return SpecificityTier.SPECIFIC_DOMAIN

    # Structural lobe/subdomain of a single curated multi-lobe domain (e.g. Protein kinase N-lobe)
    if _is_curated_subdomain_lobe(annot):
        return SpecificityTier.SPECIFIC_DOMAIN

    # Supra-domain full-protein family matches or tiny sub-domain motifs belong to Tier 3 (Context)
    if (
        atype in {"family", "conserved_site", "motif", "active_site", "binding_site"}
        or annot.overlap_category in {OverlapCategory.SUPRA_DOMAIN, OverlapCategory.SUB_DOMAIN}
    ):
        # Allow a motif/region to be Tier 1 only if it genuinely covers most of the TED domain
        if atype in {"motif", "region"} and annot.iou >= 0.60 and annot.ted_coverage >= 0.65:
            return SpecificityTier.SPECIFIC_DOMAIN
        return SpecificityTier.PROTEIN_FAMILY_CONTEXT

    if atype == "homologous_superfamily" or db in _SUPERFAMILY_DBS:
        if annot.iou >= 0.30 or annot.ted_coverage >= 0.50:
            return SpecificityTier.SUPERFAMILY

    return SpecificityTier.PROTEIN_FAMILY_CONTEXT


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
        )

    # Ensure every member has inferred_cath_label populated (even in offline unit tests)
    for m in members:
        if not m.domain.inferred_cath_label:
            inferred_code, _num_code, est_src = infer_cath_classification(
                member=m.domain,
                overlapping_annotations=m.overlapping_annotations,
            )
            m.domain.inferred_cath_label = inferred_code
            m.domain.cath_estimation_source = est_src

    # Determine dominant CATH label (including T-prefixed estimated codes) and consensus relative position
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
        elif source_db in {"cdd", "pfam", "smart", "prosite", "profile", "uniprot_feature"}:
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

        # 1. Overlapping sequence/structure annotations
        for annot in m.overlapping_annotations:
            tier = determine_annotation_tier(annot)
            ipr_group_id = annot.integrated_ipr or (
                annot.accession if annot.source_db == "interpro" and annot.accession.startswith("IPR") else None
            )

            disp = annot.name
            # Replace raw G3DSA/SSF legacy PDB titles with the curated integrated InterPro superfamily name
            if annot.integrated_ipr_name and (
                annot.source_db in {"cathgene3d", "ssf"}
                or any(
                    tok in disp
                    for tok in ("; domain", "domain 1", "domain 2", "Chain ", "Virus protein", "subunit ")
                )
            ):
                disp = annot.integrated_ipr_name

            eff_iou = annot.iou
            if _is_curated_subdomain_lobe(annot):
                sub_pos = _compute_subdomain_position(m.domain.chopping, annot.chopping)
                base_dom = re.sub(r"\s+domain$", "", disp, flags=re.IGNORECASE).strip()
                disp = f"{base_dom} {sub_pos}"
                gkey = f"SUBLOBE:{ipr_group_id or normalize_domain_name(base_dom)}:{sub_pos}"
                eff_iou = round(0.55 * annot.ted_coverage + 0.45 * annot.annot_coverage, 4)
            else:
                if ipr_group_id and tier in {SpecificityTier.SPECIFIC_DOMAIN, SpecificityTier.SUPERFAMILY}:
                    gkey = f"IPR:{ipr_group_id}:{tier.value}"
                else:
                    gkey = f"NORM:{normalize_domain_name(disp)}:{tier.value}"

                if tier == SpecificityTier.PROTEIN_FAMILY_CONTEXT and "domain" not in disp.lower():
                    # Qualify supra-domain family hit by the TED domain's relative position
                    base_fam = re.sub(r"\s+family$", "", disp, flags=re.IGNORECASE)
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

        # 2. CATH resolution evidence (H-level superfamily and T-level fold)
        if m.cath_info:
            if m.cath_info.superfamily_name:
                ipr_key = m.cath_info.integrated_ipr or f"CATH_H:{m.cath_info.cath_label}"
                gkey = f"IPR:{ipr_key}:{SpecificityTier.SUPERFAMILY.value}"
                _add_hit(
                    group_key=gkey,
                    display_name=m.cath_info.superfamily_name,
                    short_name=None,
                    tier=SpecificityTier.SUPERFAMILY,
                    source_db="cath_h",
                    accession=f"CATH:{m.cath_info.cath_label}",
                    member_id=member_id,
                    iou=0.85,
                    ted_cov=0.90,
                    annot_cov=0.85,
                    is_reviewed=is_rev,
                )
            if m.cath_info.topology_name and m.cath_info.topology_code:
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
                rec = m.protein_context.recommended_name.strip()
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
        p_quality = 0.25 if low_info else 1.0
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
                rationale += f" Estimated CATH classification: {dominant_cath}{topo_str}."
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
        cath_summary=cath_summary,
        consensus_relative_position=consensus_pos,
        taxonomic_summary=taxonomic_summary,
        top_candidates=candidates[:10],
        members=members,
    )
