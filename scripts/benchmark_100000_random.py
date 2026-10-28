"""Expand the TED domain benchmark to 100,000 unique TED domains (100 batches of 1,000 domains)
and perform comprehensive statistical and blind cross-validation evaluation across all 100,000 domains."""

from __future__ import annotations

import asyncio
import csv
import json
import random
import re
import sqlite3
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import httpx

from ted_namer.aggregator import aggregate_cluster_evidence, normalize_domain_name
from ted_namer.cache import ApiCache
from ted_namer.cath_classifier import _G3DSA_RE, classify_cath_with_rules
from ted_namer.clients.ted import TedClient
from ted_namer.export import to_delimited as to_tsv_string
from ted_namer.models import (
    ClusterNamingResult,
    EnrichedTedDomain,
    OverlapCategory,
    OverlappingAnnotation,
    ProteinContext,
    SpecificityTier,
    TedDomainMember,
)
from ted_namer.pipeline import TedNamingPipeline

# Broad taxonomic & length queries across all kingdoms of life (Swiss-Prot + TrEMBL)
# Prioritizing multi-domain proteins (length >= 220) alongside single-domain proteins (100-220)
_UNIPROT_100K_QUERIES = [
    # Reviewed Swiss-Prot across Kingdoms and Length Bins
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606) AND (length:[100 TO 300])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606) AND (length:[301 TO 500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606) AND (length:[501 TO 800])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606) AND (length:[801 TO 1300])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606) AND (length:[1301 TO 2700])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090) AND (length:[100 TO 350])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090) AND (length:[351 TO 650])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090) AND (length:[651 TO 1200])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090) AND (length:[1201 TO 2700])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10116) AND (length:[200 TO 1500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9913) AND (length:[200 TO 1500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7955) AND (length:[200 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7227) AND (length:[200 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:6239) AND (length:[200 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702) AND (length:[150 TO 450])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702) AND (length:[451 TO 900])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702) AND (length:[901 TO 2500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:39947) AND (length:[200 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:559292) AND (length:[200 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:284812) AND (length:[200 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4751) AND (length:[300 TO 900])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4751) AND (length:[901 TO 2500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:33090) AND (length:[300 TO 900])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:33208) AND (length:[400 TO 900])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:33208) AND (length:[901 TO 2500])",
    # Reviewed Bacteria & Archaea across Length Bins
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2) AND (length:[150 TO 300])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2) AND (length:[301 TO 450])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2) AND (length:[451 TO 650])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2) AND (length:[651 TO 1000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2) AND (length:[1001 TO 2500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:1224) AND (length:[250 TO 600])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:1224) AND (length:[601 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:1239) AND (length:[250 TO 600])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:1239) AND (length:[601 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:201174) AND (length:[250 TO 1500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2157) AND (length:[150 TO 400])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2157) AND (length:[401 TO 2000])",
    # Unreviewed TrEMBL across Diverse Taxa
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:1224) AND (length:[250 TO 800])",
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:1239) AND (length:[250 TO 800])",
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:201174) AND (length:[250 TO 800])",
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:976) AND (length:[250 TO 800])",
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:1117) AND (length:[250 TO 800])",
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:28890) AND (length:[250 TO 800])",
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:4751) AND (length:[300 TO 1000])",
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:33090) AND (length:[300 TO 1000])",
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:33208) AND (length:[300 TO 1000])",
]

_NEXT_LINK_RE = re.compile(r'<([^>]+)>;\s*rel="next"')


def result_to_compact(r: ClusterNamingResult) -> dict[str, Any]:
    """Convert a ClusterNamingResult into a compact 1-line JSONL dictionary retaining all validation fields."""
    m = r.members[0]
    return {
        "cluster_id": r.cluster_id,
        "ted_id": m.domain.ted_id,
        "uniprot_acc": m.domain.uniprot_acc,
        "chopping": m.domain.chopping,
        "nres_domain": m.domain.nres_domain,
        "num_helix": m.domain.num_helix,
        "num_strand": m.domain.num_strand,
        "num_turn": m.domain.num_turn,
        "norm_rg": m.domain.norm_rg,
        "packing_density": m.domain.packing_density,
        "plddt": m.domain.plddt,
        "cath_label": m.domain.cath_label,
        "inferred_cath_label": m.domain.inferred_cath_label,
        "cath_estimation_source": m.domain.cath_estimation_source,
        "consensus_name": r.consensus_name,
        "short_name": r.short_name,
        "specificity_tier": r.specificity_tier.value,
        "confidence_score": r.confidence_score,
        "dominant_cath_label": r.dominant_cath_label,
        "top_dbs": r.top_candidates[0].source_databases if r.top_candidates else [],
        "top_accs": r.top_candidates[0].supporting_accessions if r.top_candidates else [],
        "is_reviewed": bool(m.protein_context and m.protein_context.is_reviewed),
        "protein_length": m.protein_context.protein_length if m.protein_context else 0,
        "recommended_name": m.protein_context.recommended_name if m.protein_context else None,
        "similarity_families": m.protein_context.similarity_families if m.protein_context else [],
        "relative_position": m.relative_position,
        "cath_h_name": m.cath_info.superfamily_name if m.cath_info else None,
        "cath_h_ipr": m.cath_info.integrated_ipr if m.cath_info else None,
        "cath_t_name": m.cath_info.topology_name if m.cath_info else None,
        "cath_t_code": m.cath_info.topology_code if m.cath_info else None,
        "annots": [
            {
                "db": a.source_db,
                "acc": a.accession,
                "name": a.name,
                "ipr": a.integrated_ipr,
                "ipr_name": a.integrated_ipr_name,
                "type": a.annotation_type,
                "chop": a.chopping,
                "iou": a.iou,
                "tcov": a.ted_coverage,
                "acov": a.annot_coverage,
                "cat": a.overlap_category.value,
            }
            for a in m.overlapping_annotations
        ],
    }


def harvest_cached_ted_domains(cache_path: Path, seen_ted_ids: set[str]) -> list[TedDomainMember]:
    """Extract all TED domains already cached in ted_uniprot_summary_v2 that are not in seen_ted_ids."""
    harvested: list[TedDomainMember] = []
    conn = sqlite3.connect(cache_path)
    for (payload_str,) in conn.execute(
        "SELECT payload FROM api_cache WHERE namespace = 'ted_uniprot_summary_v2'"
    ):
        try:
            items = json.loads(payload_str)
            if not isinstance(items, list):
                continue
            for raw_item in items:
                dom = TedDomainMember.model_validate(raw_item)
                if dom.ted_id not in seen_ted_ids:
                    seen_ted_ids.add(dom.ted_id)
                    harvested.append(dom)
        except Exception:
            continue
    conn.close()
    return harvested


async def sample_remaining_ted_domains(
    cache: ApiCache,
    seen_ted_ids: set[str],
    target_new_count: int,
    seed: int = 20261028,
) -> list[TedDomainMember]:
    """Sample remaining unique TED domains by querying diverse UniProt buckets and collecting all TED domains per protein."""
    rng = random.Random(seed)
    headers = {"User-Agent": "ted-namer-benchmark-100000/0.1.0"}
    candidate_uniprots: list[str] = []

    limits = httpx.Limits(max_connections=120, max_keepalive_connections=60)
    async with httpx.AsyncClient(headers=headers, follow_redirects=True, timeout=30.0, limits=limits) as http:
        q_sem = asyncio.Semaphore(14)

        async def _run_query(q: str) -> list[str]:
            accs_local: list[str] = []
            url = "https://rest.uniprot.org/uniprotkb/search"
            params: dict[str, str] | None = {
                "query": q,
                "format": "json",
                "size": "500",
                "fields": "accession",
            }
            async with q_sem:
                for _page in range(3):  # Up to 3 pages (1,500 accessions) per query = ~65,000 candidates
                    if not url:
                        break
                    try:
                        resp = await http.get(url, params=params)
                        params = None
                        if resp.status_code == 200:
                            results = resp.json().get("results") or []
                            accs_local.extend(
                                r["primaryAccession"] for r in results if r.get("primaryAccession")
                            )
                            link_hdr = resp.headers.get("Link", "")
                            m = _NEXT_LINK_RE.search(link_hdr)
                            url = m.group(1) if m else ""
                        else:
                            break
                    except Exception:
                        break
            return accs_local

        query_lists = await asyncio.gather(*(_run_query(q) for q in _UNIPROT_100K_QUERIES))
        for accs in query_lists:
            rng.shuffle(accs)
            candidate_uniprots.extend(accs)

        unique_uniprots = list(dict.fromkeys(candidate_uniprots))
        rng.shuffle(unique_uniprots)
        print(f"Collected {len(unique_uniprots)} unique candidate UniProt accessions.", flush=True)

        ted_client = TedClient(http_client=http, cache=cache)
        sem = asyncio.Semaphore(60)

        async def _get_ted(acc: str) -> list[TedDomainMember]:
            async with sem:
                return await ted_client.get_domains_for_uniprot(acc)

        sampled: list[TedDomainMember] = []
        batch_size = 250
        for i in range(0, len(unique_uniprots), batch_size):
            batch = unique_uniprots[i : i + batch_size]
            batch_res = await asyncio.gather(*(_get_ted(a) for a in batch))
            for dom_list in batch_res:
                if not dom_list:
                    continue
                # Include all TED domains on each sampled protein (up to 6) to maximize protein-level API reuse
                for d in dom_list[:6]:
                    if d.ted_id not in seen_ted_ids:
                        seen_ted_ids.add(d.ted_id)
                        sampled.append(d)
            cache.flush()
            if len(sampled) >= target_new_count:
                break
            if (i // batch_size) % 4 == 0:
                print(
                    f"  Sampled {len(sampled)}/{target_new_count} additional TED domains "
                    f"(from {i + len(batch)} UniProt proteins)...",
                    flush=True,
                )

    return sampled[:target_new_count]


def evaluate_100k_jsonl(jsonl_path: Path) -> dict[str, Any]:
    """Run full statistical & 5-fold CV accuracy evaluation across all 100,000 domains in jsonl_path."""
    print(f"\nStreaming {jsonl_path} for statistical and 5-fold cross-validation evaluation...", flush=True)
    records: list[dict[str, Any]] = []
    with open(jsonl_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))

    n_total = len(records)
    print(f"Loaded {n_total} compact domain records.", flush=True)

    # 1. Compute overall and per-10k-decile statistics
    def _subset_stats(sub: list[dict[str, Any]]) -> dict[str, Any]:
        n = len(sub)
        if n == 0:
            return {}
        tiers = Counter(r["specificity_tier"] for r in sub)
        sources = Counter(r.get("cath_estimation_source") or "-" for r in sub)
        mean_conf = sum(r["confidence_score"] for r in sub) / n
        tier_confs: dict[str, list[float]] = defaultdict(list)
        for r in sub:
            tier_confs[r["specificity_tier"]].append(r["confidence_score"])

        official_cath = sources.get("ted_official", 0)
        t_level_cath = sum(v for k, v in sources.items() if k.startswith("ted_t_"))
        g3dsa_cath = sources.get("g3dsa_overlap", 0)
        profile_cath = sources.get("profile_cath_map", 0)
        direct_rules = sum(v for k, v in sources.items() if k.startswith("rule_"))

        return {
            "count": n,
            "mean_confidence": round(mean_conf, 4),
            "mean_conf_by_tier": {k: round(sum(v) / len(v), 4) for k, v in tier_confs.items()},
            "tier_specific_domain": round(tiers.get("specific_domain", 0) / n, 5),
            "tier_superfamily": round(tiers.get("superfamily", 0) / n, 5),
            "tier_specific_or_superfamily": round(
                (tiers.get("specific_domain", 0) + tiers.get("superfamily", 0)) / n, 5
            ),
            "tier_family_context": round(tiers.get("protein_family_context", 0) / n, 5),
            "tier_topology_fold": round(tiers.get("topology_fold", 0) / n, 5),
            "cath_official": round(official_cath / n, 5),
            "cath_t_level": round(t_level_cath / n, 5),
            "cath_profile_map": round(profile_cath / n, 5),
            "cath_g3dsa": round(g3dsa_cath / n, 5),
            "cath_direct_3d_rules": round(direct_rules / n, 5),
            "conf_below_50_pct": round(sum(1 for r in sub if r["confidence_score"] < 0.50) / n, 5),
            "low_plddt_pct": round(sum(1 for r in sub if r.get("plddt", 100.0) < 60.0) / n, 5),
            "sublobe_count": sum(1 for r in sub if "subdomain" in r["consensus_name"].lower()),
        }

    deciles: dict[str, Any] = {}
    for d_idx in range(n_total // 10000):
        start = d_idx * 10000
        end = (d_idx + 1) * 10000
        deciles[f"Decile_{d_idx + 1}_{start + 1}_to_{end}"] = _subset_stats(records[start:end])

    overall_stats = _subset_stats(records)

    # 2. Blind Hold-Out CATH Accuracy & 5-Fold Cross-Validation on all Official 4-Level CATH Domains
    official_entries: list[tuple[list[str], dict[str, Any]]] = []
    for r in records:
        raw = (r.get("cath_label") or "-").split(",")[0].strip()
        parts = [p.strip() for p in raw.split(".") if p.strip() and p.strip() != "-"]
        if len(parts) >= 4:
            official_entries.append((parts[:4], r))

    n_off = len(official_entries)
    print(f"Running 5-fold CV & blind hold-out CATH evaluation on {n_off} official 4-level CATH domains...", flush=True)

    num_folds = 5
    fold_maps: list[dict[str, tuple[list[str], int]]] = []
    for fold_idx in range(num_folds):
        acc_to_cath: dict[str, Counter[str]] = defaultdict(Counter)
        for i, (tp, rec) in enumerate(official_entries):
            if i % num_folds == fold_idx:
                continue
            cath_str = ".".join(tp)
            for a in rec.get("annots", []):
                if a["db"] == "cathgene3d" or str(a["acc"]).upper().startswith("G3DSA:"):
                    continue
                if a["iou"] >= 0.40 or (a["tcov"] >= 0.65 and a["acov"] >= 0.45):
                    for key in (a["acc"], a.get("ipr")):
                        if key:
                            acc_to_cath[str(key)][cath_str] += 1

        fmap: dict[str, tuple[list[str], int]] = {}
        for k, cnts in acc_to_cath.items():
            tot = sum(cnts.values())
            top_cath, top_n = cnts.most_common(1)[0]
            h_purity = top_n / tot
            cat_cnts = Counter(".".join(c.split(".")[:3]) for c, n in cnts.items() for _ in range(n))
            top_cat, top_cat_n = cat_cnts.most_common(1)[0]
            cat_purity = top_cat_n / tot
            if cat_purity < 0.85:
                continue
            k_depth = 3 if (tot >= 2 and cat_purity >= 0.90) else 1
            code = top_cath if h_purity >= 0.85 else f"{top_cat}.10"
            fmap[k] = (code.split("."), k_depth)
        fold_maps.append(fmap)

    cas_a_c = cas_a_ca = cas_a_cat = cas_a_cath = 0
    cas_a_direct_claimed: Counter[int] = Counter()
    cas_a_direct_ok: Counter[int] = Counter()
    cas_b_c = cas_b_ca = cas_b_cat = cas_b_cath = 0
    cas_b_prof_hits = cas_b_prof_c = cas_b_prof_ca = cas_b_prof_cat = cas_b_prof_cath = 0
    cas_c_c = cas_c_ca = cas_c_cat = 0
    cas_c_pure_ab_tot = cas_c_pure_ab_ok = 0

    for i, (tp, rec) in enumerate(official_entries):
        fmap = fold_maps[i % num_folds]
        annots = rec.get("annots", [])

        # Check G3DSA
        best_g3d = None
        best_g3d_iou = best_g3d_tcov = 0.0
        best_g3d_sc = -1.0
        for a in annots:
            if a["db"] == "cathgene3d" or str(a["acc"]).upper().startswith("G3DSA:"):
                if a["iou"] >= 0.25 or a["tcov"] >= 0.45:
                    m_re = _G3DSA_RE.match(str(a["acc"]).strip())
                    if m_re:
                        sc = a["iou"] + 0.5 * a["tcov"]
                        if sc > best_g3d_sc:
                            best_g3d_sc = sc
                            best_g3d_iou = a["iou"]
                            best_g3d_tcov = a["tcov"]
                            best_g3d = [m_re.group(1), m_re.group(2), m_re.group(3), m_re.group(4)]

        # Check 5-fold CV Profile-to-CATH
        best_prof = None
        best_prof_depth = 1
        best_prof_iou = 0.0
        best_prof_sc = -1.0
        for a in annots:
            if a["db"] == "cathgene3d" or str(a["acc"]).upper().startswith("G3DSA:"):
                continue
            if a["iou"] >= 0.40 or (a["tcov"] >= 0.65 and a["acov"] >= 0.45):
                for key in (a["acc"], a.get("ipr")):
                    if key and str(key) in fmap:
                        p_parts, p_depth = fmap[str(key)]
                        db_boost = 1.15 if a["db"] == "ssf" else 1.0
                        sc = (a["iou"] + 0.3 * a["tcov"]) * db_boost * (1.1 if p_depth >= 3 else 1.0)
                        if sc > best_prof_sc:
                            best_prof_sc = sc
                            best_prof = p_parts
                            best_prof_depth = p_depth
                            best_prof_iou = a["iou"]

        # Pure SSE rules
        dom_obj = TedDomainMember(
            ted_id=rec["ted_id"],
            uniprot_acc=rec["uniprot_acc"],
            chopping=rec["chopping"],
            nres_domain=rec["nres_domain"],
            num_helix=rec.get("num_helix", 0),
            num_strand=rec.get("num_strand", 0),
            num_turn=rec.get("num_turn", 0),
            norm_rg=rec.get("norm_rg", 0.3),
            packing_density=rec.get("packing_density", 10.0),
            plddt=rec.get("plddt", 85.0),
            cath_label="-",
        )
        rule_parts, rule_depth, _ = classify_cath_with_rules(dom_obj, pdb_text=None)
        if rule_parts[0] == tp[0]:
            cas_c_c += 1
        if rule_parts[:2] == tp[:2]:
            cas_c_ca += 1
        if rule_parts[:3] == tp[:3]:
            cas_c_cat += 1
        if dom_obj.num_strand == 0 or dom_obj.num_helix == 0:
            cas_c_pure_ab_tot += 1
            if rule_parts[0] == tp[0]:
                cas_c_pure_ab_ok += 1

        # Cascade B (No TED CATH, No G3DSA)
        if best_prof:
            pred_b = best_prof
            cas_b_prof_hits += 1
            if pred_b[0] == tp[0]:
                cas_b_prof_c += 1
            if pred_b[:2] == tp[:2]:
                cas_b_prof_ca += 1
            if pred_b[:3] == tp[:3]:
                cas_b_prof_cat += 1
            if pred_b[:4] == tp[:4]:
                cas_b_prof_cath += 1
        else:
            pred_b = rule_parts
        if pred_b[0] == tp[0]:
            cas_b_c += 1
        if pred_b[:2] == tp[:2]:
            cas_b_ca += 1
        if pred_b[:3] == tp[:3]:
            cas_b_cat += 1
        if pred_b[:4] == tp[:4]:
            cas_b_cath += 1

        # Cascade A (Full fallback cascade when TED CATH is masked)
        if best_g3d:
            pred_a = best_g3d
            k_depth = 4 if (best_g3d_iou >= 0.65 and best_g3d_tcov >= 0.75) else (3 if (best_g3d_iou >= 0.40 or best_g3d_tcov >= 0.65) else 1)
        elif best_prof:
            pred_a = best_prof
            k_depth = best_prof_depth if best_prof_iou >= 0.45 else 1
        else:
            pred_a = rule_parts
            k_depth = rule_depth

        if pred_a[0] == tp[0]:
            cas_a_c += 1
        if pred_a[:2] == tp[:2]:
            cas_a_ca += 1
        if pred_a[:3] == tp[:3]:
            cas_a_cat += 1
        if pred_a[:4] == tp[:4]:
            cas_a_cath += 1

        for d_lvl in range(1, k_depth + 1):
            cas_a_direct_claimed[d_lvl] += 1
            if pred_a[:d_lvl] == tp[:d_lvl]:
                cas_a_direct_ok[d_lvl] += 1

    # 3. Leave-Pfam-Out Naming Concordance across all eligible domains in the 100k dataset
    lodo_tot = lodo_exact = lodo_semantic = 0
    pdb_noise_re = re.compile(r"(;\s*domain\s+\d+|\.\s*Chain\s+[A-Z0-9]+|Vaccinia Virus protein VP39)", re.IGNORECASE)
    pdb_noise_cnt = 0

    for rec in records:
        if pdb_noise_re.search(rec["consensus_name"]):
            pdb_noise_cnt += 1
        if rec["specificity_tier"] == "specific_domain" and "pfam" in rec.get("top_dbs", []):
            annots = rec.get("annots", [])
            non_pfam = [a for a in annots if a["db"] != "pfam" and a["type"] == "domain" and a["iou"] >= 0.35]
            if non_pfam:
                lodo_tot += 1
                # Reconstruct EnrichedTedDomain without pfam
                dom_obj = TedDomainMember(
                    ted_id=rec["ted_id"],
                    uniprot_acc=rec["uniprot_acc"],
                    chopping=rec["chopping"],
                    nres_domain=rec["nres_domain"],
                    num_helix=rec.get("num_helix", 0),
                    num_strand=rec.get("num_strand", 0),
                    num_turn=rec.get("num_turn", 0),
                    norm_rg=rec.get("norm_rg", 0.3),
                    packing_density=rec.get("packing_density", 10.0),
                    plddt=rec.get("plddt", 85.0),
                    cath_label=rec.get("cath_label", "-"),
                    inferred_cath_label=rec.get("inferred_cath_label"),
                    cath_estimation_source=rec.get("cath_estimation_source"),
                )
                ov_list = [
                    OverlappingAnnotation(
                        source_db=a["db"],
                        accession=a["acc"],
                        name=a["name"],
                        integrated_ipr=a.get("ipr"),
                        integrated_ipr_name=a.get("ipr_name"),
                        annotation_type=a["type"],
                        chopping=a["chop"],
                        iou=a["iou"],
                        ted_coverage=a["tcov"],
                        annot_coverage=a["acov"],
                        overlap_category=OverlapCategory(a["cat"]),
                    )
                    for a in annots
                    if a["db"] != "pfam"
                ]
                p_ctx = ProteinContext(
                    uniprot_acc=rec["uniprot_acc"],
                    recommended_name=rec.get("recommended_name"),
                    is_reviewed=rec.get("is_reviewed", False),
                    protein_length=rec.get("protein_length", 0),
                    similarity_families=rec.get("similarity_families", []),
                )
                em = EnrichedTedDomain(
                    domain=dom_obj,
                    protein_context=p_ctx,
                    overlapping_annotations=ov_list,
                    relative_position=rec.get("relative_position", "domain"),
                )
                res_ho = aggregate_cluster_evidence(rec["cluster_id"], [em])
                norm_full = normalize_domain_name(rec["consensus_name"])
                norm_ho = normalize_domain_name(res_ho.consensus_name)
                iprs_full = {a for a in rec.get("top_accs", []) if a.startswith("IPR")}
                iprs_ho = (
                    {a for a in res_ho.top_candidates[0].supporting_accessions if a.startswith("IPR")}
                    if res_ho.top_candidates
                    else set()
                )
                if norm_full == norm_ho or (iprs_full and (iprs_full & iprs_ho)):
                    lodo_exact += 1
                    lodo_semantic += 1
                else:
                    t1 = {w for w in norm_full.split() if len(w) >= 3}
                    t2 = {w for w in norm_ho.split() if len(w) >= 3}
                    if t1 and t2 and (t1 & t2):
                        lodo_semantic += 1

    return {
        "overall_100k_stats": overall_stats,
        "deciles_10k": deciles,
        "cath_holdout_validation": {
            "n_official": n_off,
            "cascade_a_full": {
                "c_acc": round(cas_a_c / max(1, n_off), 5),
                "ca_acc": round(cas_a_ca / max(1, n_off), 5),
                "cat_acc": round(cas_a_cat / max(1, n_off), 5),
                "cath_acc": round(cas_a_cath / max(1, n_off), 5),
                "direct_level_precision": {
                    str(d): round(cas_a_direct_ok[d] / max(1, cas_a_direct_claimed[d]), 5)
                    for d in sorted(cas_a_direct_claimed)
                },
                "direct_level_counts": dict(cas_a_direct_claimed),
            },
            "cascade_b_no_g3dsa": {
                "c_acc": round(cas_b_c / max(1, n_off), 5),
                "ca_acc": round(cas_b_ca / max(1, n_off), 5),
                "cat_acc": round(cas_b_cat / max(1, n_off), 5),
                "cath_acc": round(cas_b_cath / max(1, n_off), 5),
                "profile_map_coverage": round(cas_b_prof_hits / max(1, n_off), 5),
                "profile_map_c_acc": round(cas_b_prof_c / max(1, cas_b_prof_hits), 5),
                "profile_map_ca_acc": round(cas_b_prof_ca / max(1, cas_b_prof_hits), 5),
                "profile_map_cat_acc": round(cas_b_prof_cat / max(1, cas_b_prof_hits), 5),
                "profile_map_cath_acc": round(cas_b_prof_cath / max(1, cas_b_prof_hits), 5),
            },
            "cascade_c_pure_3d_sse": {
                "c_acc": round(cas_c_c / max(1, n_off), 5),
                "pure_alpha_or_beta_c_acc": round(cas_c_pure_ab_ok / max(1, cas_c_pure_ab_tot), 5),
                "ca_acc": round(cas_c_ca / max(1, n_off), 5),
                "cat_acc": round(cas_c_cat / max(1, n_off), 5),
            },
        },
        "naming_validation": {
            "lodo_n": lodo_tot,
            "lodo_exact_concordance": round(lodo_exact / max(1, lodo_tot), 5),
            "lodo_semantic_concordance": round(lodo_semantic / max(1, lodo_tot), 5),
            "legacy_pdb_noise_count": pdb_noise_cnt,
            "legacy_pdb_noise_pct": round(pdb_noise_cnt / max(1, n_total), 6),
        },
    }


async def main() -> None:
    examples_dir = Path("examples")
    examples_dir.mkdir(exist_ok=True)
    cache_path = Path(".ted_namer_cache.sqlite")
    cache = ApiCache(db_path=cache_path)

    jsonl_path = examples_dir / "random_100000_compact.jsonl"
    tsv_path = examples_dir / "random_100000_results.tsv"

    seen_ted_ids: set[str] = set()
    existing_count = 0

    # 1. Seed jsonl_path and tsv_path from random_20000_results.json if jsonl_path does not exist yet
    if not jsonl_path.exists():
        b20k_path = examples_dir / "random_20000_results.json"
        if b20k_path.exists():
            print("Seeding first 20,000 domains from examples/random_20000_results.json...", flush=True)
            raw_20k = json.loads(b20k_path.read_text(encoding="utf-8"))
            res_20k = [ClusterNamingResult.model_validate(x) for x in raw_20k]
            with open(jsonl_path, "w", encoding="utf-8") as jf:
                for r in res_20k:
                    c = result_to_compact(r)
                    seen_ted_ids.add(c["ted_id"])
                    jf.write(json.dumps(c) + "\n")
            tsv_path.write_text(to_tsv_string(res_20k), encoding="utf-8")
            existing_count = len(res_20k)
            del raw_20k, res_20k
            print(f"Seeded {existing_count} domains into {jsonl_path} and {tsv_path}.", flush=True)
    else:
        with open(jsonl_path, encoding="utf-8") as jf:
            for line in jf:
                if line.strip():
                    obj = json.loads(line)
                    seen_ted_ids.add(obj["ted_id"])
                    existing_count += 1
        print(f"Resuming from {existing_count} already processed domains in {jsonl_path}.", flush=True)

    needed = 100000 - existing_count
    if needed > 0:
        sampled_100k_pool_path = examples_dir / "sampled_domains_pool_100000.json"
        pool_domains: list[TedDomainMember] = []
        if sampled_100k_pool_path.exists():
            raw_pool = json.loads(sampled_100k_pool_path.read_text(encoding="utf-8"))
            for x in raw_pool:
                d = TedDomainMember.model_validate(x)
                if d.ted_id not in seen_ted_ids:
                    seen_ted_ids.add(d.ted_id)
                    pool_domains.append(d)
            print(f"Loaded {len(pool_domains)} pre-sampled unseen TED domains from {sampled_100k_pool_path}.", flush=True)

        if len(pool_domains) < needed:
            # First harvest already-cached TED domains from ted_uniprot_summary_v2
            harvested = harvest_cached_ted_domains(cache_path, seen_ted_ids)
            pool_domains.extend(harvested)
            print(
                f"Harvested {len(harvested)} already-cached unseen TED domains from SQLite cache "
                f"(pool now has {len(pool_domains)}/{needed}).",
                flush=True,
            )

        if len(pool_domains) < needed:
            still_needed = needed - len(pool_domains)
            print(f"Sampling {still_needed} additional unique TED domains from UniProt/TED...", flush=True)
            t_samp = time.perf_counter()
            newly_sampled = await sample_remaining_ted_domains(
                cache=cache,
                seen_ted_ids=seen_ted_ids,
                target_new_count=still_needed,
            )
            pool_domains.extend(newly_sampled)
            print(
                f"Sampled {len(newly_sampled)} new TED domains in {time.perf_counter() - t_samp:.1f}s "
                f"(total pool: {len(pool_domains)}).",
                flush=True,
            )
            sampled_100k_pool_path.write_text(
                json.dumps([d.model_dump() for d in pool_domains]),
                encoding="utf-8",
            )

        target_domains = pool_domains[:needed]
        pipeline = TedNamingPipeline(
            use_cache=True,
            concurrency=85,
            use_llm=False,
            fetch_pdb=False,
        )

        batch_size = 2000
        num_chunks = (len(target_domains) + batch_size - 1) // batch_size
        start_total = existing_count

        for c_idx in range(num_chunks):
            chunk_doms = target_domains[c_idx * batch_size : (c_idx + 1) * batch_size]
            cur_start = start_total + c_idx * batch_size + 1
            cur_end = start_total + c_idx * batch_size + len(chunk_doms)
            b_label = f"Domains {cur_start:,}-{cur_end:,} / 100,000"
            print(f"\n=======================================================", flush=True)
            print(f"Processing {b_label} ({len(chunk_doms)} domains)...", flush=True)
            print(f"=======================================================", flush=True)

            pre_map = {
                f"Q_{dom.uniprot_acc}_{dom.short_ted_num}": [dom]
                for dom in chunk_doms
            }
            c_t0 = time.perf_counter()
            last_log = 0

            def _prog(done: int, total: int) -> None:
                nonlocal last_log
                if done - last_log >= 500 or done == total:
                    last_log = done
                    el = time.perf_counter() - c_t0
                    rate = done / max(0.1, el)
                    print(
                        f"  [{b_label}] {done:4d}/{total:4d} ({done/total*100:5.1f}%) | "
                        f"{el:5.1f}s | {rate:5.1f} dom/s",
                        flush=True,
                    )

            chunk_results = await pipeline.name_clusters(
                cluster_map={},
                pre_resolved_domains=pre_map,
                progress_callback=_prog,
            )
            pipeline.cache.flush()
            c_el = time.perf_counter() - c_t0
            print(
                f"  Finished {b_label} in {c_el:.1f}s ({len(chunk_results)/max(0.1, c_el):.1f} dom/s).",
                flush=True,
            )

            # Append to compact JSONL and TSV
            with open(jsonl_path, "a", encoding="utf-8") as jf:
                for r in chunk_results:
                    jf.write(json.dumps(result_to_compact(r)) + "\n")

            tsv_chunk = to_tsv_string(chunk_results)
            tsv_lines = tsv_chunk.splitlines()
            if tsv_path.exists() and tsv_path.stat().st_size > 0:
                with open(tsv_path, "a", encoding="utf-8") as tf:
                    for ln in tsv_lines[1:]:
                        tf.write(ln + "\n")
            else:
                tsv_path.write_text(tsv_chunk, encoding="utf-8")

    # Run full statistical & accuracy evaluation across all 100,000 domains
    eval_summary = evaluate_100k_jsonl(jsonl_path)
    out_summary_path = examples_dir / "benchmark_100000_summary.json"
    out_summary_path.write_text(json.dumps(eval_summary, indent=2), encoding="utf-8")

    print("\n" + "=" * 100, flush=True)
    print("100,000 TED DOMAINS STATISTICAL & ACCURACY EVALUATION SUMMARY", flush=True)
    print("=" * 100, flush=True)
    print(json.dumps(eval_summary, indent=2), flush=True)


if __name__ == "__main__":
    asyncio.run(main())
