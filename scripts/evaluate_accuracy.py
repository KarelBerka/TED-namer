"""Comprehensive blind (hold-out & 5-fold CV) accuracy evaluation for TED-namer on 10,000 domains."""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from ted_namer.aggregator import aggregate_cluster_evidence, normalize_domain_name
from ted_namer.cache import ApiCache
from ted_namer.cath_classifier import (
    _G3DSA_RE,
    classify_cath_with_rules,
    format_t_prefixed_cath,
    infer_cath_classification,
)
from ted_namer.export import save_output
from ted_namer.models import ClusterNamingResult, EnrichedTedDomain, SpecificityTier, TedDomainMember


def evaluate_cath_holdout(results: list[ClusterNamingResult], cache: ApiCache) -> dict[str, any]:
    """Run blind hold-out and 5-fold cross-validation on all official 4-level CATH domains."""
    official_entries: list[tuple[list[str], dict, TedDomainMember]] = []
    for r in results:
        for m in r.members:
            dom = m.domain
            raw = (dom.cath_label or "-").split(",")[0].strip()
            parts = [p.strip() for p in raw.split(".") if p.strip() and p.strip() != "-"]
            if len(parts) >= 4:
                official_entries.append((parts[:4], m.model_dump(), dom.model_copy(deep=True)))

    n_total = len(official_entries)
    print(f" Evaluating CATH hold-out accuracy on {n_total} domains with official 4-level TED CATH...")

    # Pre-build 5-fold out-of-sample profile->CATH maps so there is zero data leakage
    num_folds = 5
    fold_maps: list[dict[str, tuple[list[str], int]]] = []
    for fold_idx in range(num_folds):
        acc_to_cath: dict[str, Counter[str]] = defaultdict(Counter)
        for i, (tp, m_dict, _dom) in enumerate(official_entries):
            if i % num_folds == fold_idx:
                continue  # held-out test fold
            cath_str = ".".join(tp)
            for a in m_dict.get("overlapping_annotations", []):
                if a.get("source_db") == "cathgene3d" or str(a.get("accession", "")).upper().startswith("G3DSA:"):
                    continue
                if a.get("iou", 0) >= 0.40 or (a.get("ted_coverage", 0) >= 0.65 and a.get("annot_coverage", 0) >= 0.45):
                    for key in (a.get("accession"), a.get("integrated_ipr")):
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

    # Metrics for Cascade A (TED CATH masked; G3DSA -> 5-fold CV Profile-to-CATH -> 3D/SSE rules)
    cas_a_c, cas_a_ca, cas_a_cat, cas_a_cath = 0, 0, 0, 0
    cas_a_by_src = Counter()
    cas_a_src_cat_ok = Counter()
    cas_a_direct_claimed = Counter()
    cas_a_direct_ok = Counter()

    # Metrics for Cascade B (Both TED CATH and G3DSA masked; 5-fold CV Profile-to-CATH -> 3D/SSE rules)
    cas_b_c, cas_b_ca, cas_b_cat, cas_b_cath = 0, 0, 0, 0
    cas_b_prof_hits = 0
    cas_b_prof_c, cas_b_prof_ca, cas_b_prof_cat, cas_b_prof_cath = 0, 0, 0, 0

    # Metrics for Cascade C (Pure 3D/SSE rules only, zero sequence/profile annotations)
    cas_c_c, cas_c_ca, cas_c_cat = 0, 0, 0
    cas_c_pure_alpha_beta_tot, cas_c_pure_alpha_beta_ok = 0, 0

    for i, (tp, m_dict, dom) in enumerate(official_entries):
        fmap = fold_maps[i % num_folds]
        pdb_text = cache.get("ted_pdb", dom.ted_id)
        annots = m_dict.get("overlapping_annotations", [])

        # --- Check G3DSA ---
        best_g3d = None
        best_g3d_iou = 0.0
        best_g3d_tcov = 0.0
        best_g3d_sc = -1.0
        for a in annots:
            if a.get("source_db") == "cathgene3d" or str(a.get("accession", "")).upper().startswith("G3DSA:"):
                if a.get("iou", 0) >= 0.25 or a.get("ted_coverage", 0) >= 0.45:
                    m_re = _G3DSA_RE.match(str(a.get("accession", "")).strip())
                    if m_re:
                        sc = a.get("iou", 0) + 0.5 * a.get("ted_coverage", 0)
                        if sc > best_g3d_sc:
                            best_g3d_sc = sc
                            best_g3d_iou = a.get("iou", 0)
                            best_g3d_tcov = a.get("ted_coverage", 0)
                            best_g3d = [m_re.group(1), m_re.group(2), m_re.group(3), m_re.group(4)]

        # --- Check 5-fold CV Profile-to-CATH ---
        best_prof = None
        best_prof_depth = 1
        best_prof_iou = 0.0
        best_prof_sc = -1.0
        for a in annots:
            if a.get("source_db") == "cathgene3d" or str(a.get("accession", "")).upper().startswith("G3DSA:"):
                continue
            if a.get("iou", 0) >= 0.40 or (a.get("ted_coverage", 0) >= 0.65 and a.get("annot_coverage", 0) >= 0.45):
                for key in (a.get("accession"), a.get("integrated_ipr")):
                    if key and str(key) in fmap:
                        p_parts, p_depth = fmap[str(key)]
                        db_boost = 1.15 if a.get("source_db") == "ssf" else 1.0
                        sc = (a.get("iou", 0) + 0.3 * a.get("ted_coverage", 0)) * db_boost * (1.1 if p_depth >= 3 else 1.0)
                        if sc > best_prof_sc:
                            best_prof_sc = sc
                            best_prof = p_parts
                            best_prof_depth = p_depth
                            best_prof_iou = a.get("iou", 0)

        # --- Check Pure 3D/SSE Rules ---
        rule_parts, rule_depth, _rule_name = classify_cath_with_rules(dom, pdb_text=pdb_text)
        if rule_parts[0] == tp[0]:
            cas_c_c += 1
        if rule_parts[:2] == tp[:2]:
            cas_c_ca += 1
        if rule_parts[:3] == tp[:3]:
            cas_c_cat += 1
        if dom.num_strand == 0 or dom.num_helix == 0:
            cas_c_pure_alpha_beta_tot += 1
            if rule_parts[0] == tp[0]:
                cas_c_pure_alpha_beta_ok += 1

        # --- Cascade B (No TED CATH, No G3DSA) ---
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

        # --- Cascade A (Full fallback cascade when TED CATH is masked) ---
        if best_g3d:
            pred_a = best_g3d
            src_a = "g3dsa_overlap"
            k_depth = 4 if (best_g3d_iou >= 0.65 and best_g3d_tcov >= 0.75) else (3 if (best_g3d_iou >= 0.40 or best_g3d_tcov >= 0.65) else 1)
        elif best_prof:
            pred_a = best_prof
            src_a = "profile_cath_map"
            k_depth = best_prof_depth if best_prof_iou >= 0.45 else 1
        else:
            pred_a = rule_parts
            src_a = "rule_3d_sse"
            k_depth = rule_depth

        cas_a_by_src[src_a] += 1
        if pred_a[0] == tp[0]:
            cas_a_c += 1
        if pred_a[:2] == tp[:2]:
            cas_a_ca += 1
        if pred_a[:3] == tp[:3]:
            cas_a_cat += 1
            cas_a_src_cat_ok[src_a] += 1
        if pred_a[:4] == tp[:4]:
            cas_a_cath += 1

        for d_lvl in range(1, k_depth + 1):
            cas_a_direct_claimed[d_lvl] += 1
            if pred_a[:d_lvl] == tp[:d_lvl]:
                cas_a_direct_ok[d_lvl] += 1

    return {
        "n_official": n_total,
        "cascade_a": {
            "c_acc": cas_a_c / n_total,
            "ca_acc": cas_a_ca / n_total,
            "cat_acc": cas_a_cat / n_total,
            "cath_acc": cas_a_cath / n_total,
            "by_src": dict(cas_a_by_src),
            "src_cat_acc": {k: cas_a_src_cat_ok[k] / max(1, v) for k, v in cas_a_by_src.items()},
            "direct_level_precision": {
                d: cas_a_direct_ok[d] / max(1, cas_a_direct_claimed[d])
                for d in sorted(cas_a_direct_claimed)
            },
            "direct_level_counts": dict(cas_a_direct_claimed),
        },
        "cascade_b_no_g3dsa": {
            "c_acc": cas_b_c / n_total,
            "ca_acc": cas_b_ca / n_total,
            "cat_acc": cas_b_cat / n_total,
            "cath_acc": cas_b_cath / n_total,
            "profile_map_coverage": cas_b_prof_hits / n_total,
            "profile_map_c_acc": cas_b_prof_c / max(1, cas_b_prof_hits),
            "profile_map_ca_acc": cas_b_prof_ca / max(1, cas_b_prof_hits),
            "profile_map_cat_acc": cas_b_prof_cat / max(1, cas_b_prof_hits),
            "profile_map_cath_acc": cas_b_prof_cath / max(1, cas_b_prof_hits),
        },
        "cascade_c_pure_3d_sse": {
            "c_acc": cas_c_c / n_total,
            "pure_alpha_or_beta_c_acc": cas_c_pure_alpha_beta_ok / max(1, cas_c_pure_alpha_beta_tot),
            "ca_acc": cas_c_ca / n_total,
            "cat_acc": cas_c_cat / n_total,
        },
    }


def evaluate_naming_concordance_and_refresh(results: list[ClusterNamingResult]) -> dict[str, any]:
    """Evaluate leave-one-database-out naming concordance and refresh all 10,000 results with the improved aggregator."""
    # 1. Leave-Pfam-Out Cross-Validation on domains that have both Pfam and another Tier 1/2 DB (CDD/SMART/PROSITE/InterPro)
    lodo_tot = 0
    lodo_exact_ipr_or_norm = 0
    lodo_token_overlap = 0

    # 2. Refresh all 10,000 results using the updated CATH classifier and aggregator
    refreshed_results: list[ClusterNamingResult] = []
    legacy_pdb_before = 0
    legacy_pdb_after = 0
    sublobe_after = 0

    pdb_noise_re = re.compile(r"(;\s*domain\s+\d+|\.\s*Chain\s+[A-Z0-9]+|Vaccinia Virus protein VP39)", re.IGNORECASE)

    for r in results:
        if pdb_noise_re.search(r.consensus_name):
            legacy_pdb_before += 1

        # Clone members and re-run CATH inference + cluster aggregation
        new_members: list[EnrichedTedDomain] = []
        for m in r.members:
            mc = m.model_copy(deep=True)
            # Also fix legacy G3DSA superfamily_name on cath_info if integrated_ipr_name is available on annotations
            if mc.cath_info and mc.cath_info.integrated_ipr:
                for a in mc.overlapping_annotations:
                    if a.integrated_ipr == mc.cath_info.integrated_ipr and a.integrated_ipr_name:
                        mc.cath_info.superfamily_name = a.integrated_ipr_name
                        break
            inferred_code, _num_code, est_src = infer_cath_classification(
                member=mc.domain,
                overlapping_annotations=mc.overlapping_annotations,
            )
            mc.domain.inferred_cath_label = inferred_code
            mc.domain.cath_estimation_source = est_src
            new_members.append(mc)

        new_res = aggregate_cluster_evidence(r.cluster_id, new_members)
        refreshed_results.append(new_res)

        if pdb_noise_re.search(new_res.consensus_name):
            legacy_pdb_after += 1
        if "subdomain" in new_res.consensus_name.lower():
            sublobe_after += 1

        # Check Leave-Pfam-Out on domains where Pfam was a supporting database of a Tier 1 domain
        if new_res.specificity_tier == SpecificityTier.SPECIFIC_DOMAIN and new_res.top_candidates:
            top_c = new_res.top_candidates[0]
            if "pfam" in top_c.source_databases and len(new_members) == 1:
                m0 = new_members[0]
                non_pfam = [
                    a for a in m0.overlapping_annotations
                    if a.source_db not in {"pfam"} and a.annotation_type == "domain" and a.iou >= 0.35
                ]
                if non_pfam:
                    lodo_tot += 1
                    m_holdout = m0.model_copy(deep=True)
                    m_holdout.overlapping_annotations = [
                        a for a in m0.overlapping_annotations if a.source_db != "pfam"
                    ]
                    res_holdout = aggregate_cluster_evidence(r.cluster_id, [m_holdout])
                    norm_full = normalize_domain_name(new_res.consensus_name)
                    norm_ho = normalize_domain_name(res_holdout.consensus_name)
                    iprs_full = {a for a in top_c.supporting_accessions if a.startswith("IPR")}
                    iprs_ho = (
                        {a for a in res_holdout.top_candidates[0].supporting_accessions if a.startswith("IPR")}
                        if res_holdout.top_candidates
                        else set()
                    )
                    if norm_full == norm_ho or (iprs_full and (iprs_full & iprs_ho)):
                        lodo_exact_ipr_or_norm += 1
                        lodo_token_overlap += 1
                    else:
                        t1 = {w for w in norm_full.split() if len(w) >= 3}
                        t2 = {w for w in norm_ho.split() if len(w) >= 3}
                        if t1 and t2 and (t1 & t2):
                            lodo_token_overlap += 1

    # Compute refreshed 10,000 summary stats
    n = len(refreshed_results)
    tiers = Counter(r.specificity_tier.value for r in refreshed_results)
    sources = Counter(r.members[0].domain.cath_estimation_source for r in refreshed_results if r.members)
    mean_conf = sum(r.confidence_score for r in refreshed_results) / n
    tier_confs = defaultdict(list)
    for r in refreshed_results:
        tier_confs[r.specificity_tier.value].append(r.confidence_score)

    return {
        "refreshed_results": refreshed_results,
        "lodo_n": lodo_tot,
        "lodo_exact_concordance": lodo_exact_ipr_or_norm / max(1, lodo_tot),
        "lodo_semantic_concordance": lodo_token_overlap / max(1, lodo_tot),
        "legacy_pdb_noise_before": legacy_pdb_before,
        "legacy_pdb_noise_after": legacy_pdb_after,
        "sublobe_domains_resolved": sublobe_after,
        "mean_confidence": mean_conf,
        "mean_conf_by_tier": {k: sum(v) / len(v) for k, v in tier_confs.items()},
        "tiers": {k: v / n for k, v in tiers.items()},
        "cath_sources": {k: v / n for k, v in sources.items()},
        "conf_below_50_pct": sum(1 for r in refreshed_results if r.confidence_score < 0.50) / n,
    }


def main() -> None:
    cache = ApiCache()
    results_10k_path = Path("examples/random_10000_results.json")
    results_20k_path = Path("examples/random_20000_results.json")
    target_path = results_20k_path if results_20k_path.exists() else results_10k_path

    print(f"Loading {target_path}...")
    raw_data = json.loads(target_path.read_text(encoding="utf-8"))
    results = [ClusterNamingResult.model_validate(x) for x in raw_data]

    cath_eval = evaluate_cath_holdout(results, cache)
    naming_eval = evaluate_naming_concordance_and_refresh(results)

    refreshed = naming_eval.pop("refreshed_results")
    print(f"\nSaving refreshed {len(refreshed)} domain results...")
    save_output(refreshed, target_path)
    save_output(refreshed, target_path.with_suffix(".tsv"))
    if len(refreshed) >= 20000:
        save_output(refreshed[:10000], Path("examples/random_10000_results.json"))
        save_output(refreshed[:10000], Path("examples/random_10000_results.tsv"))
        save_output(refreshed[10000:20000], Path("examples/random_10000_batches_11_to_20.json"))
        save_output(refreshed[10000:20000], Path("examples/random_10000_batches_11_to_20.tsv"))

    summary = {"cath_holdout_validation": cath_eval, "naming_validation_and_refresh": naming_eval}
    Path("examples/accuracy_validation_20000_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print("\n" + "=" * 90)
    print(f"ACCURACY & CROSS-VALIDATION SUMMARY ({len(refreshed)} TED DOMAINS)")
    print("=" * 90)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
