"""Offline re-evaluation of all 100,000 TED domains using the upgraded naming hygiene,
InterPro annotation_type rules, curated unintegrated CATH H-names, '?' wildcard notation,
and structural quality flags — computing Before vs. After validity statistics and 5-fold CV."""

from __future__ import annotations

import csv
import io
import json
import re
import time
from collections import Counter
from pathlib import Path
from typing import Any

from ted_namer.aggregator import (
    _LEGACY_PDB_HEADER_RE,
    aggregate_cluster_evidence,
    is_bare_positional_name,
    is_low_information_name,
    is_sentence_like_description,
)
from ted_namer.cache import ApiCache
from ted_namer.cath_classifier import (
    CURATED_UNINTEGRATED_CATH_H_NAMES,
    infer_cath_classification,
    is_estimated_cath_code,
    is_novel_fold_source,
    strip_t_prefixes,
)
from ted_namer.export import to_delimited
from ted_namer.models import (
    CathResolution,
    ClusterNamingResult,
    EnrichedTedDomain,
    OverlapCategory,
    OverlappingAnnotation,
    ProteinContext,
    TedDomainMember,
)
from scripts.benchmark_100000_random import evaluate_100k_jsonl, result_to_compact

_BRACKET_RE = re.compile(r"\[[^\]]+\]")
_DOUBLE_DOM_RE = re.compile(r"\b(domain\s+domain|domain-like\s+domain)\b", re.IGNORECASE)
_ANY_DOUBLE_DOM_RE = re.compile(r"\bdomain\b.*\bdomain\b", re.IGNORECASE)


def audit_names(rows: list[dict[str, Any]], label_key: str = "consensus_name", cath_key: str = "dominant_cath_label") -> dict[str, Any]:
    """Compute validity & hygiene audit counts across 100k records."""
    n = len(rows)
    legacy_pdb = 0
    sentence_like = 0
    bare_pos = 0
    brackets = 0
    low_info = 0
    adjacent_double_dom = 0
    any_double_dom = 0
    gt90 = 0
    t_prefix = 0
    wildcard_q = 0
    unknown_domain = 0
    residual_legacy = Counter()

    for r in rows:
        nm = r.get(label_key) or ""
        cath = r.get(cath_key) or ""
        if _LEGACY_PDB_HEADER_RE.search(nm):
            legacy_pdb += 1
            residual_legacy[nm] += 1
        if is_sentence_like_description(nm):
            sentence_like += 1
        if is_bare_positional_name(nm):
            bare_pos += 1
        if _BRACKET_RE.search(nm):
            brackets += 1
        if is_low_information_name(nm):
            low_info += 1
        if "unknown" in nm.lower():
            unknown_domain += 1
        if _DOUBLE_DOM_RE.search(nm):
            adjacent_double_dom += 1
        if _ANY_DOUBLE_DOM_RE.search(nm):
            any_double_dom += 1
        if len(nm) > 90:
            gt90 += 1
        if any(p.strip().upper().startswith("T") for p in cath.split(".") if p.strip()):
            t_prefix += 1
        if "?" in cath:
            wildcard_q += 1

    return {
        "total": n,
        "legacy_pdb_headers": legacy_pdb,
        "sentence_like_descriptions": sentence_like,
        "bare_positional_names": bare_pos,
        "bracketed_gene_or_ec_tags": brackets,
        "low_info_or_unknown": low_info,
        "unknown_keyword": unknown_domain,
        "adjacent_double_domain": adjacent_double_dom,
        "any_double_domain": any_double_dom,
        "length_gt_90": gt90,
        "legacy_t_prefix_cath": t_prefix,
        "wildcard_q_cath": wildcard_q,
        "top_residual_legacy_names": residual_legacy.most_common(10),
    }


def _reconstruct_enriched(rec: dict[str, Any], pdb_cache: ApiCache) -> EnrichedTedDomain:
    dom = TedDomainMember(
        ted_id=rec["ted_id"],
        uniprot_acc=rec["uniprot_acc"],
        chopping=rec["chopping"],
        nres_domain=rec["nres_domain"],
        num_helix=rec.get("num_helix", 0),
        num_strand=rec.get("num_strand", 0),
        num_turn=rec.get("num_turn", 0),
        norm_rg=rec.get("norm_rg"),
        packing_density=rec.get("packing_density"),
        plddt=rec.get("plddt", 0.0),
        cath_label=rec.get("cath_label", "-"),
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
        for a in rec.get("annots", [])
    ]

    pdb_text = pdb_cache.get("ted_domain_pdb", dom.ted_id)
    inferred_cath, numeric_hint, est_source = infer_cath_classification(
        member=dom,
        overlapping_annotations=ov_list,
        pdb_text=pdb_text,
    )
    dom.inferred_cath_label = inferred_cath
    dom.cath_estimation_source = est_source
    dom.cath_fold_hint = numeric_hint if is_estimated_cath_code(inferred_cath) else None
    dom.novel_fold_candidate = is_novel_fold_source(est_source)

    confirmed_parts = [p for p in inferred_cath.split(".") if p and p != "?" and not p.upper().startswith("T")]
    confirmed_depth = len(confirmed_parts)

    cath_info = None
    if inferred_cath != "-":
        topo_code = ".".join(confirmed_parts[:3]) if confirmed_depth >= 3 else None
        topo_name = None
        if topo_code and strip_t_prefixes(rec.get("cath_t_code") or "") == topo_code:
            topo_name = rec.get("cath_t_name")

        can_have_h = (confirmed_depth >= 4) or (
            est_source in {"g3dsa_overlap", "ted_t_plus_g3dsa", "ted_t_plus_g3dsa_h", "profile_cath_map", "ted_t_plus_profile"}
            and confirmed_depth >= 3
        )
        sf_name = rec.get("cath_h_name") if can_have_h else None
        sf_ipr = rec.get("cath_h_ipr") if can_have_h else None
        if can_have_h and not sf_ipr and numeric_hint in CURATED_UNINTEGRATED_CATH_H_NAMES:
            sf_name = CURATED_UNINTEGRATED_CATH_H_NAMES[numeric_hint]

        cath_info = CathResolution(
            cath_label=inferred_cath,
            numeric_cath_code=strip_t_prefixes(inferred_cath),
            raw_ted_cath_label=dom.cath_label,
            is_estimated=is_estimated_cath_code(inferred_cath) or est_source != "ted_official",
            estimation_source=est_source,
            assignment_level="H" if confirmed_depth >= 4 and est_source == "ted_official" else "inferred",
            superfamily_name=sf_name,
            integrated_ipr=sf_ipr,
            topology_code=topo_code,
            topology_name=topo_name,
            fold_hint=dom.cath_fold_hint,
            novel_fold_candidate=dom.novel_fold_candidate,
        )

    p_ctx = ProteinContext(
        uniprot_acc=rec["uniprot_acc"],
        recommended_name=rec.get("recommended_name"),
        is_reviewed=rec.get("is_reviewed", False),
        protein_length=rec.get("protein_length", 0),
        similarity_families=rec.get("similarity_families", []),
    )

    return EnrichedTedDomain(
        domain=dom,
        cath_info=cath_info,
        protein_context=p_ctx,
        overlapping_annotations=ov_list,
        relative_position=rec.get("relative_position", "domain"),
    )


def main() -> None:
    t0 = time.perf_counter()
    examples_dir = Path("examples")
    jsonl_path = examples_dir / "random_100000_compact.jsonl"
    tmp_jsonl = examples_dir / "random_100000_compact.tmp.jsonl"
    tsv_path = examples_dir / "random_100000_results.tsv"

    cache = ApiCache()
    before_audit_rows: list[dict[str, Any]] = []
    after_compact_rows: list[dict[str, Any]] = []
    quality_flag_counts: Counter[str] = Counter()
    novel_fold_count = 0

    print(f"Streaming and re-evaluating {jsonl_path}...", flush=True)
    header_written = False

    with (
        open(jsonl_path, encoding="utf-8") as fin,
        open(tmp_jsonl, "w", encoding="utf-8") as jout,
        open(tsv_path, "w", encoding="utf-8", newline="\n") as tout,
    ):
        batch_results: list[ClusterNamingResult] = []
        for idx, line in enumerate(fin, start=1):
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            before_audit_rows.append(
                {
                    "consensus_name": rec.get("consensus_name", ""),
                    "dominant_cath_label": rec.get("dominant_cath_label", ""),
                }
            )

            em = _reconstruct_enriched(rec, cache)
            res = aggregate_cluster_evidence(rec["cluster_id"], [em])
            if res.novel_fold_candidate:
                novel_fold_count += 1
            for qf in res.structure_quality_flags:
                quality_flag_counts[qf] += 1

            compact = result_to_compact(res)
            jout.write(json.dumps(compact) + "\n")
            after_compact_rows.append(
                {
                    "consensus_name": compact["consensus_name"],
                    "dominant_cath_label": compact["dominant_cath_label"],
                }
            )
            batch_results.append(res)

            if len(batch_results) >= 5000:
                tsv_chunk = to_delimited(batch_results, delimiter="\t")
                lines = tsv_chunk.splitlines(keepends=True)
                if not header_written:
                    tout.writelines(lines)
                    header_written = True
                else:
                    tout.writelines(lines[1:])
                batch_results.clear()

            if idx % 20000 == 0:
                print(f"  Re-evaluated {idx:,}/100,000 domains in {time.perf_counter() - t0:.1f}s...", flush=True)

        if batch_results:
            tsv_chunk = to_delimited(batch_results, delimiter="\t")
            lines = tsv_chunk.splitlines(keepends=True)
            if not header_written:
                tout.writelines(lines)
            else:
                tout.writelines(lines[1:])

    tmp_jsonl.replace(jsonl_path)
    print(f"Finished refreshing TSV & compact JSONL in {time.perf_counter() - t0:.1f}s.", flush=True)

    before_audit = audit_names(before_audit_rows)
    after_audit = audit_names(after_compact_rows)

    eval_summary = evaluate_100k_jsonl(jsonl_path)
    eval_summary["validity_audit"] = {
        "before": before_audit,
        "after": after_audit,
        "novel_fold_candidate_count": novel_fold_count,
        "novel_fold_candidate_pct": round(novel_fold_count / max(1, len(after_compact_rows)), 5),
        "structure_quality_flag_counts": dict(quality_flag_counts),
    }

    summary_path = examples_dir / "benchmark_100000_summary.json"
    summary_path.write_text(json.dumps(eval_summary, indent=2), encoding="utf-8")

    print("\n" + "=" * 95)
    print("BEFORE VS AFTER VALIDITY HYGIENE AUDIT (100,000 TED DOMAINS)")
    print("=" * 95)
    print(json.dumps(eval_summary["validity_audit"], indent=2))
    print("\n" + "=" * 95)
    print("REFRESHED 100,000-DOMAIN ACCURACY & DECILES SUMMARY")
    print("=" * 95)
    print(json.dumps({k: v for k, v in eval_summary.items() if k != "deciles_10k"}, indent=2))


if __name__ == "__main__":
    main()
