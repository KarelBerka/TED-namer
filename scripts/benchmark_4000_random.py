"""Sample and test 4,000 new random TED domains (4x 1,000 batches) and compare distributions."""

from __future__ import annotations

import asyncio
import json
import random
import time
from collections import Counter
from pathlib import Path

import httpx

from ted_namer.cache import ApiCache
from ted_namer.clients.ted import TedClient
from ted_namer.export import save_output
from ted_namer.models import ClusterNamingResult, TedDomainMember
from ted_namer.pipeline import TedNamingPipeline

# Diverse taxonomic buckets covering all major kingdoms of life
_UNIPROT_SAMPLE_QUERIES = [
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606)",       # Human
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090)",      # Mouse
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7955)",       # Zebrafish
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7227)",       # Drosophila
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:6239)",       # C. elegans
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83333)",      # E. coli K-12
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:224308)",     # B. subtilis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:559292)",     # S. cerevisiae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4896)",       # S. pombe
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702)",       # A. thaliana
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4577)",       # Zea mays (Maize)
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83332)",      # M. tuberculosis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2157)",       # Archaea reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2)",          # Bacteria reviewed
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2)",         # Bacteria TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2759)",      # Eukaryota TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2157)",      # Archaea TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:10239)",     # Viruses
]


async def fetch_4000_new_ted_domains(
    cache: ApiCache,
    existing_ted_ids: set[str],
    target_count: int = 4000,
    seed: int = 777,
) -> list[TedDomainMember]:
    """Fetch 4,000 new random TED domains distinct from previously sampled ones."""
    rng = random.Random(seed)
    headers = {"User-Agent": "ted-namer-benchmark-4000/0.1.0"}
    candidate_uniprots: list[str] = []

    async with httpx.AsyncClient(headers=headers, follow_redirects=True, timeout=30.0) as http:
        for q in _UNIPROT_SAMPLE_QUERIES:
            url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": q,
                "format": "json",
                "size": "400",
                "fields": "accession",
            }
            try:
                resp = await http.get(url, params=params)
                if resp.status_code == 200:
                    results = resp.json().get("results") or []
                    accs = [r["primaryAccession"] for r in results if r.get("primaryAccession")]
                    rng.shuffle(accs)
                    candidate_uniprots.extend(accs[:250])
            except Exception as exc:
                print(f"Warning fetching UniProt query '{q}': {exc}")

        # Deduplicate and shuffle candidate accessions
        seen_u: set[str] = set()
        unique_uniprots: list[str] = []
        for u in candidate_uniprots:
            if u not in seen_u:
                seen_u.add(u)
                unique_uniprots.append(u)
        rng.shuffle(unique_uniprots)
        print(f"Collected {len(unique_uniprots)} candidate UniProt accessions across {len(_UNIPROT_SAMPLE_QUERIES)} buckets.")

        ted_client = TedClient(http_client=http, cache=cache)
        sem = asyncio.Semaphore(16)

        async def _get_ted(acc: str) -> list[TedDomainMember]:
            async with sem:
                return await ted_client.get_domains_for_uniprot(acc)

        sampled_domains: list[TedDomainMember] = []
        seen_new_ids: set[str] = set(existing_ted_ids)
        batch_size = 40

        for i in range(0, len(unique_uniprots), batch_size):
            batch = unique_uniprots[i : i + batch_size]
            batch_res = await asyncio.gather(*(_get_ted(a) for a in batch))
            for dom_list in batch_res:
                if not dom_list:
                    continue
                k = min(len(dom_list), rng.choice([1, 1, 2, 2, 3]))
                chosen = rng.sample(dom_list, k=k)
                for d in chosen:
                    if d.ted_id not in seen_new_ids:
                        seen_new_ids.add(d.ted_id)
                        sampled_domains.append(d)
            if len(sampled_domains) >= target_count:
                break
            if len(sampled_domains) % 400 < 50 or i % 160 == 0:
                print(f"  Sampled {len(sampled_domains)}/{target_count} new unique TED domains so far...")

    return sampled_domains[:target_count]


def compute_batch_stats(results: list[ClusterNamingResult]) -> dict[str, any]:
    """Compute summary metrics for a list of named results."""
    n = len(results)
    if n == 0:
        return {}
    avg_conf = sum(r.confidence_score for r in results) / n
    tiers = Counter(r.specificity_tier.value for r in results)
    sources = Counter(
        r.members[0].domain.cath_estimation_source if r.members else "-" for r in results
    )
    conf_below_50 = sum(1 for r in results if r.confidence_score < 0.50)
    low_plddt = sum(1 for r in results if r.members and r.members[0].domain.plddt < 60)

    # Group CATH sources into high-level categories
    official_cath = sources.get("ted_official", 0)
    t_level_cath = sources.get("ted_t_plus_sse", 0) + sources.get("ted_t_plus_g3dsa", 0)
    g3dsa_cath = sources.get("g3dsa_overlap", 0)
    direct_3d_sse_rules = sum(v for k, v in sources.items() if k.startswith("rule_"))

    return {
        "count": n,
        "mean_confidence": avg_conf,
        "tier_specific_domain": tiers.get("specific_domain", 0) / n,
        "tier_family_context": tiers.get("protein_family_context", 0) / n,
        "tier_superfamily": tiers.get("superfamily", 0) / n,
        "tier_topology_fold": tiers.get("topology_fold", 0) / n,
        "cath_official": official_cath / n,
        "cath_t_level": t_level_cath / n,
        "cath_direct_3d_rules": direct_3d_sse_rules / n,
        "cath_g3dsa": g3dsa_cath / n,
        "conf_below_50_count": conf_below_50,
        "conf_below_50_pct": conf_below_50 / n,
        "low_plddt_count": low_plddt,
        "low_plddt_pct": low_plddt / n,
    }


async def main() -> None:
    examples_dir = Path("examples")
    examples_dir.mkdir(exist_ok=True)
    cache = ApiCache()

    # Load previously sampled 1,000 domains (Batch 1) to avoid overlap
    existing_ted_ids: set[str] = set()
    b1_results_path = examples_dir / "random_1000_results.json"
    batch_1_results: list[ClusterNamingResult] = []
    if b1_results_path.exists():
        with open(b1_results_path, encoding="utf-8") as f:
            b1_data = json.load(f)
            batch_1_results = [ClusterNamingResult.model_validate(x) for x in b1_data]
            for r in batch_1_results:
                for m in r.members:
                    existing_ted_ids.add(m.domain.ted_id)
        print(f"Loaded Batch 1: {len(batch_1_results)} existing TED domains.")

    t0 = time.perf_counter()
    print("\nSampling 4,000 NEW random TED domains across diverse taxa (distinct from Batch 1)...")
    new_domains = await fetch_4000_new_ted_domains(cache=cache, existing_ted_ids=existing_ted_ids, target_count=4000, seed=999)
    print(f"Successfully collected {len(new_domains)} new random TED domains in {time.perf_counter() - t0:.1f}s.")

    # Split into 4 distinct batches of 1,000 domains (Batches 2, 3, 4, 5)
    batches = [
        ("Batch 2 (1,001-2,000)", new_domains[0:1000]),
        ("Batch 3 (2,001-3,000)", new_domains[1000:2000]),
        ("Batch 4 (3,001-4,000)", new_domains[2000:3000]),
        ("Batch 5 (4,001-5,000)", new_domains[3000:4000]),
    ]

    pipeline = TedNamingPipeline(use_cache=True, concurrency=20, use_llm=False)
    all_new_results: list[ClusterNamingResult] = []
    batch_stats_list: list[tuple[str, dict[str, any]]] = []

    if batch_1_results:
        batch_stats_list.append(("Batch 1 (1-1,000)", compute_batch_stats(batch_1_results)))

    # Process each 1,000-domain batch sequentially and record statistics
    for b_name, b_domains in batches:
        print(f"\n=======================================================")
        print(f"Running pipeline on {b_name} ({len(b_domains)} domains)...")
        print(f"=======================================================")
        cmap = {
            f"Q_{dom.uniprot_acc}_{dom.short_ted_num}": [dom.ted_id]
            for dom in b_domains
        }

        b_start = time.perf_counter()
        last_log = 0

        def _prog(done: int, total: int) -> None:
            nonlocal last_log
            if done - last_log >= 100 or done == total:
                last_log = done
                el = time.perf_counter() - b_start
                rate = done / max(0.1, el)
                print(f"  [{b_name}] {done:4d}/{total:4d} ({done/total*100:5.1f}%) | {el:5.1f}s | {rate:4.1f} dom/s")

        b_res = await pipeline.name_clusters(cluster_map=cmap, group_by="cluster_id", progress_callback=_prog)
        b_elapsed = time.perf_counter() - b_start
        print(f"  Finished {b_name} in {b_elapsed:.1f}s ({b_elapsed/max(1, len(b_res)):.3f}s/query).")

        all_new_results.extend(b_res)
        b_stats = compute_batch_stats(b_res)
        batch_stats_list.append((b_name, b_stats))

        # Checkpoint incremental saves
        save_output(all_new_results, examples_dir / "random_4000_results.json")
        save_output(all_new_results, examples_dir / "random_4000_results.tsv")

    # Combine Batch 1 (1,000) + Batches 2-5 (4,000) = 5,000 Total Domains!
    total_5000_results = batch_1_results + all_new_results
    total_5000_stats = compute_batch_stats(total_5000_results)

    print("\nSaving complete 5,000 and 4,000 domain outputs...")
    save_output(total_5000_results, examples_dir / "random_5000_results.json")
    save_output(total_5000_results, examples_dir / "random_5000_results.tsv")
    save_output(total_5000_results, examples_dir / "random_5000_report.md")
    save_output(all_new_results, examples_dir / "random_4000_report.md")

    # Print comprehensive comparison table across all 5 batches and the 5,000-domain total
    print("\n" + "=" * 95)
    print("STATISTICAL COMPARISON ACROSS 5 INDEPENDENT BATCHES (5,000 TED DOMAINS TOTAL)")
    print("=" * 95)
    header = (
        f"{'Metric':<32} | "
        + " | ".join(f"{b[0][:7]:^9}" for b in batch_stats_list)
        + f" | {'TOTAL 5K':^9}"
    )
    print(header)
    print("-" * len(header))

    metrics_to_show = [
        ("Domain Count", "count", "{:^9d}"),
        ("Mean Confidence", "mean_confidence", "{:^9.3f}"),
        ("Tier 1: Specific Domain", "tier_specific_domain", "{:^9.1%}"),
        ("Tier 2: Superfamily", "tier_superfamily", "{:^9.1%}"),
        ("Tier 3: Family Context", "tier_family_context", "{:^9.1%}"),
        ("Tier 4: Topology Fold", "tier_topology_fold", "{:^9.1%}"),
        ("CATH: Official (4-level)", "cath_official", "{:^9.1%}"),
        ("CATH: TED T-level (.T10)", "cath_t_level", "{:^9.1%}"),
        ("CATH: Direct 3D/SSE Rules", "cath_direct_3d_rules", "{:^9.1%}"),
        ("CATH: Gene3D Overlap", "cath_g3dsa", "{:^9.1%}"),
        ("Hard Cases (Conf < 0.50)", "conf_below_50_pct", "{:^9.1%}"),
        ("Low pLDDT (< 60)", "low_plddt_pct", "{:^9.1%}"),
    ]

    for label, key, fmt in metrics_to_show:
        row_vals = []
        for _, s in batch_stats_list:
            val = s.get(key, 0)
            row_vals.append(fmt.format(val))
        tot_val = total_5000_stats.get(key, 0)
        tot_str = fmt.format(tot_val)
        print(f"{label:<32} | " + " | ".join(row_vals) + f" | {tot_str}")

    print("=" * 95)
    print("Benchmark complete! All tables and reports successfully generated.")


if __name__ == "__main__":
    asyncio.run(main())
