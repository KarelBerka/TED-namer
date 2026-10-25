"""Sample and test 5,000 new random TED domains (5x 1,000 batches: 6-10) and compare across all 10,000 domains."""

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

# Diverse taxonomic buckets covering all major kingdoms of life, model organisms, and unreviewed sets
_UNIPROT_SAMPLE_QUERIES = [
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606)",       # Human
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090)",      # Mouse
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10116)",      # Rat
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9031)",       # Chicken (Gallus)
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7955)",       # Zebrafish
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7227)",       # Drosophila
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:6239)",       # C. elegans
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83333)",      # E. coli K-12
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:224308)",     # B. subtilis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:1883)",       # Streptomyces
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:1117)",       # Cyanobacteria
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:976)",        # Bacteroidetes
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:559292)",     # S. cerevisiae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4896)",       # S. pombe
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4751)",       # Fungi
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702)",       # A. thaliana
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4577)",       # Zea mays (Maize)
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:39947)",      # Oryza sativa (Rice)
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83332)",      # M. tuberculosis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2157)",       # Archaea reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2)",          # Bacteria reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:5794)",       # Plasmodium
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:44689)",      # Dictyostelium
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2)",         # Bacteria TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2759)",      # Eukaryota TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2157)",      # Archaea TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:10239)",     # Viruses
]


async def fetch_5000_new_ted_domains(
    cache: ApiCache,
    existing_ted_ids: set[str],
    target_count: int = 5000,
    seed: int = 12345,
) -> list[TedDomainMember]:
    """Fetch 5,000 new random TED domains distinct from existing ones."""
    rng = random.Random(seed)
    headers = {"User-Agent": "ted-namer-benchmark-10000/0.1.0"}
    candidate_uniprots: list[str] = []

    async with httpx.AsyncClient(headers=headers, follow_redirects=True, timeout=30.0) as http:
        for q in _UNIPROT_SAMPLE_QUERIES:
            url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": q,
                "format": "json",
                "size": "500",
                "fields": "accession",
            }
            try:
                resp = await http.get(url, params=params)
                if resp.status_code == 200:
                    results = resp.json().get("results") or []
                    accs = [r["primaryAccession"] for r in results if r.get("primaryAccession")]
                    rng.shuffle(accs)
                    candidate_uniprots.extend(accs[:350])
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
        sem = asyncio.Semaphore(20)

        async def _get_ted(acc: str) -> list[TedDomainMember]:
            async with sem:
                return await ted_client.get_domains_for_uniprot(acc)

        sampled_domains: list[TedDomainMember] = []
        seen_new_ids: set[str] = set(existing_ted_ids)
        batch_size = 50

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
            if len(sampled_domains) % 500 < 60 or i % 150 == 0:
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
        "tier_superfamily": tiers.get("superfamily", 0) / n,
        "tier_family_context": tiers.get("protein_family_context", 0) / n,
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

    # Load previously sampled 5,000 domains (Batches 1-5) to avoid any overlap
    existing_ted_ids: set[str] = set()
    b5k_results_path = examples_dir / "random_5000_results.json"
    batch_1_to_5_results: list[ClusterNamingResult] = []
    if b5k_results_path.exists():
        with open(b5k_results_path, encoding="utf-8") as f:
            b5k_data = json.load(f)
            batch_1_to_5_results = [ClusterNamingResult.model_validate(x) for x in b5k_data]
            for r in batch_1_to_5_results:
                for m in r.members:
                    existing_ted_ids.add(m.domain.ted_id)
        print(f"Loaded existing 5,000 domains from Batches 1 to 5 ({len(existing_ted_ids)} unique TED IDs).")

    t0 = time.perf_counter()
    print("\nSampling 5,000 NEW random TED domains across diverse taxa (distinct from Batches 1-5)...")
    new_domains = await fetch_5000_new_ted_domains(
        cache=cache,
        existing_ted_ids=existing_ted_ids,
        target_count=5000,
        seed=424242,
    )
    print(f"Successfully collected {len(new_domains)} new random TED domains in {time.perf_counter() - t0:.1f}s.")

    # Split into 5 distinct batches of 1,000 domains (Batches 6, 7, 8, 9, 10)
    batches = [
        ("Batch 6 (5,001-6,000)", new_domains[0:1000]),
        ("Batch 7 (6,001-7,000)", new_domains[1000:2000]),
        ("Batch 8 (7,001-8,000)", new_domains[2000:3000]),
        ("Batch 9 (8,001-9,000)", new_domains[3000:4000]),
        ("Batch 10 (9,001-10,000)", new_domains[4000:5000]),
    ]

    pipeline = TedNamingPipeline(use_cache=True, concurrency=25, use_llm=False)
    all_new_results: list[ClusterNamingResult] = []
    all_batch_stats: list[tuple[str, dict[str, any]]] = []

    # Record Batches 1 to 5 stats
    if len(batch_1_to_5_results) >= 5000:
        all_batch_stats.append(("Batch 1 (1-1k)", compute_batch_stats(batch_1_to_5_results[0:1000])))
        all_batch_stats.append(("Batch 2 (1-2k)", compute_batch_stats(batch_1_to_5_results[1000:2000])))
        all_batch_stats.append(("Batch 3 (2-3k)", compute_batch_stats(batch_1_to_5_results[2000:3000])))
        all_batch_stats.append(("Batch 4 (3-4k)", compute_batch_stats(batch_1_to_5_results[3000:4000])))
        all_batch_stats.append(("Batch 5 (4-5k)", compute_batch_stats(batch_1_to_5_results[4000:5000])))

    # Process each new 1,000-domain batch sequentially and record statistics
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
        all_batch_stats.append((b_name, b_stats))

        # Checkpoint incremental saves for the new batches (5,000 new domains)
        save_output(all_new_results, examples_dir / "random_5000_additional_results.json")
        save_output(all_new_results, examples_dir / "random_5000_additional_results.tsv")

    # Combine all 10 batches: Batches 1-5 (5,000) + Batches 6-10 (5,000) = 10,000 Total Domains!
    total_10000_results = batch_1_to_5_results + all_new_results
    total_10000_stats = compute_batch_stats(total_10000_results)

    print("\nSaving complete 10,000 domain outputs...")
    save_output(total_10000_results, examples_dir / "random_10000_results.json")
    save_output(total_10000_results, examples_dir / "random_10000_results.tsv")
    save_output(total_10000_results, examples_dir / "random_10000_report.md")

    # Print comprehensive comparison table across all 10 batches and 10,000-domain total
    print("\n" + "=" * 135)
    print("STATISTICAL COMPARISON ACROSS 10 INDEPENDENT BATCHES (10,000 TED DOMAINS TOTAL)")
    print("=" * 135)
    header = (
        f"{'Metric':<28} | "
        + " | ".join(f"{b[0][:8]:^8}" for b in all_batch_stats)
        + f" | {'TOTAL 10K':^9}"
    )
    print(header)
    print("-" * len(header))

    metrics_to_show = [
        ("Domain Count", "count", "{:^8d}"),
        ("Mean Confidence", "mean_confidence", "{:^8.3f}"),
        ("Tier 1: Specific Domain", "tier_specific_domain", "{:^8.1%}"),
        ("Tier 2: Superfamily", "tier_superfamily", "{:^8.1%}"),
        ("Tier 3: Family Context", "tier_family_context", "{:^8.1%}"),
        ("Tier 4: Topology Fold", "tier_topology_fold", "{:^8.1%}"),
        ("CATH: Official (4-level)", "cath_official", "{:^8.1%}"),
        ("CATH: TED T-level (.T10)", "cath_t_level", "{:^8.1%}"),
        ("CATH: Direct 3D/SSE Rules", "cath_direct_3d_rules", "{:^8.1%}"),
        ("CATH: Gene3D Overlap", "cath_g3dsa", "{:^8.1%}"),
        ("Hard Cases (Conf < 0.50)", "conf_below_50_pct", "{:^8.1%}"),
        ("Low pLDDT (< 60)", "low_plddt_pct", "{:^8.1%}"),
    ]

    for label, key, fmt in metrics_to_show:
        row_vals = []
        for _, s in all_batch_stats:
            val = s.get(key, 0)
            row_vals.append(fmt.format(val))
        tot_val = total_10000_stats.get(key, 0)
        tot_str = fmt.format(tot_val)
        print(f"{label:<28} | " + " | ".join(row_vals) + f" | {tot_str}")

    print("=" * 135)
    print("Benchmark complete! 10,000 TED domains benchmark generated successfully.")


if __name__ == "__main__":
    asyncio.run(main())
