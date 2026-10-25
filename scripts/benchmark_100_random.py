"""Sample 100 random TED domain queries across diverse taxa and test ted-namer on them."""

from __future__ import annotations

import asyncio
import random
import time
from collections import Counter
from pathlib import Path

import httpx

from ted_namer.cache import ApiCache
from ted_namer.clients.ted import TedClient
from ted_namer.export import save_output
from ted_namer.models import TedDomainMember
from ted_namer.pipeline import TedNamingPipeline

# Diverse query buckets across tree of life and annotation status to get a representative random sample
_UNIPROT_SAMPLE_QUERIES = [
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606)",       # Human reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090)",      # Mouse reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83333)",      # E. coli K-12
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:224308)",     # B. subtilis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:559292)",     # S. cerevisiae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702)",       # A. thaliana
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2157)",       # Archaea reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2)",          # Bacteria reviewed
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2)",         # Bacteria TrEMBL (unreviewed)
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2759)",      # Eukaryota TrEMBL (unreviewed)
]


async def fetch_random_ted_domains(
    cache: ApiCache,
    target_count: int = 100,
    seed: int = 2026,
) -> list[TedDomainMember]:
    """Fetch random UniProt proteins with AlphaFold models and sample 100 TED domains."""
    rng = random.Random(seed)
    headers = {"User-Agent": "ted-namer-benchmark/0.1.0"}
    candidate_uniprots: list[str] = []

    async with httpx.AsyncClient(headers=headers, follow_redirects=True, timeout=20.0) as http:
        for q in _UNIPROT_SAMPLE_QUERIES:
            url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": q,
                "format": "json",
                "size": "100",
                "fields": "accession",
            }
            try:
                resp = await http.get(url, params=params)
                if resp.status_code == 200:
                    results = resp.json().get("results") or []
                    accs = [r["primaryAccession"] for r in results if r.get("primaryAccession")]
                    rng.shuffle(accs)
                    candidate_uniprots.extend(accs[:45])
            except Exception as exc:
                print(f"Warning fetching UniProt bucket '{q}': {exc}")

        # Deduplicate while preserving order, then shuffle across buckets
        seen_u: set[str] = set()
        unique_uniprots: list[str] = []
        for u in candidate_uniprots:
            if u not in seen_u:
                seen_u.add(u)
                unique_uniprots.append(u)
        rng.shuffle(unique_uniprots)

        ted_client = TedClient(http_client=http, cache=cache)
        sem = asyncio.Semaphore(8)

        async def _get_ted(acc: str) -> list[TedDomainMember]:
            async with sem:
                return await ted_client.get_domains_for_uniprot(acc)

        all_sampled_domains: list[TedDomainMember] = []
        batch_size = 24
        for i in range(0, len(unique_uniprots), batch_size):
            batch = unique_uniprots[i : i + batch_size]
            batch_res = await asyncio.gather(*(_get_ted(a) for a in batch))
            for dom_list in batch_res:
                if not dom_list:
                    continue
                # Pick 1 or 2 random TED domains per protein so 100 queries span ~75-90 distinct proteins
                k = min(len(dom_list), rng.choice([1, 1, 2]))
                chosen = rng.sample(dom_list, k=k)
                all_sampled_domains.extend(chosen)
            if len(all_sampled_domains) >= target_count:
                break

    return all_sampled_domains[:target_count]


async def main() -> None:
    t0 = time.perf_counter()
    cache = ApiCache()
    print("Sampling 100 random TED domains across diverse taxa & review statuses...")
    sampled_domains = await fetch_random_ted_domains(cache=cache, target_count=100, seed=2026)
    print(f"Collected {len(sampled_domains)} random TED domains in {time.perf_counter() - t0:.1f}s.")

    examples_dir = Path("examples")
    examples_dir.mkdir(exist_ok=True)
    queries_tsv = examples_dir / "random_100_queries.tsv"

    lines = ["cluster_id\tted_id"]
    cluster_map: dict[str, list[str]] = {}
    for idx, dom in enumerate(sampled_domains, start=1):
        qid = f"Q{idx:03d}_{dom.uniprot_acc}_{dom.short_ted_num}"
        lines.append(f"{qid}\t{dom.ted_id}")
        cluster_map[qid] = [dom.ted_id]
    queries_tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Saved {len(sampled_domains)} queries to {queries_tsv}")

    t1 = time.perf_counter()
    pipeline = TedNamingPipeline(
        use_cache=True,
        concurrency=10,
        use_llm=False,
    )
    results = await pipeline.name_clusters(cluster_map=cluster_map, group_by="cluster_id")
    elapsed = time.perf_counter() - t1
    print(f"Named {len(results)} TED domain queries in {elapsed:.1f}s ({elapsed / max(1, len(results)):.2f}s/query).")

    save_output(results, examples_dir / "random_100_results.json")
    save_output(results, examples_dir / "random_100_results.tsv")
    save_output(results, examples_dir / "random_100_report.md")

    tier_counts = Counter(r.specificity_tier.value for r in results)
    cath_est_counts = Counter(
        r.members[0].domain.cath_estimation_source if r.members else "-" for r in results
    )
    avg_conf = sum(r.confidence_score for r in results) / len(results) if results else 0.0

    print("\n=== SUMMARY ACROSS 100 RANDOM TED QUERIES ===")
    print(f"Total Queries Named: {len(results)}")
    print(f"Mean Confidence Score: {avg_conf:.3f}")
    print("Specificity Tier Breakdown:")
    for tier, cnt in tier_counts.most_common():
        print(f"  - {tier}: {cnt} ({cnt / len(results):.0%})")
    print("CATH Classification Source Breakdown:")
    for src, cnt in cath_est_counts.most_common():
        print(f"  - {src}: {cnt} ({cnt / len(results):.0%})")


if __name__ == "__main__":
    asyncio.run(main())
