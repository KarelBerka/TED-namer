"""Sample 1,000 random TED domain queries across diverse taxa and test ted-namer on them."""

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

# Diverse query buckets across the tree of life (reviewed Swiss-Prot and unreviewed TrEMBL)
_UNIPROT_SAMPLE_QUERIES = [
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606)",       # Human reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090)",      # Mouse reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7955)",       # Zebrafish reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7227)",       # Drosophila reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:6239)",       # C. elegans reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83333)",      # E. coli K-12
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:224308)",     # B. subtilis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:559292)",     # S. cerevisiae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702)",       # A. thaliana
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2157)",       # Archaea reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:2)",          # Bacteria reviewed
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2)",         # Bacteria TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2759)",      # Eukaryota TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:2157)",      # Archaea TrEMBL
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4896)",       # S. pombe reviewed
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83332)",      # M. tuberculosis reviewed
]


async def fetch_random_ted_domains(
    cache: ApiCache,
    target_count: int = 1000,
    seed: int = 2026,
) -> list[TedDomainMember]:
    """Fetch random UniProt proteins with AlphaFold models and sample 1,000 TED domains."""
    rng = random.Random(seed)
    headers = {"User-Agent": "ted-namer-benchmark-1000/0.1.0"}
    candidate_uniprots: list[str] = []

    async with httpx.AsyncClient(headers=headers, follow_redirects=True, timeout=25.0) as http:
        for q in _UNIPROT_SAMPLE_QUERIES:
            url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                "query": q,
                "format": "json",
                "size": "150",
                "fields": "accession",
            }
            try:
                resp = await http.get(url, params=params)
                if resp.status_code == 200:
                    results = resp.json().get("results") or []
                    accs = [r["primaryAccession"] for r in results if r.get("primaryAccession")]
                    rng.shuffle(accs)
                    candidate_uniprots.extend(accs[:80])
            except Exception as exc:
                print(f"Warning fetching UniProt bucket '{q}': {exc}")

        # Deduplicate while preserving shuffled diversity
        seen_u: set[str] = set()
        unique_uniprots: list[str] = []
        for u in candidate_uniprots:
            if u not in seen_u:
                seen_u.add(u)
                unique_uniprots.append(u)
        rng.shuffle(unique_uniprots)
        print(f"Collected {len(unique_uniprots)} candidate UniProt accessions across {len(_UNIPROT_SAMPLE_QUERIES)} buckets.")

        ted_client = TedClient(http_client=http, cache=cache)
        sem = asyncio.Semaphore(12)

        async def _get_ted(acc: str) -> list[TedDomainMember]:
            async with sem:
                return await ted_client.get_domains_for_uniprot(acc)

        all_sampled_domains: list[TedDomainMember] = []
        seen_ted_ids: set[str] = set()
        batch_size = 30

        for i in range(0, len(unique_uniprots), batch_size):
            batch = unique_uniprots[i : i + batch_size]
            batch_res = await asyncio.gather(*(_get_ted(a) for a in batch))
            for dom_list in batch_res:
                if not dom_list:
                    continue
                # Pick up to 2 random domains per protein for high diversity
                k = min(len(dom_list), rng.choice([1, 1, 2, 2]))
                chosen = rng.sample(dom_list, k=k)
                for d in chosen:
                    if d.ted_id not in seen_ted_ids:
                        seen_ted_ids.add(d.ted_id)
                        all_sampled_domains.append(d)
            if len(all_sampled_domains) >= target_count:
                break
            if (len(all_sampled_domains) % 150 < 30) or i % 120 == 0:
                print(f"  Sampled {len(all_sampled_domains)}/{target_count} unique TED domains so far...")

    return all_sampled_domains[:target_count]


async def main() -> None:
    t0 = time.perf_counter()
    cache = ApiCache()
    print("Sampling 1,000 random TED domains across diverse taxa & review statuses...")
    sampled_domains = await fetch_random_ted_domains(cache=cache, target_count=1000, seed=42)
    print(f"Collected {len(sampled_domains)} random TED domains in {time.perf_counter() - t0:.1f}s.")

    examples_dir = Path("examples")
    examples_dir.mkdir(exist_ok=True)
    queries_tsv = examples_dir / "random_1000_queries.tsv"

    lines = ["cluster_id\tted_id"]
    cluster_map: dict[str, list[str]] = {}
    for idx, dom in enumerate(sampled_domains, start=1):
        qid = f"Q{idx:04d}_{dom.uniprot_acc}_{dom.short_ted_num}"
        lines.append(f"{qid}\t{dom.ted_id}")
        cluster_map[qid] = [dom.ted_id]
    queries_tsv.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Saved {len(sampled_domains)} queries to {queries_tsv}")

    t1 = time.perf_counter()
    pipeline = TedNamingPipeline(
        use_cache=True,
        concurrency=15,
        use_llm=False,
    )

    last_print = 0
    t_start = time.perf_counter()

    def _progress(done: int, total: int) -> None:
        nonlocal last_print
        if done - last_print >= 50 or done == total:
            last_print = done
            elapsed_now = time.perf_counter() - t_start
            pct = (done / max(1, total)) * 100
            rate = done / max(0.1, elapsed_now)
            print(f"[{done:4d}/{total:4d}] ({pct:5.1f}%) Enriched TED domains... (elapsed {elapsed_now:5.1f}s, {rate:4.1f} dom/s)")

    print(f"\nStarting pipeline run on {len(cluster_map)} TED domain queries...")
    results = await pipeline.name_clusters(
        cluster_map=cluster_map,
        group_by="cluster_id",
        progress_callback=_progress,
    )
    elapsed = time.perf_counter() - t1
    print(f"\nNamed {len(results)} TED domain queries in {elapsed:.1f}s ({elapsed / max(1, len(results)):.3f}s/query).")

    print("Saving output files...")
    save_output(results, examples_dir / "random_1000_results.json")
    save_output(results, examples_dir / "random_1000_results.tsv")
    save_output(results, examples_dir / "random_1000_report.md")
    print("Output files saved successfully.")

    tier_counts = Counter(r.specificity_tier.value for r in results)
    cath_est_counts = Counter(
        r.members[0].domain.cath_estimation_source if r.members else "-" for r in results
    )
    avg_conf = sum(r.confidence_score for r in results) / len(results) if results else 0.0

    print("\n=== SUMMARY ACROSS 1,000 RANDOM TED QUERIES ===")
    print(f"Total Queries Named: {len(results)}")
    print(f"Mean Confidence Score: {avg_conf:.3f}")
    print("Specificity Tier Breakdown:")
    for tier, cnt in tier_counts.most_common():
        print(f"  - {tier}: {cnt} ({cnt / len(results):.1%})")
    print("CATH Classification Source Breakdown:")
    for src, cnt in cath_est_counts.most_common():
        print(f"  - {src}: {cnt} ({cnt / len(results):.1%})")


if __name__ == "__main__":
    asyncio.run(main())
