"""Sample and benchmark an additional 10x 1,000 random TED domains (Batches 11-20, 10,000 new domains)
and evaluate accuracy across Batches 1-10 (10k), Batches 11-20 (10k out-of-sample), and all 20,000 TED domains."""

from __future__ import annotations

import asyncio
import json
import random
import re
import time
from collections import Counter
from pathlib import Path

import httpx

from ted_namer.cache import ApiCache
from ted_namer.clients.ted import TedClient
from ted_namer.export import save_output
from ted_namer.models import ClusterNamingResult, TedDomainMember
from ted_namer.pipeline import TedNamingPipeline

# Expanded diverse taxonomic & length buckets across all kingdoms of life (Swiss-Prot + TrEMBL)
_UNIPROT_SAMPLE_QUERIES = [
    # Vertebrates & Model Metazoans (Reviewed & Unreviewed)
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606) AND (length:[100 TO 400])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606) AND (length:[401 TO 800])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9606) AND (length:[801 TO 2500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090) AND (length:[100 TO 500])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10090) AND (length:[501 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:10116)",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9913)",       # Bos taurus
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9823)",       # Sus scrofa
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:9031)",       # Gallus gallus
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:8364)",       # Xenopus tropicalis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7955)",       # Danio rerio
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:7227)",       # Drosophila
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:6239)",       # C. elegans
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:6656)",       # Arthropoda
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:6447)",       # Mollusca
    # Plants & Algae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702) AND (length:[100 TO 450])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3702) AND (length:[451 TO 2000])",
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:39947)",      # Oryza sativa
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4577)",       # Zea mays
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:4081)",       # Solanum lycopersicum
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3847)",       # Glycine max
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:3055)",       # Chlamydomonas
    # Fungi & Protists
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:559292)",     # S. cerevisiae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:284812)",     # S. pombe
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:237561)",     # Candida albicans
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:162425)",     # Aspergillus nidulans
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:367110)",     # Neurospora crassa
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:44689)",      # Dictyostelium
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:5833)",       # Plasmodium falciparum
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:5691)",       # Trypanosoma brucei
    # Bacteria (Diverse Phyla)
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83333)",      # E. coli K-12
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:99287)",      # Salmonella typhimurium
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:208964)",     # Pseudomonas aeruginosa
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:243277)",     # Vibrio cholerae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:224308)",     # B. subtilis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:93061)",      # S. aureus
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:171101)",     # S. pneumoniae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:83332)",      # M. tuberculosis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:246196)",     # M. smegmatis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:100226)",     # Streptomyces coelicolor
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:196627)",     # Corynebacterium glutamicum
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:1148)",       # Synechocystis
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:226186)",     # Bacteroides thetaiotaomicron
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:85962)",      # Helicobacter pylori
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:192222)",     # Campylobacter jejuni
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:243230)",     # Deinococcus radiodurans
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:243274)",     # Thermotoga maritima
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:224324)",     # Aquifex aeolicus
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:272634)",     # Mycoplasma pneumoniae
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:272561)",     # Chlamydia trachomatis
    # Archaea (Diverse Phyla)
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:243232)",     # Methanocaldococcus jannaschii
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:186497)",     # Pyrococcus furiosus
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:273057)",     # Saccharolobus solfataricus
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:64091)",      # Halobacterium salinarum
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:224325)",     # Archaeoglobus fulgidus
    "(database:alphafolddb) AND (reviewed:true) AND (taxonomy_id:188937)",     # Methanosarcina acetivorans
    # Unreviewed TrEMBL across Kingdoms & Length Bins
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:1224) AND (length:[100 TO 400])",   # Proteobacteria TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:1239) AND (length:[100 TO 400])",   # Firmicutes TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:201174) AND (length:[100 TO 400])", # Actinobacteria TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:976) AND (length:[100 TO 400])",    # Bacteroidetes TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:1117) AND (length:[100 TO 400])",   # Cyanobacteria TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:28890) AND (length:[100 TO 400])",  # Euryarchaeota TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:28889) AND (length:[100 TO 400])",  # Crenarchaeota TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:4751) AND (length:[100 TO 500])",   # Fungi TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:33090) AND (length:[100 TO 500])",  # Viridiplantae TrEMBL
    "(database:alphafolddb) AND (reviewed:false) AND (taxonomy_id:33208) AND (length:[100 TO 500])",  # Metazoa TrEMBL
]

_NEXT_LINK_RE = re.compile(r'<([^>]+)>;\s*rel="next"')


async def fetch_10000_new_ted_domains(
    cache: ApiCache,
    existing_ted_ids: set[str],
    target_count: int = 10000,
    seed: int = 20261027,
) -> list[TedDomainMember]:
    """Fetch 10,000 new random TED domains strictly disjoint from existing_ted_ids."""
    rng = random.Random(seed)
    headers = {"User-Agent": "ted-namer-benchmark-20000/0.1.0"}
    candidate_uniprots: list[str] = []

    async with httpx.AsyncClient(headers=headers, follow_redirects=True, timeout=30.0) as http:
        q_sem = asyncio.Semaphore(12)

        async def _run_uniprot_query(q: str) -> list[str]:
            accs_local: list[str] = []
            url = "https://rest.uniprot.org/uniprotkb/search"
            params: dict[str, str] | None = {
                "query": q,
                "format": "json",
                "size": "500",
                "fields": "accession",
            }
            async with q_sem:
                for _page in range(2):  # Up to 2 pages (1,000 accessions) per query
                    if not url:
                        break
                    try:
                        resp = await http.get(url, params=params)
                        params = None  # Subsequent next-link URLs already include query params
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

        query_lists = await asyncio.gather(*(_run_uniprot_query(q) for q in _UNIPROT_SAMPLE_QUERIES))
        for accs in query_lists:
            rng.shuffle(accs)
            candidate_uniprots.extend(accs)

        # Deduplicate and shuffle candidate accessions
        seen_u: set[str] = set()
        unique_uniprots: list[str] = []
        for u in candidate_uniprots:
            if u not in seen_u:
                seen_u.add(u)
                unique_uniprots.append(u)
        rng.shuffle(unique_uniprots)
        print(
            f"Collected {len(unique_uniprots)} unique candidate UniProt accessions "
            f"across {len(_UNIPROT_SAMPLE_QUERIES)} queries."
        )

        ted_client = TedClient(http_client=http, cache=cache)
        sem = asyncio.Semaphore(35)

        async def _get_ted(acc: str) -> list[TedDomainMember]:
            async with sem:
                return await ted_client.get_domains_for_uniprot(acc)

        sampled_domains: list[TedDomainMember] = []
        seen_ids: set[str] = set(existing_ted_ids)
        batch_size = 100

        for i in range(0, len(unique_uniprots), batch_size):
            batch = unique_uniprots[i : i + batch_size]
            batch_res = await asyncio.gather(*(_get_ted(a) for a in batch))
            for dom_list in batch_res:
                if not dom_list:
                    continue
                k = min(len(dom_list), rng.choice([1, 2, 2, 3]))
                chosen = rng.sample(dom_list, k=k)
                for d in chosen:
                    if d.ted_id not in seen_ids:
                        seen_ids.add(d.ted_id)
                        sampled_domains.append(d)
            if len(sampled_domains) >= target_count:
                break
            if (i // batch_size) % 5 == 0:
                print(f"  Sampled {len(sampled_domains)}/{target_count} new unique TED domains...")

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

    official_cath = sources.get("ted_official", 0)
    t_level_cath = sum(v for k, v in sources.items() if k.startswith("ted_t_"))
    g3dsa_cath = sources.get("g3dsa_overlap", 0)
    profile_cath = sources.get("profile_cath_map", 0)
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
        "cath_profile_map": profile_cath / n,
        "cath_g3dsa": g3dsa_cath / n,
        "cath_direct_3d_rules": direct_3d_sse_rules / n,
        "conf_below_50_pct": conf_below_50 / n,
        "low_plddt_pct": low_plddt / n,
    }


async def main() -> None:
    examples_dir = Path("examples")
    examples_dir.mkdir(exist_ok=True)
    cache = ApiCache()

    # Load existing 10,000 domains (Batches 1-10)
    existing_ted_ids: set[str] = set()
    b10k_path = examples_dir / "random_10000_results.json"
    batch_1_to_10_results: list[ClusterNamingResult] = []
    if b10k_path.exists():
        raw_10k = json.loads(b10k_path.read_text(encoding="utf-8"))
        batch_1_to_10_results = [ClusterNamingResult.model_validate(x) for x in raw_10k]
        for r in batch_1_to_10_results:
            for m in r.members:
                existing_ted_ids.add(m.domain.ted_id)
        print(f"Loaded existing 10,000 domains from Batches 1-10 ({len(existing_ted_ids)} unique TED IDs).")

    # Check if we already have partial checkpoints for Batches 11-20
    ckpt_path = examples_dir / "random_10000_batches_11_to_20.json"
    all_new_results: list[ClusterNamingResult] = []
    if ckpt_path.exists():
        raw_ckpt = json.loads(ckpt_path.read_text(encoding="utf-8"))
        all_new_results = [ClusterNamingResult.model_validate(x) for x in raw_ckpt]
        for r in all_new_results:
            for m in r.members:
                existing_ted_ids.add(m.domain.ted_id)
        print(f"Resuming from checkpoint with {len(all_new_results)} already processed domains in Batches 11-20.")

    remaining_needed = 10000 - len(all_new_results)
    if remaining_needed > 0:
        sampled_list_path = examples_dir / "sampled_domains_batches_11_to_20.json"
        if sampled_list_path.exists():
            raw_sampled = json.loads(sampled_list_path.read_text(encoding="utf-8"))
            all_sampled_10k = [TedDomainMember.model_validate(x) for x in raw_sampled]
            new_domains = [d for d in all_sampled_10k if d.ted_id not in existing_ted_ids][:remaining_needed]
            print(f"Loaded {len(new_domains)} remaining pre-sampled TED domains from {sampled_list_path}.", flush=True)
        else:
            t0 = time.perf_counter()
            print(f"\nSampling {remaining_needed} NEW random TED domains (strictly disjoint from existing)...", flush=True)
            new_domains = await fetch_10000_new_ted_domains(
                cache=cache,
                existing_ted_ids=existing_ted_ids,
                target_count=remaining_needed,
                seed=20261027,
            )
            sampled_list_path.write_text(
                json.dumps([d.model_dump() for d in new_domains], indent=2), encoding="utf-8"
            )
            print(f"Collected {len(new_domains)} new unique TED domains in {time.perf_counter() - t0:.1f}s.", flush=True)

        pipeline = TedNamingPipeline(use_cache=True, concurrency=60, use_llm=False)
        start_batch_idx = 11 + (len(all_new_results) // 1000)
        num_new_batches = len(new_domains) // 1000

        for b_offset in range(num_new_batches):
            b_num = start_batch_idx + b_offset
            b_domains = new_domains[b_offset * 1000 : (b_offset + 1) * 1000]
            b_name = f"Batch {b_num} ({(b_num - 1) * 1000 + 1}-{b_num * 1000})"
            print(f"\n=======================================================", flush=True)
            print(f"Running pipeline on {b_name} ({len(b_domains)} domains)...", flush=True)
            print(f"=======================================================", flush=True)
            cmap = {
                f"Q_{dom.uniprot_acc}_{dom.short_ted_num}": [dom.ted_id]
                for dom in b_domains
            }

            b_start = time.perf_counter()
            last_log = 0

            def _prog(done: int, total: int) -> None:
                nonlocal last_log
                if done - last_log >= 200 or done == total:
                    last_log = done
                    el = time.perf_counter() - b_start
                    rate = done / max(0.1, el)
                    print(f"  [{b_name}] {done:4d}/{total:4d} ({done/total*100:5.1f}%) | {el:5.1f}s | {rate:4.1f} dom/s", flush=True)

            b_res = await pipeline.name_clusters(cluster_map=cmap, group_by="cluster_id", progress_callback=_prog)
            b_elapsed = time.perf_counter() - b_start
            print(f"  Finished {b_name} in {b_elapsed:.1f}s ({len(b_res)/max(0.1, b_elapsed):.1f} dom/s).", flush=True)

            all_new_results.extend(b_res)
            save_output(all_new_results, ckpt_path)
            save_output(all_new_results, examples_dir / "random_10000_batches_11_to_20.tsv")

    # Combine all 20 batches (20,000 domains)
    total_20000_results = batch_1_to_10_results + all_new_results
    save_output(total_20000_results, examples_dir / "random_20000_results.json")
    save_output(total_20000_results, examples_dir / "random_20000_results.tsv")

    stats_1_10 = compute_batch_stats(batch_1_to_10_results)
    stats_11_20 = compute_batch_stats(all_new_results)
    stats_total_20k = compute_batch_stats(total_20000_results)

    summary = {
        "batches_1_to_10_10k": stats_1_10,
        "batches_11_to_20_10k_new": stats_11_20,
        "total_20k": stats_total_20k,
        "per_batch_11_to_20": {
            f"Batch_{11 + i}": compute_batch_stats(all_new_results[i * 1000 : (i + 1) * 1000])
            for i in range(min(10, len(all_new_results) // 1000))
        },
    }
    (examples_dir / "benchmark_20000_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("\n" + "=" * 100)
    print("20,000 TED DOMAINS BENCHMARK SUMMARY (Batches 1-10 vs Batches 11-20 vs Total 20,000)")
    print("=" * 100)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    asyncio.run(main())
