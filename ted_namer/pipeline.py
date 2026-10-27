"""End-to-end async pipeline for enriching TED domains and naming domain clusters."""

from __future__ import annotations

import asyncio
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Callable, Optional

import httpx

from ted_namer.aggregator import aggregate_cluster_evidence
from ted_namer.cache import ApiCache
from ted_namer.cath_classifier import infer_cath_classification
from ted_namer.clients.cath import CathClient
from ted_namer.clients.interpro import InterProClient
from ted_namer.clients.ted import TedClient
from ted_namer.clients.uniprot import UniProtClient
from ted_namer.llm import refine_with_llm
from ted_namer.models import ClusterNamingResult, EnrichedTedDomain, TedDomainMember
from ted_namer.overlap import compute_relative_position


class TedNamingPipeline:
    """Async pipeline coordinating TED, InterPro, CATH, UniProt, cluster aggregation, and LLM synthesis."""

    def __init__(
        self,
        cache_path: str | Path = ".ted_namer_cache.sqlite",
        use_cache: bool = True,
        concurrency: int = 6,
        use_llm: bool = True,
        llm_provider: str = "auto",
        llm_model: Optional[str] = None,
        llm_base_url: Optional[str] = None,
        llm_api_key: Optional[str] = None,
    ) -> None:
        self.cache = ApiCache(db_path=cache_path, enabled=use_cache)
        self.concurrency = max(1, concurrency)
        self.use_llm = use_llm
        self.llm_provider = llm_provider
        self.llm_model = llm_model
        self.llm_base_url = llm_base_url
        self.llm_api_key = llm_api_key

    async def enrich_ted_domain(
        self,
        member: TedDomainMember,
        interpro_client: InterProClient,
        cath_client: CathClient,
        uniprot_client: UniProtClient,
        semaphore: asyncio.Semaphore,
        ted_client: Optional[TedClient] = None,
    ) -> EnrichedTedDomain:
        """Fetch and compute all overlapping annotations, CATH resolution (including 3D & rule-based), and UniProt context."""
        async with semaphore:
            ipr_task = interpro_client.get_overlapping_annotations(
                uniprot_acc=member.uniprot_acc,
                ted_chopping=member.chopping,
            )
            uni_task = uniprot_client.get_context_and_features(
                uniprot_acc=member.uniprot_acc,
                ted_chopping=member.chopping,
            )
            ipr_annots, (prot_ctx, uni_feats) = await asyncio.gather(ipr_task, uni_task)

            all_annots = list(ipr_annots) + list(uni_feats)
            all_annots.sort(key=lambda a: (a.iou, a.ted_coverage), reverse=True)

            # First attempt CATH classification using TED official, G3DSA overlap, and Profile-to-CATH map
            inferred_cath, _numeric_cath, est_source = infer_cath_classification(
                member=member,
                overlapping_annotations=all_annots,
                pdb_text=None,
            )

            # Only fetch 3D PDB C-alpha trace when sequence/structural profiles could not resolve CATH
            if ted_client is not None and est_source.startswith("rule_"):
                pdb_text = await ted_client.get_domain_pdb(member.ted_id)
                if pdb_text:
                    inferred_cath, _numeric_cath, est_source = infer_cath_classification(
                        member=member,
                        overlapping_annotations=all_annots,
                        pdb_text=pdb_text,
                    )

            member.inferred_cath_label = inferred_cath
            member.cath_estimation_source = est_source

            cath_info = await cath_client.resolve_cath_label(
                cath_label=inferred_cath,
                raw_ted_cath_label=member.cath_label,
                estimation_source=est_source,
                overlapping_annotations=all_annots,
            )

        prot_len = prot_ctx.protein_length if prot_ctx else 0
        rel_pos = compute_relative_position(member.chopping, prot_len)

        return EnrichedTedDomain(
            domain=member,
            cath_info=cath_info,
            protein_context=prot_ctx,
            overlapping_annotations=all_annots,
            relative_position=rel_pos,
        )

    async def name_clusters(
        self,
        cluster_map: dict[str, list[str]],
        group_by: str = "cluster_id",
        max_members_per_cluster: int = 50,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> list[ClusterNamingResult]:
        """Resolve TED identifiers in cluster_map, enrich members, and compute cluster names."""
        semaphore = asyncio.Semaphore(self.concurrency)
        headers = {"User-Agent": "ted-namer/0.1.0 (https://ted.cathdb.info)"}

        limits = httpx.Limits(max_connections=max(100, self.concurrency * 3), max_keepalive_connections=50)
        async with httpx.AsyncClient(headers=headers, follow_redirects=True, limits=limits) as http:
            ted_client = TedClient(http_client=http, cache=self.cache)
            interpro_client = InterProClient(http_client=http, cache=self.cache)
            cath_client = CathClient(http_client=http, cache=self.cache)
            uniprot_client = UniProtClient(http_client=http, cache=self.cache)

            # 1. Resolve all identifiers to TedDomainMember objects concurrently
            resolved_clusters: dict[str, list[TedDomainMember]] = defaultdict(list)

            async def _resolve_one(ident: str) -> list[TedDomainMember]:
                async with semaphore:
                    return await ted_client.resolve_identifier(ident)

            async def _resolve_cluster_item(raw_cid: str, id_list: list[str]) -> tuple[str, list[list[TedDomainMember]]]:
                sampled_ids = id_list[:max_members_per_cluster]
                nested = await asyncio.gather(*(_resolve_one(i) for i in sampled_ids))
                return raw_cid, list(nested)

            cluster_resolved_items = await asyncio.gather(
                *(_resolve_cluster_item(raw_cid, id_list) for raw_cid, id_list in cluster_map.items())
            )
            for raw_cid, nested_members in cluster_resolved_items:
                for dom_list in nested_members:
                    for dom in dom_list:
                        if group_by == "cath_label":
                            target_cid = (
                                f"CATH_{dom.cath_label}"
                                if dom.cath_label and dom.cath_label != "-"
                                else f"UNASSIGNED_{dom.uniprot_acc}_{dom.short_ted_num}"
                            )
                        elif group_by == "single":
                            target_cid = "cluster_1"
                        else:
                            target_cid = raw_cid
                        # Deduplicate by ted_id within cluster
                        if all(existing.ted_id != dom.ted_id for existing in resolved_clusters[target_cid]):
                            resolved_clusters[target_cid].append(dom)

            # 2. Enrich all unique TED domains across clusters
            unique_domains: dict[str, TedDomainMember] = {}
            for dom_list in resolved_clusters.values():
                for d in dom_list:
                    unique_domains[d.ted_id] = d

            completed_count = 0
            total_enrich = len(unique_domains)

            async def _enrich_tracked(dom: TedDomainMember) -> tuple[str, EnrichedTedDomain]:
                nonlocal completed_count
                enriched = await self.enrich_ted_domain(
                    member=dom,
                    interpro_client=interpro_client,
                    cath_client=cath_client,
                    uniprot_client=uniprot_client,
                    semaphore=semaphore,
                    ted_client=ted_client,
                )
                completed_count += 1
                if progress_callback:
                    progress_callback(completed_count, total_enrich)
                return dom.ted_id, enriched

            enriched_items = await asyncio.gather(*(_enrich_tracked(d) for d in unique_domains.values()))
            enriched_lookup: dict[str, EnrichedTedDomain] = dict(enriched_items)

            # 3. Aggregate per cluster and optionally refine with LLM
            results: list[ClusterNamingResult] = []
            for cid, dom_list in resolved_clusters.items():
                enriched_members = [enriched_lookup[d.ted_id] for d in dom_list if d.ted_id in enriched_lookup]
                res = aggregate_cluster_evidence(cluster_id=cid, members=enriched_members)
                if self.use_llm:
                    res = await refine_with_llm(
                        result=res,
                        http_client=http,
                        provider=self.llm_provider,
                        model=self.llm_model,
                        base_url=self.llm_base_url,
                        api_key=self.llm_api_key,
                    )
                results.append(res)

        return results


def load_clusters_from_file(path: str | Path) -> dict[str, list[str]]:
    """Load cluster definitions from a JSON, TSV, CSV, or plain text file of TED/UniProt IDs."""
    file_path = Path(path)
    text = file_path.read_text(encoding="utf-8").strip()
    if not text:
        return {}

    suffix = file_path.suffix.lower()
    if suffix == ".json":
        data = json.loads(text)
        if isinstance(data, dict):
            return {str(k): [str(x) for x in v] for k, v in data.items() if isinstance(v, list)}
        if isinstance(data, list):
            clusters: dict[str, list[str]] = defaultdict(list)
            for row in data:
                if isinstance(row, dict):
                    cid = str(row.get("cluster_id") or row.get("cluster") or "cluster_1")
                    tid = str(row.get("ted_id") or row.get("uniprot_acc") or row.get("id") or "")
                    if tid:
                        clusters[cid].append(tid)
                elif isinstance(row, str):
                    clusters["cluster_1"].append(row)
            return dict(clusters)

    # Delimited TSV/CSV or line-by-line list
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if not lines:
        return {}

    delimiter = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    if delimiter is None:
        return {"cluster_1": lines}

    reader = csv.reader(lines, delimiter=delimiter)
    rows = list(reader)
    if not rows:
        return {}

    header = [h.strip().lower() for h in rows[0]]
    has_header = any(h in {"cluster_id", "cluster", "ted_id", "uniprot_acc", "uniprot", "domain_id"} for h in header)
    start_idx = 1 if has_header else 0

    cid_col = 0
    id_col = 1 if len(rows[0]) > 1 else 0
    if has_header:
        for idx, col_name in enumerate(header):
            if col_name in {"cluster_id", "cluster", "group", "family"}:
                cid_col = idx
            elif col_name in {"ted_id", "domain_id", "uniprot_acc", "uniprot", "id"}:
                id_col = idx

    clusters = defaultdict(list)
    for row in rows[start_idx:]:
        if not row:
            continue
        if len(row) == 1:
            clusters["cluster_1"].append(row[0].strip())
        else:
            cid = row[cid_col].strip() or "cluster_1"
            ident = row[id_col].strip()
            if ident:
                clusters[cid].append(ident)
    return dict(clusters)
