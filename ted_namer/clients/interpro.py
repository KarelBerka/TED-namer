"""Client for EBI InterPro API to retrieve sequence and structural domain annotations."""

from __future__ import annotations

import asyncio
from typing import Any, Optional

import httpx

from ted_namer.cache import ApiCache
from ted_namer.models import OverlapCategory, OverlappingAnnotation, SegmentRange
from ted_namer.overlap import compute_overlap_metrics, format_segments

INTERPRO_API_BASE = "https://www.ebi.ac.uk/interpro/api"


class InterProClient:
    """Fetches all InterPro and member-database matches for a UniProt protein and computes TED overlap."""

    def __init__(self, http_client: httpx.AsyncClient, cache: Optional[ApiCache] = None) -> None:
        self.http = http_client
        self.cache = cache

    async def get_protein_entries(self, uniprot_acc: str) -> list[dict[str, Any]]:
        """Retrieve all InterPro + member DB entries for a UniProt protein."""
        acc = uniprot_acc.strip().upper()
        if self.cache:
            cached = self.cache.get("interpro_protein_all", acc)
            if cached is not None:
                return cached

        url = f"{INTERPRO_API_BASE}/entry/all/protein/uniprot/{acc}?page_size=200"
        results: list[dict[str, Any]] = []
        for attempt in range(3):
            try:
                resp = await self.http.get(url, headers={"Accept": "application/json"}, timeout=20.0)
                if resp.status_code == 404:
                    return []
                if resp.status_code in {429, 500, 502, 503, 504}:
                    retry_after = float(resp.headers.get("Retry-After", 1.0 * (attempt + 1)))
                    await asyncio.sleep(min(retry_after, 5.0))
                    continue
                if resp.status_code == 200:
                    payload = resp.json()
                    if isinstance(payload, dict):
                        results = payload.get("results") or []
                    break
            except Exception:
                if attempt < 2:
                    await asyncio.sleep(0.5 * (attempt + 1))
                else:
                    results = []

        if self.cache and results:
            self.cache.set("interpro_protein_all", acc, results)
        return results

    async def get_overlapping_annotations(
        self,
        uniprot_acc: str,
        ted_chopping: str,
    ) -> list[OverlappingAnnotation]:
        """Compute residue overlap between a TED domain chopping and all InterPro/member DB hits on the protein."""
        raw_entries = await self.get_protein_entries(uniprot_acc)
        if not raw_entries:
            return []

        # Build lookup map for integrated InterPro entries present in the same payload
        ipr_name_map: dict[str, str] = {}
        for item in raw_entries:
            meta = item.get("metadata") or {}
            if meta.get("source_database") == "interpro":
                ipr_acc = meta.get("accession")
                ipr_name = meta.get("name")
                if isinstance(ipr_name, dict):
                    ipr_name = ipr_name.get("name")
                if ipr_acc and ipr_name:
                    ipr_name_map[str(ipr_acc)] = str(ipr_name)

        overlapping: list[OverlappingAnnotation] = []
        for item in raw_entries:
            meta = item.get("metadata") or {}
            source_db = str(meta.get("source_database") or "").lower()
            accession = str(meta.get("accession") or "")
            raw_name = meta.get("name")
            if isinstance(raw_name, dict):
                short_name = raw_name.get("short")
                name = raw_name.get("name") or short_name or ""
            else:
                name = str(raw_name or "")
                short_name = None

            integrated_ipr = meta.get("integrated")
            integrated_ipr_str = str(integrated_ipr) if integrated_ipr else None
            integrated_ipr_name = ipr_name_map.get(integrated_ipr_str) if integrated_ipr_str else None

            if not name and integrated_ipr_name:
                name = integrated_ipr_name
            if not name:
                continue

            annot_type = str(meta.get("type") or "domain").lower()
            proteins_list = item.get("proteins") or []

            for prot in proteins_list:
                locations = prot.get("entry_protein_locations") or []
                for loc in locations:
                    fragments = loc.get("fragments") or []
                    segs: list[SegmentRange] = []
                    for frag in fragments:
                        if isinstance(frag, dict) and frag.get("start") and frag.get("end"):
                            segs.append(SegmentRange(start=int(frag["start"]), end=int(frag["end"])))
                    if not segs:
                        continue
                    annot_chopping = format_segments(segs)
                    iou, ted_cov, annot_cov, category = compute_overlap_metrics(ted_chopping, segs)
                    if category == OverlapCategory.DISJOINT:
                        continue
                    overlapping.append(
                        OverlappingAnnotation(
                            source_db=source_db,
                            accession=accession,
                            name=name,
                            short_name=short_name,
                            integrated_ipr=integrated_ipr_str,
                            integrated_ipr_name=integrated_ipr_name,
                            annotation_type=annot_type,
                            chopping=annot_chopping,
                            iou=iou,
                            ted_coverage=ted_cov,
                            annot_coverage=annot_cov,
                            overlap_category=category,
                        )
                    )

        overlapping.sort(key=lambda a: (a.iou, a.ted_coverage), reverse=True)
        return overlapping
