"""Client for UniProtKB REST API to fetch protein metadata, domain features, family comments, and GO terms."""

from __future__ import annotations

import asyncio
import re
from typing import Any, Optional

import httpx

from ted_namer.cache import ApiCache
from ted_namer.models import OverlapCategory, OverlappingAnnotation, ProteinContext, SegmentRange
from ted_namer.overlap import compute_overlap_metrics

UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"

_RELEVANT_FEATURE_TYPES = {
    "domain": "domain",
    "repeat": "repeat",
    "zinc finger": "domain",
    "dna binding": "domain",
    "region": "region",
    "motif": "motif",
    "coiled coil": "region",
}

_BELONGS_TO_RE = re.compile(r"Belongs to the\s+(.+?)(?:\.\s|$)", re.IGNORECASE)


class UniProtClient:
    """Fetches UniProtKB JSON records and extracts protein context and residue-overlapping features."""

    def __init__(self, http_client: httpx.AsyncClient, cache: Optional[ApiCache] = None) -> None:
        self.http = http_client
        self.cache = cache

    async def prefetch_uniprot_records(self, accessions: list[str], semaphore: Optional[asyncio.Semaphore] = None) -> None:
        """Batch-fetch up to 100 UniProtKB JSON records per HTTP request and store in cache."""
        if not self.cache:
            return
        missing = [
            acc.strip().upper()
            for acc in accessions
            if acc and self.cache.get("uniprot_record", acc.strip().upper()) is None
        ]
        if not missing:
            return

        # Deduplicate while preserving order
        unique_missing = list(dict.fromkeys(missing))
        chunk_size = 100

        async def _fetch_chunk(chunk: list[str]) -> None:
            url = f"{UNIPROT_API_BASE}/accessions"
            params = {
                "accessions": ",".join(chunk),
                "format": "json",
            }
            for attempt in range(3):
                try:
                    if semaphore:
                        async with semaphore:
                            resp = await self.http.get(url, params=params, headers={"Accept": "application/json"}, timeout=25.0)
                    else:
                        resp = await self.http.get(url, params=params, headers={"Accept": "application/json"}, timeout=25.0)
                    if resp.status_code == 200:
                        results = resp.json().get("results") or []
                        for rec in results:
                            if isinstance(rec, dict) and rec.get("primaryAccession"):
                                self.cache.set("uniprot_record", str(rec["primaryAccession"]).upper(), rec)
                        return
                    if resp.status_code in {429, 500, 502, 503, 504}:
                        await asyncio.sleep(1.0 * (attempt + 1))
                except Exception:
                    if attempt < 2:
                        await asyncio.sleep(0.5 * (attempt + 1))

        chunks = [unique_missing[i : i + chunk_size] for i in range(0, len(unique_missing), chunk_size)]
        await asyncio.gather(*(_fetch_chunk(c) for c in chunks))

    async def get_uniprot_record(self, uniprot_acc: str) -> Optional[dict[str, Any]]:
        acc = uniprot_acc.strip().upper()
        if self.cache:
            cached = self.cache.get("uniprot_record", acc)
            if cached is not None:
                return cached

        url = f"{UNIPROT_API_BASE}/{acc}.json"
        for attempt in range(3):
            try:
                resp = await self.http.get(url, headers={"Accept": "application/json"}, timeout=15.0)
                if resp.status_code == 404:
                    return None
                if resp.status_code in {429, 500, 502, 503, 504}:
                    retry_after = float(resp.headers.get("Retry-After", 1.0 * (attempt + 1)))
                    await asyncio.sleep(min(retry_after, 5.0))
                    continue
                if resp.status_code == 200:
                    record = resp.json()
                    if self.cache and isinstance(record, dict):
                        self.cache.set("uniprot_record", acc, record)
                    return record
            except Exception:
                if attempt < 2:
                    await asyncio.sleep(0.5 * (attempt + 1))
        return None

    async def get_context_and_features(
        self,
        uniprot_acc: str,
        ted_chopping: str,
    ) -> tuple[Optional[ProteinContext], list[OverlappingAnnotation]]:
        """Extract both ProteinContext and residue-overlapping UniProt features for a TED domain."""
        record = await self.get_uniprot_record(uniprot_acc)
        if not record:
            return None, []

        acc = uniprot_acc.strip().upper()
        entry_type = str(record.get("entryType") or "")
        is_reviewed = "reviewed" in entry_type.lower() and "unreviewed" not in entry_type.lower()
        entry_name = record.get("uniProtkbId")
        seq_len = int((record.get("sequence") or {}).get("length") or 0)
        organism = (record.get("organism") or {}).get("scientificName")

        # Extract recommended or submitted protein name
        prot_desc = record.get("proteinDescription") or {}
        rec_name = (
            ((prot_desc.get("recommendedName") or {}).get("fullName") or {}).get("value")
        )
        if not rec_name:
            sub_names = prot_desc.get("submissionNames") or []
            if sub_names and isinstance(sub_names[0], dict):
                rec_name = ((sub_names[0].get("fullName") or {}).get("value"))

        # Gene name
        genes = record.get("genes") or []
        gene_name = None
        if genes and isinstance(genes[0], dict):
            gene_name = ((genes[0].get("geneName") or {}).get("value"))

        # SIMILARITY comments -> protein families
        similarity_families: list[str] = []
        for comment in record.get("comments") or []:
            if isinstance(comment, dict) and comment.get("commentType") == "SIMILARITY":
                for txt in comment.get("texts") or []:
                    val = str(txt.get("value") or "")
                    m = _BELONGS_TO_RE.search(val)
                    if m:
                        fam = m.group(1).strip().rstrip(".")
                        similarity_families.append(fam)
                    elif val:
                        similarity_families.append(val.rstrip("."))

        # GO terms from cross-references
        go_mf: list[str] = []
        go_bp: list[str] = []
        for xref in record.get("uniProtKBCrossReferences") or []:
            if isinstance(xref, dict) and xref.get("database") == "GO":
                for prop in xref.get("properties") or []:
                    if prop.get("key") == "GoTerm":
                        val = str(prop.get("value") or "")
                        if val.startswith("F:"):
                            go_mf.append(val[2:].strip())
                        elif val.startswith("P:"):
                            go_bp.append(val[2:].strip())

        context = ProteinContext(
            uniprot_acc=acc,
            entry_name=entry_name,
            recommended_name=rec_name,
            gene_name=gene_name,
            is_reviewed=is_reviewed,
            protein_length=seq_len,
            organism=organism,
            similarity_families=similarity_families,
            go_molecular_functions=go_mf[:8],
            go_biological_processes=go_bp[:8],
        )

        # Extract residue-overlapping UniProt features
        overlapping_features: list[OverlappingAnnotation] = []
        for feat in record.get("features") or []:
            if not isinstance(feat, dict):
                continue
            raw_ftype = str(feat.get("type") or "").strip().lower()
            mapped_type = _RELEVANT_FEATURE_TYPES.get(raw_ftype)
            if not mapped_type:
                continue
            desc = str(feat.get("description") or "").strip()
            if not desc or desc.lower() in {"disordered", "interaction", " compositional bias"}:
                continue
            # Clean trailing notes like "; in form X"
            desc_clean = desc.split(";")[0].strip()
            loc = feat.get("location") or {}
            start_val = (loc.get("start") or {}).get("value")
            end_val = (loc.get("end") or {}).get("value")
            if not start_val or not end_val:
                continue

            seg = SegmentRange(start=int(start_val), end=int(end_val))
            iou, ted_cov, annot_cov, category = compute_overlap_metrics(ted_chopping, [seg])
            if category == OverlapCategory.DISJOINT:
                continue

            # Append 'domain' if feature is a Domain/Repeat and doesn't already include domain/repeat
            display_name = desc_clean
            if raw_ftype == "domain" and not any(
                w in display_name.lower() for w in ("domain", "repeat", "region", "finger", "fold", "lobe")
            ):
                display_name = f"{display_name} domain"
            elif raw_ftype == "dna binding" and "dna" not in display_name.lower():
                display_name = f"{display_name} DNA-binding domain"

            overlapping_features.append(
                OverlappingAnnotation(
                    source_db="uniprot_feature",
                    accession=f"UP_FEAT:{acc}:{start_val}-{end_val}",
                    name=display_name,
                    annotation_type=mapped_type,
                    chopping=f"{start_val}-{end_val}",
                    iou=iou,
                    ted_coverage=ted_cov,
                    annot_coverage=annot_cov,
                    overlap_category=category,
                )
            )

        return context, overlapping_features
