"""Client for TED (The Encyclopedia of Domains) REST API (ted.cathdb.info)."""

from __future__ import annotations

import asyncio
import re
from typing import Optional

import httpx

from ted_namer.cache import ApiCache
from ted_namer.models import TedDomainMember

TED_API_BASE = "https://ted.cathdb.info/api/v1"

# Matches e.g. "AF-P04637-F1-model_v4_TED01", "P04637_TED01", or "P04637-TED01"
_TED_ID_RE = re.compile(
    r"^(?:AF-)?([A-NR-Z0-9]{6,10})(?:-F\d+-model_v\d+)?[_-](TED\d+)$",
    re.IGNORECASE,
)
_UNIPROT_RE = re.compile(r"^([A-NR-Z0-9]{6,10})$", re.IGNORECASE)


def parse_ted_identifier(identifier: str) -> tuple[str, Optional[str]]:
    """Parse a TED ID or UniProt accession into (uniprot_acc, optional_ted_suffix e.g. 'TED01')."""
    clean = identifier.strip()
    m = _TED_ID_RE.match(clean)
    if m:
        return m.group(1).upper(), m.group(2).upper()
    m_uni = _UNIPROT_RE.match(clean)
    if m_uni:
        return m_uni.group(1).upper(), None
    # Fallback if custom format containing _TEDxx
    if "_TED" in clean.upper():
        parts = clean.upper().rsplit("_TED", 1)
        ted_suffix = f"TED{parts[1]}"
        uni_part = parts[0].replace("AF-", "").split("-F")[0]
        return uni_part, ted_suffix
    return clean.upper(), None


class TedClient:
    """Async client for querying TED domain summaries by UniProt accession or TED ID."""

    def __init__(self, http_client: httpx.AsyncClient, cache: Optional[ApiCache] = None) -> None:
        self.http = http_client
        self.cache = cache

    async def get_domains_for_uniprot(self, uniprot_acc: str) -> list[TedDomainMember]:
        """Retrieve all TED domains for a given UniProt accession."""
        acc = uniprot_acc.strip().upper()
        if self.cache:
            cached = self.cache.get("ted_uniprot_summary_v2", acc)
            if cached is not None:
                return [TedDomainMember.model_validate(item) for item in cached]

        url = f"{TED_API_BASE}/uniprot/summary/{acc}"
        data: list[dict] = []
        for attempt in range(3):
            try:
                resp = await self.http.get(url, timeout=20.0)
                if resp.status_code == 404:
                    return []
                if resp.status_code in {429, 500, 502, 503, 504}:
                    await asyncio.sleep(0.5 * (attempt + 1))
                    continue
                resp.raise_for_status()
                payload = resp.json()
                data = payload.get("data", []) if isinstance(payload, dict) else []
                break
            except Exception:
                if attempt < 2:
                    await asyncio.sleep(0.5 * (attempt + 1))
                else:
                    return []

        domains: list[TedDomainMember] = []
        for item in data:
            domains.append(
                TedDomainMember(
                    ted_id=item.get("ted_id", ""),
                    uniprot_acc=item.get("uniprot_acc", acc),
                    chopping=item.get("chopping", ""),
                    nres_domain=int(item.get("nres_domain") or 0),
                    num_segments=int(item.get("num_segments") or 1),
                    plddt=float(item.get("plddt") or 0.0),
                    num_helix_strand_turn=int(item.get("num_helix_strand_turn") or 0),
                    num_helix=int(item.get("num_helix") or 0),
                    num_strand=int(item.get("num_strand") or 0),
                    num_helix_strand=int(item.get("num_helix_strand") or 0),
                    num_turn=int(item.get("num_turn") or 0),
                    cath_label=str(item.get("cath_label") or "-"),
                    cath_assignment_level=str(item.get("cath_assignment_level") or "-"),
                    cath_assignment_method=str(item.get("cath_assignment_method") or "-"),
                    packing_density=item.get("packing_density"),
                    norm_rg=item.get("norm_rg"),
                    tax_common_name=item.get("tax_common_name"),
                    tax_scientific_name=item.get("tax_scientific_name"),
                    tax_lineage=item.get("tax_lineage"),
                )
            )

        if self.cache and domains:
            self.cache.set("ted_uniprot_summary_v2", acc, [d.model_dump() for d in domains])
        return domains

    async def get_domain_pdb(self, ted_id: str) -> Optional[str]:
        """Retrieve the 3D coordinate PDB text for a single TED domain from /api/v1/files/{ted_id}.pdb."""
        clean_id = ted_id.strip()
        if not clean_id:
            return None
        if self.cache:
            cached = self.cache.get("ted_domain_pdb", clean_id)
            if cached is not None:
                return str(cached) if cached else None

        url = f"{TED_API_BASE}/files/{clean_id}.pdb"
        for attempt in range(2):
            try:
                resp = await self.http.get(url, timeout=15.0)
                if resp.status_code == 404:
                    return None
                if resp.status_code in {429, 500, 502, 503, 504}:
                    await asyncio.sleep(0.4 * (attempt + 1))
                    continue
                if resp.status_code == 200 and "ATOM" in resp.text:
                    pdb_text = resp.text
                    if self.cache:
                        # Store only ATOM C-alpha lines to keep cache compact and fast
                        ca_lines = "\n".join(
                            line
                            for line in pdb_text.splitlines()
                            if line.startswith("ATOM  ") and line[12:16].strip() == "CA"
                        )
                        self.cache.set("ted_domain_pdb", clean_id, ca_lines)
                        return ca_lines
                    return pdb_text
            except Exception:
                if attempt == 0:
                    await asyncio.sleep(0.4)
        return None

    async def resolve_identifier(self, identifier: str) -> list[TedDomainMember]:
        """Resolve either a specific TED ID (e.g. AF-P04637-F1-model_v4_TED01) or a full UniProt ID."""
        uniprot_acc, ted_suffix = parse_ted_identifier(identifier)
        all_domains = await self.get_domains_for_uniprot(uniprot_acc)
        if not ted_suffix:
            return all_domains
        matched = [
            d
            for d in all_domains
            if d.ted_id.upper().endswith(f"_{ted_suffix}") or d.short_ted_num.upper() == ted_suffix
        ]
        return matched
