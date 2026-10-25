"""Client for resolving CATH superfamily (H-level) and topology (T-level) codes into human-readable names."""

from __future__ import annotations

from typing import Optional

import httpx

from ted_namer.cache import ApiCache
from ted_namer.cath_classifier import is_estimated_cath_code, strip_t_prefixes
from ted_namer.models import CathResolution

CATH_API_BASE = "https://www.cathdb.info/version/v4_3_0/api/rest"
INTERPRO_API_BASE = "https://www.ebi.ac.uk/interpro/api"


class CathClient:
    """Resolves CATH codes (both official '2.60.40.720' and T-estimated '3.T40.T50.T10') via InterPro and CATH APIs."""

    def __init__(self, http_client: httpx.AsyncClient, cache: Optional[ApiCache] = None) -> None:
        self.http = http_client
        self.cache = cache

    async def resolve_cath_label(
        self,
        cath_label: str,
        raw_ted_cath_label: str = "-",
        estimation_source: str = "ted_official",
    ) -> Optional[CathResolution]:
        """Resolve a CATH code (official or 'T'-prefixed estimate) to its hierarchy of names."""
        raw_code = (cath_label or "").strip()
        if not raw_code or raw_code == "-":
            return None

        display_code = raw_code.split(",")[0].strip()
        numeric_code = strip_t_prefixes(display_code)
        if not numeric_code or numeric_code == "-":
            return None

        cache_key = f"{display_code}|{estimation_source}"
        if self.cache:
            cached = self.cache.get("cath_resolution_v2", cache_key)
            if cached is not None:
                return CathResolution.model_validate(cached)

        display_parts = [p for p in display_code.split(".") if p]
        num_parts = [p for p in numeric_code.split(".") if p]
        depth = len(num_parts)
        estimated = is_estimated_cath_code(display_code) or estimation_source != "ted_official"

        assignment_level = (
            "inferred"
            if estimated
            else ("H" if depth >= 4 else ("T" if depth == 3 else "-"))
        )
        t_code_numeric = ".".join(num_parts[:3]) if depth >= 3 else None
        t_code_display = ".".join(display_parts[:3]) if len(display_parts) >= 3 else t_code_numeric
        a_code_numeric = ".".join(num_parts[:2]) if depth >= 2 else None
        c_code_numeric = num_parts[0] if depth >= 1 else None

        resolution = CathResolution(
            cath_label=display_code,
            numeric_cath_code=numeric_code,
            raw_ted_cath_label=raw_ted_cath_label,
            is_estimated=estimated,
            estimation_source=estimation_source,
            assignment_level=assignment_level,
            topology_code=t_code_display,
        )

        # 1. Populate H-level superfamily name only when H-code comes from official TED or G3DSA match
        if depth >= 4 and estimation_source in {"ted_official", "g3dsa_overlap", "ted_t_plus_g3dsa"}:
            h_code_numeric = ".".join(num_parts[:4])
            await self._populate_h_level_from_interpro(h_code_numeric, resolution)
            await self._populate_funfams(h_code_numeric, resolution)

        # 2. Populate T, A, C names from CATH cathtree using the underlying numeric codes
        if t_code_numeric:
            resolution.topology_name = await self._fetch_cathtree_name(t_code_numeric, 3)
        if a_code_numeric:
            resolution.architecture_name = await self._fetch_cathtree_name(a_code_numeric, 2)
        if c_code_numeric:
            resolution.class_name = await self._fetch_cathtree_name(c_code_numeric, 1)

        if self.cache:
            self.cache.set("cath_resolution_v2", cache_key, resolution.model_dump())

        return resolution

    async def _populate_h_level_from_interpro(self, h_code: str, resolution: CathResolution) -> None:
        url = f"{INTERPRO_API_BASE}/entry/cathgene3d/G3DSA:{h_code}"
        try:
            resp = await self.http.get(url, headers={"Accept": "application/json"}, timeout=15.0)
            if resp.status_code != 200:
                return
            metadata = resp.json().get("metadata", {})
            name_obj = metadata.get("name") or {}
            direct_name = name_obj.get("name") if isinstance(name_obj, dict) else name_obj
            integrated_ipr = metadata.get("integrated")
            rep_struct = metadata.get("representative_structure") or {}
            if isinstance(rep_struct, dict) and rep_struct.get("name"):
                resolution.representative_pdb_title = rep_struct.get("name")

            if isinstance(integrated_ipr, str) and integrated_ipr:
                resolution.integrated_ipr = integrated_ipr

            is_pdb_header_style = bool(
                direct_name and any(tok in str(direct_name) for tok in ("; Chain", "; domain", "subunit A;"))
            )
            if direct_name and not is_pdb_header_style:
                resolution.superfamily_name = str(direct_name).rstrip(".")
                return

            if resolution.integrated_ipr:
                ipr_url = f"{INTERPRO_API_BASE}/entry/interpro/{resolution.integrated_ipr}"
                ipr_resp = await self.http.get(ipr_url, headers={"Accept": "application/json"}, timeout=15.0)
                if ipr_resp.status_code == 200:
                    ipr_meta = ipr_resp.json().get("metadata", {})
                    ipr_name_obj = ipr_meta.get("name") or {}
                    ipr_name = ipr_name_obj.get("name") if isinstance(ipr_name_obj, dict) else ipr_name_obj
                    if ipr_name:
                        resolution.superfamily_name = str(ipr_name).rstrip(".")
                        return

            if direct_name:
                resolution.superfamily_name = str(direct_name).rstrip(".")
        except Exception:
            return

    async def _fetch_cathtree_name(self, cath_id: str, depth: int) -> Optional[str]:
        if self.cache:
            cached_name = self.cache.get("cath_tree_node", f"{cath_id}:{depth}")
            if cached_name is not None:
                return cached_name or None

        url = f"{CATH_API_BASE}/cathtree/from_cath_id_to_depth/{cath_id}/{depth}"
        try:
            resp = await self.http.get(url, headers={"Accept": "application/json"}, timeout=12.0)
            if resp.status_code != 200:
                return None
            data = resp.json()
            name = (data.get("name") or "").strip().rstrip(".")
            if self.cache:
                self.cache.set("cath_tree_node", f"{cath_id}:{depth}", name)
            return name or None
        except Exception:
            return None

    async def _populate_funfams(self, h_code: str, resolution: CathResolution) -> None:
        url = f"{CATH_API_BASE}/superfamily/{h_code}/funfam"
        try:
            resp = await self.http.get(url, headers={"Accept": "application/json"}, timeout=12.0)
            if resp.status_code != 200:
                return
            data = resp.json().get("data", [])
            if not isinstance(data, list):
                return
            seen: list[str] = []
            for item in data:
                if not isinstance(item, dict):
                    continue
                ff_name = (item.get("name") or "").strip()
                if ff_name and ff_name.lower() not in {"uncharacterized protein", "unassigned"}:
                    if ff_name not in seen:
                        seen.append(ff_name)
                if len(seen) >= 5:
                    break
            resolution.funfam_names = seen
        except Exception:
            return
