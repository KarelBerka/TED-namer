"""Optional LLM consensus naming synthesizer supporting OpenAI, Gemini, and OpenAI-compatible APIs."""

from __future__ import annotations

import json
import os
import re
from typing import Any, Optional

import httpx

from ted_namer.aggregator import make_short_label
from ted_namer.models import ClusterNamingResult

_SYSTEM_PROMPT = """You are an expert structural bioinformatician curating standardized domain names for TED (The Encyclopedia of Domains, ted.cathdb.info).
Given a JSON evidence summary for a cluster of homologous TED structural domains (including overlapping CDD, Pfam, InterPro, SMART, PROSITE, CATH H/T-level classifications, and UniProt protein context), synthesize the single most accurate, concise, and standardized biological domain name.

Guidelines:
1. Prefer specific functional/structural domain names (e.g., "p53 DNA-binding domain", "P-loop NTPase domain") when high-IoU sequence/structure domain hits exist.
2. Avoid redundant or noisy suffixes (use "domain", "repeat", or "fold" appropriately; do not output ALL-CAPS database strings).
3. If the cluster has no specific domain match but has a CATH superfamily or topology, or is an unannotated region of a known protein family, combine the protein/superfamily context with the domain's relative position or fold (e.g., "BRCA1 C-terminal BRCT domain" or "Immunoglobulin-like beta-sandwich domain").
4. Output strictly valid JSON with keys:
   - "consensus_name": string (concise canonical name, 2-6 words)
   - "short_name": string (compact identifier <= 24 chars, e.g. "p53_DNA-bd")
   - "confidence_score": number between 0.0 and 1.0
   - "rationale": string (1-2 sentences citing the key database accessions, cluster frequency, and residue IoU)
"""


def build_evidence_dossier(result: ClusterNamingResult) -> dict[str, Any]:
    """Create a compact JSON dossier of cluster evidence for the LLM prompt."""
    return {
        "cluster_id": result.cluster_id,
        "member_count": result.member_count,
        "dominant_cath_label": result.dominant_cath_label,
        "cath_summary": result.cath_summary.model_dump() if result.cath_summary else None,
        "consensus_relative_position": result.consensus_relative_position,
        "taxonomic_summary": result.taxonomic_summary,
        "rule_based_baseline": {
            "consensus_name": result.consensus_name,
            "short_name": result.short_name,
            "confidence_score": result.confidence_score,
            "specificity_tier": result.specificity_tier.value,
        },
        "top_candidates": [
            {
                "canonical_name": c.canonical_name,
                "tier": c.specificity_tier.value,
                "databases": c.source_databases,
                "accessions": c.supporting_accessions[:5],
                "cluster_frequency": c.cluster_frequency,
                "mean_iou": c.mean_iou,
                "mean_ted_coverage": c.mean_ted_coverage,
                "composite_score": c.composite_score,
            }
            for c in result.top_candidates[:6]
        ],
    }


async def refine_with_llm(
    result: ClusterNamingResult,
    http_client: httpx.AsyncClient,
    provider: str = "auto",
    model: Optional[str] = None,
    base_url: Optional[str] = None,
    api_key: Optional[str] = None,
) -> ClusterNamingResult:
    """Refine a ClusterNamingResult using an LLM if credentials are available; otherwise return rule-based result."""
    resolved_provider, resolved_key, resolved_url, resolved_model = _resolve_llm_config(
        provider=provider,
        api_key=api_key,
        base_url=base_url,
        model=model,
    )
    if not resolved_provider:
        return result

    dossier = build_evidence_dossier(result)
    user_prompt = f"Cluster Evidence Dossier:\n{json.dumps(dossier, indent=2)}"

    try:
        if resolved_provider == "gemini":
            parsed = await _call_gemini(
                http_client=http_client,
                api_key=resolved_key or "",
                model=resolved_model or "gemini-2.5-flash",
                user_prompt=user_prompt,
            )
        else:
            parsed = await _call_openai_compatible(
                http_client=http_client,
                api_key=resolved_key or "",
                base_url=resolved_url or "https://api.openai.com/v1",
                model=resolved_model or "gpt-4o-mini",
                user_prompt=user_prompt,
            )

        if parsed and isinstance(parsed, dict) and parsed.get("consensus_name"):
            result.consensus_name = str(parsed["consensus_name"]).strip()
            raw_short = str(parsed.get("short_name") or "").strip()
            result.short_name = raw_short[:24] if raw_short else make_short_label(result.consensus_name)
            if parsed.get("confidence_score") is not None:
                result.confidence_score = round(min(1.0, max(0.0, float(parsed["confidence_score"]))), 3)
            if parsed.get("rationale"):
                result.rationale = str(parsed["rationale"]).strip()
            result.naming_method = f"llm_consensus ({resolved_provider})"
    except Exception:
        # Gracefully retain rule-based result on any network/API error
        pass

    return result


def _resolve_llm_config(
    provider: str,
    api_key: Optional[str],
    base_url: Optional[str],
    model: Optional[str],
) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    prov = (provider or "auto").lower().strip()
    if prov in {"none", "off", "false", "rule_based"}:
        return None, None, None, None

    gemini_key = api_key or os.environ.get("GEMINI_API_KEY") or os.environ.get("GOOGLE_API_KEY")
    openai_key = api_key or os.environ.get("OPENAI_API_KEY")

    if prov == "gemini" or (prov == "auto" and gemini_key and not openai_key):
        if not gemini_key:
            return None, None, None, None
        return "gemini", gemini_key, base_url, model or "gemini-2.5-flash"

    if prov == "openai" or (prov == "auto" and openai_key):
        if not openai_key and not base_url:
            return None, None, None, None
        return "openai", openai_key, base_url or "https://api.openai.com/v1", model or "gpt-4o-mini"

    if base_url:
        return "openai_compatible", api_key or "local", base_url, model or "llama3"

    return None, None, None, None


async def _call_openai_compatible(
    http_client: httpx.AsyncClient,
    api_key: str,
    base_url: str,
    model: str,
    user_prompt: str,
) -> Optional[dict[str, Any]]:
    url = f"{base_url.rstrip('/')}/chat/completions"
    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": _SYSTEM_PROMPT},
            {"role": "user", "content": user_prompt},
        ],
        "temperature": 0.1,
        "response_format": {"type": "json_object"},
    }
    resp = await http_client.post(url, headers=headers, json=payload, timeout=25.0)
    resp.raise_for_status()
    content = resp.json()["choices"][0]["message"]["content"]
    return _extract_json(content)


async def _call_gemini(
    http_client: httpx.AsyncClient,
    api_key: str,
    model: str,
    user_prompt: str,
) -> Optional[dict[str, Any]]:
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={api_key}"
    payload = {
        "systemInstruction": {"parts": [{"text": _SYSTEM_PROMPT}]},
        "contents": [{"parts": [{"text": user_prompt}]}],
        "generationConfig": {
            "temperature": 0.1,
            "responseMimeType": "application/json",
        },
    }
    resp = await http_client.post(url, json=payload, timeout=25.0)
    resp.raise_for_status()
    candidates = resp.json().get("candidates") or []
    if not candidates:
        return None
    text = candidates[0]["content"]["parts"][0]["text"]
    return _extract_json(text)


def _extract_json(text: str) -> Optional[dict[str, Any]]:
    clean = text.strip()
    clean = re.sub(r"^```(?:json)?\s*", "", clean)
    clean = re.sub(r"\s*```$", "", clean)
    return json.loads(clean)
