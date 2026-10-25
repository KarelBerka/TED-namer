"""Live end-to-end integration tests against TED, InterPro, CATH, and UniProt REST APIs."""

import pytest

from ted_namer.models import SpecificityTier
from ted_namer.pipeline import TedNamingPipeline


@pytest.mark.asyncio
async def test_live_p53_homologous_cluster_naming(tmp_path) -> None:
    """Verify naming of the homologous p53/p63/p73 DNA-binding TED domain cluster across 3 human paralogs."""
    cache_file = tmp_path / "test_cache.sqlite"
    pipeline = TedNamingPipeline(
        cache_path=cache_file,
        use_cache=True,
        use_llm=False,
    )

    # P04637 (p53 TED01), Q9H3D4 (p63 TED01), O15350 (p73 TED01) are all DNA-binding domains (CATH 2.60.40.720)
    clusters = {
        "p53_family_DBD": [
            "AF-P04637-F1-model_v4_TED01",
            "AF-Q9H3D4-F1-model_v4_TED01",
            "AF-O15350-F1-model_v4_TED01",
        ]
    }

    results = await pipeline.name_clusters(clusters, group_by="cluster_id")
    assert len(results) == 1
    res = results[0]
    assert res.cluster_id == "p53_family_DBD"
    assert res.member_count == 3
    assert res.dominant_cath_label == "2.60.40.720"
    assert res.specificity_tier == SpecificityTier.SPECIFIC_DOMAIN
    assert "dna-binding" in res.consensus_name.lower() or "p53" in res.consensus_name.lower()
    assert res.confidence_score >= 0.85
    assert res.cath_summary is not None
    assert res.cath_summary.topology_name == "Immunoglobulin-like"
