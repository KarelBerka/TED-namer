"""Unit tests for export formats (TSV/CSV/JSON/MD), Typer CLI commands, and LLM evidence dossiers."""

from __future__ import annotations

import csv
import io
import json

from typer.testing import CliRunner

from ted_namer.cli import app
from ted_namer.export import save_output, to_delimited, to_json, to_markdown_report
from ted_namer.llm import _extract_json, _resolve_llm_config, build_evidence_dossier
from ted_namer.models import (
    CathResolution,
    ClusterCandidateTerm,
    ClusterNamingResult,
    EnrichedTedDomain,
    ProteinContext,
    SpecificityTier,
    TedDomainMember,
)


def _sample_result() -> ClusterNamingResult:
    dom = TedDomainMember(
        ted_id="AF-Q8U3D2-F1-model_v4_TED01",
        uniprot_acc="Q8U3D2",
        chopping="1-146_284-319",
        nres_domain=182,
        plddt=94.9,
        norm_rg=0.32,
        consensus_level="high",
        cath_label="-",
        inferred_cath_label="3.?.?.?",
        cath_fold_hint="3.20.20.10",
        novel_fold_candidate=True,
        cath_estimation_source="rule_3d_tim_barrel_3.20.20",
    )
    cath = CathResolution(
        cath_label="3.?.?.?",
        numeric_cath_code="3",
        raw_ted_cath_label="-",
        is_estimated=True,
        estimation_source="rule_3d_tim_barrel_3.20.20",
        class_name="Alpha Beta",
        fold_hint="3.20.20.10",
        novel_fold_candidate=True,
    )
    cand = ClusterCandidateTerm(
        canonical_name="Argonaute N-terminal domain",
        short_name="Argonaute_N-term_dom",
        specificity_tier=SpecificityTier.SPECIFIC_DOMAIN,
        source_databases=["uniprot_feature"],
        supporting_accessions=["UP_FEAT:Q8U3D2:1-146"],
        member_hits=1,
        cluster_frequency=1.0,
        mean_iou=0.82,
        mean_ted_coverage=0.95,
        mean_annot_coverage=0.84,
        reviewed_hit_count=1,
        composite_score=0.892,
    )
    member = EnrichedTedDomain(
        domain=dom,
        cath_info=cath,
        protein_context=ProteinContext(
            uniprot_acc="Q8U3D2",
            recommended_name="Protein argonaute",
            is_reviewed=True,
            protein_length=770,
        ),
        relative_position="N-terminal",
    )
    return ClusterNamingResult(
        cluster_id="test_cluster_1",
        member_count=1,
        consensus_name="Argonaute N-terminal domain",
        short_name="Argonaute_N-term_dom",
        confidence_score=0.892,
        specificity_tier=SpecificityTier.SPECIFIC_DOMAIN,
        naming_method="rule_based",
        rationale="Selected 'Argonaute N-terminal domain' via curated UniProt feature overlap.",
        dominant_cath_label="3.?.?.?",
        cath_fold_hint="3.20.20.10",
        novel_fold_candidate=True,
        structure_quality_flags=["ok"],
        cath_summary=cath,
        consensus_relative_position="N-terminal",
        top_candidates=[cand],
        members=[member],
    )


def test_export_formats_and_save_output(tmp_path) -> None:
    res = _sample_result()

    # TSV & CSV delimited output
    for delimiter, ext in (("\t", "tsv"), (",", "csv")):
        text = to_delimited([res], delimiter=delimiter)
        rows = list(csv.DictReader(io.StringIO(text), delimiter=delimiter))
        assert len(rows) == 1
        row = rows[0]
        assert row["cluster_id"] == "test_cluster_1"
        assert row["consensus_name"] == "Argonaute N-terminal domain"
        assert row["dominant_cath_label"] == "3.?.?.?"
        assert row["cath_fold_hint"] == "3.20.20.10"
        assert row["novel_fold_candidate"] == "true"
        assert row["structure_quality_flags"] == "ok"

        out_file = tmp_path / f"out.{ext}"
        save_output([res], out_file)
        assert out_file.exists() and "Argonaute N-terminal domain" in out_file.read_text(encoding="utf-8")

    # JSON export
    json_str = to_json([res])
    loaded = json.loads(json_str)
    assert loaded[0]["cath_fold_hint"] == "3.20.20.10"
    assert loaded[0]["novel_fold_candidate"] is True
    assert loaded[0]["structure_quality_flags"] == ["ok"]
    json_file = tmp_path / "out.json"
    save_output([res], json_file)
    assert json_file.exists()

    # Markdown report export
    md = to_markdown_report([res])
    assert "Novel Fold Candidate" in md
    assert "Structure Quality" in md
    assert "3.20.20.10" in md
    md_file = tmp_path / "out.md"
    save_output([res], md_file)
    assert md_file.exists()


def test_cli_help_and_missing_input() -> None:
    runner = CliRunner()
    help_res = runner.invoke(app, ["--help"])
    assert help_res.exit_code == 0
    assert "name" in help_res.stdout
    assert "inspect" in help_res.stdout
    assert "batch" not in help_res.stdout

    err_res = runner.invoke(app, ["name"])
    assert err_res.exit_code == 1
    assert "Provide either --input" in err_res.stdout


def test_llm_dossier_config_and_json_extraction(monkeypatch) -> None:
    res = _sample_result()
    dossier = build_evidence_dossier(res)
    assert dossier["dominant_cath_label"] == "3.?.?.?"
    assert dossier["cath_fold_hint"] == "3.20.20.10"
    assert dossier["novel_fold_candidate"] is True
    assert dossier["structure_quality_flags"] == ["ok"]

    # _resolve_llm_config without env keys
    monkeypatch.delenv("GEMINI_API_KEY", raising=False)
    monkeypatch.delenv("GOOGLE_API_KEY", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    assert _resolve_llm_config("none", None, None, None) == (None, None, None, None)
    assert _resolve_llm_config("auto", None, None, None) == (None, None, None, None)
    assert _resolve_llm_config("openai", "sk-test", None, "gpt-4o")[0] == "openai"
    assert _resolve_llm_config("gemini", "gm-test", None, "gemini-2.5-flash")[0] == "gemini"

    # _extract_json strips fenced code blocks
    raw = '```json\n{"consensus_name": "p53 DNA-binding domain", "confidence_score": 0.95}\n```'
    parsed = _extract_json(raw)
    assert parsed["consensus_name"] == "p53 DNA-binding domain"
