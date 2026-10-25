"""Export utilities for formatting ClusterNamingResult lists as JSON, TSV/CSV, and Markdown reports."""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path

from ted_namer.models import ClusterNamingResult


def to_json(results: list[ClusterNamingResult], indent: int = 2) -> str:
    """Serialize ClusterNamingResult list to formatted JSON."""
    payload = [r.model_dump() for r in results]
    return json.dumps(payload, indent=indent)


def to_delimited(results: list[ClusterNamingResult], delimiter: str = "\t") -> str:
    """Serialize ClusterNamingResult list to a tabular TSV or CSV string."""
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter=delimiter, lineterminator="\n")
    writer.writerow(
        [
            "cluster_id",
            "member_count",
            "consensus_name",
            "short_name",
            "confidence_score",
            "specificity_tier",
            "dominant_cath_label",
            "cath_topology_name",
            "consensus_relative_position",
            "supporting_databases",
            "supporting_accessions",
            "member_ted_ids",
            "rationale",
        ]
    )
    for r in results:
        top = r.top_candidates[0] if r.top_candidates else None
        dbs = ",".join(top.source_databases) if top else ""
        accs = ",".join(top.supporting_accessions[:5]) if top else ""
        topo_name = r.cath_summary.topology_name if r.cath_summary and r.cath_summary.topology_name else ""
        member_ids = ",".join(m.domain.ted_id for m in r.members[:20])
        writer.writerow(
            [
                r.cluster_id,
                r.member_count,
                r.consensus_name,
                r.short_name,
                f"{r.confidence_score:.3f}",
                r.specificity_tier.value,
                r.dominant_cath_label,
                topo_name,
                r.consensus_relative_position or "",
                dbs,
                accs,
                member_ids,
                r.rationale,
            ]
        )
    return buf.getvalue()


def to_markdown_report(results: list[ClusterNamingResult]) -> str:
    """Generate a rich Markdown report with candidate tables, member details, and ted.cathdb.info links."""
    lines: list[str] = [
        "# TED Structural Domain Cluster Naming Report",
        "",
        "| Cluster ID | Members | Consensus Domain Name | Short Label | Tier | CATH Label | Confidence |",
        "| :--- | :---: | :--- | :--- | :--- | :---: | :---: |",
    ]
    for r in results:
        lines.append(
            f"| `{r.cluster_id}` | {r.member_count} | **{r.consensus_name}** | `{r.short_name}` | "
            f"`{r.specificity_tier.value}` | `{r.dominant_cath_label}` | {r.confidence_score:.2f} |"
        )

    lines.append("")
    for r in results:
        lines.extend(
            [
                f"## Cluster `{r.cluster_id}`: {r.consensus_name} (`{r.short_name}`)",
                "",
                f"- **Confidence Score**: `{r.confidence_score:.3f}` (`{r.naming_method}`)",
                f"- **Specificity Tier**: `{r.specificity_tier.value}`",
                f"- **Dominant CATH Classification**: `{r.dominant_cath_label}`"
                + (
                    f" — Superfamily: *{r.cath_summary.superfamily_name or 'N/A'}*, "
                    f"Topology: *{r.cath_summary.topology_name or 'N/A'}*"
                    if r.cath_summary
                    else ""
                ),
                f"- **Consensus Relative Position**: `{r.consensus_relative_position}`",
                f"- **Rationale**: {r.rationale}",
                "",
            ]
        )
        if r.top_candidates:
            lines.extend(
                [
                    "### Top Ranked Candidate Annotations",
                    "",
                    "| Rank | Candidate Name | Tier | Freq | Mean IoU | TED Cov | Databases | Accessions | Score |",
                    "| :---: | :--- | :--- | :---: | :---: | :---: | :--- | :--- | :---: |",
                ]
            )
            for idx, c in enumerate(r.top_candidates[:6], start=1):
                dbs = ", ".join(c.source_databases)
                accs = ", ".join(c.supporting_accessions[:4])
                lines.append(
                    f"| {idx} | {c.canonical_name} | `{c.specificity_tier.value}` | "
                    f"{c.cluster_frequency:.0%} ({c.member_hits}/{r.member_count}) | "
                    f"{c.mean_iou:.2f} | {c.mean_ted_coverage:.2f} | {dbs} | `{accs}` | **{c.composite_score:.3f}** |"
                )
            lines.append("")

        if r.members:
            lines.extend(
                [
                    "### Cluster Member TED Domains",
                    "",
                    "| TED Domain ID | UniProt | Protein Name | Chopping | Position | pLDDT | CATH Label | Top Overlapping Hit (IoU) |",
                    "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- |",
                ]
            )
            for m in r.members:
                uni = m.domain.uniprot_acc
                ted_link = f"[{m.domain.ted_id}](https://ted.cathdb.info/uniprot/{uni})"
                prot_name = (
                    m.protein_context.recommended_name
                    if m.protein_context and m.protein_context.recommended_name
                    else "-"
                )
                top_hit = (
                    f"{m.overlapping_annotations[0].name} (`{m.overlapping_annotations[0].accession}`, IoU={m.overlapping_annotations[0].iou:.2f})"
                    if m.overlapping_annotations
                    else "None"
                )
                eff_cath = m.domain.inferred_cath_label or m.domain.cath_label
                lines.append(
                    f"| {ted_link} | `{uni}` | {prot_name} | `{m.domain.chopping}` | "
                    f"`{m.relative_position}` | {m.domain.plddt:.1f} | `{eff_cath}` | {top_hit} |"
                )
            lines.append("")

    return "\n".join(lines)


def save_output(results: list[ClusterNamingResult], output_path: str | Path) -> None:
    """Write ClusterNamingResult list to disk based on file extension (.json, .tsv, .csv, .md)."""
    path = Path(output_path)
    ext = path.suffix.lower()
    if ext == ".json":
        content = to_json(results)
    elif ext == ".csv":
        content = to_delimited(results, delimiter=",")
    elif ext in {".md", ".markdown"}:
        content = to_markdown_report(results)
    else:
        content = to_delimited(results, delimiter="\t")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
