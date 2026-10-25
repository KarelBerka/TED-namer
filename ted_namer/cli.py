"""Command-line interface for TED Domain Cluster & Family Namer (ted-namer)."""

from __future__ import annotations

import asyncio
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from ted_namer.export import save_output
from ted_namer.models import ClusterNamingResult
from ted_namer.pipeline import TedNamingPipeline, load_clusters_from_file

app = typer.Typer(
    name="ted-namer",
    help="Generate human-readable names for TED (ted.cathdb.info) structural domain clusters and families.",
    add_completion=False,
)
console = Console()


def _render_cluster_results_to_console(results: list[ClusterNamingResult]) -> None:
    summary_table = Table(
        title="TED Domain Cluster & Family Naming Summary",
        show_lines=True,
        header_style="bold cyan",
    )
    summary_table.add_column("Cluster ID", style="bold")
    summary_table.add_column("Members", justify="right")
    summary_table.add_column("Consensus Name", style="bold green")
    summary_table.add_column("Short Label", style="yellow")
    summary_table.add_column("Tier", style="magenta")
    summary_table.add_column("CATH Label", style="cyan")
    summary_table.add_column("Conf", justify="right")
    summary_table.add_column("Top Supporting Evidence")

    for r in results:
        top = r.top_candidates[0] if r.top_candidates else None
        ev_str = (
            f"{','.join(top.source_databases)} ({','.join(top.supporting_accessions[:3])}) "
            f"IoU={top.mean_iou:.2f}"
            if top
            else r.rationale
        )
        summary_table.add_row(
            r.cluster_id,
            str(r.member_count),
            r.consensus_name,
            r.short_name,
            r.specificity_tier.value,
            r.dominant_cath_label,
            f"{r.confidence_score:.2f}",
            ev_str,
        )

    console.print(summary_table)


@app.command("name")
def name_command(
    input_file: Optional[Path] = typer.Option(
        None,
        "--input",
        "-i",
        help="Path to TSV, CSV, JSON, or text file defining TED domain clusters.",
    ),
    ted_ids: Optional[str] = typer.Option(
        None,
        "--ted-ids",
        "-t",
        help="Comma-separated TED IDs or UniProt accessions (e.g. 'AF-P04637-F1-model_v4_TED01,AF-Q9H3D4-F1-model_v4_TED01').",
    ),
    cluster_id: str = typer.Option(
        "cluster_1",
        "--cluster-id",
        "-c",
        help="Cluster ID when passing --ted-ids with --group-by cluster_id.",
    ),
    group_by: str = typer.Option(
        "cluster_id",
        "--group-by",
        "-g",
        help="Grouping mode: 'cluster_id' (keep input clusters), 'cath_label' (group by CATH code), or 'single'.",
    ),
    output: Optional[list[Path]] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file path(s) (.json, .tsv, .csv, .md). Can be specified multiple times.",
    ),
    no_llm: bool = typer.Option(
        False,
        "--no-llm",
        help="Disable LLM consensus refinement and use deterministic rule-based naming only.",
    ),
    llm_provider: str = typer.Option(
        "auto",
        "--llm-provider",
        help="LLM provider ('auto', 'openai', 'gemini', 'none').",
    ),
    llm_model: Optional[str] = typer.Option(
        None,
        "--llm-model",
        help="Optional LLM model name override.",
    ),
    llm_base_url: Optional[str] = typer.Option(
        None,
        "--llm-base-url",
        help="Optional OpenAI-compatible base URL (e.g. for local Ollama / vLLM).",
    ),
    max_members: int = typer.Option(
        50,
        "--max-members",
        help="Maximum number of domain members to sample per cluster.",
    ),
    concurrency: int = typer.Option(
        6,
        "--concurrency",
        help="Maximum concurrent API requests.",
    ),
    no_cache: bool = typer.Option(
        False,
        "--no-cache",
        help="Disable local SQLite API response caching.",
    ),
) -> None:
    """Generate consensus human-readable names for TED domain clusters/families."""
    cluster_map: dict[str, list[str]] = {}
    if input_file:
        if not input_file.exists():
            console.print(f"[bold red]Error:[/bold red] Input file not found: {input_file}")
            raise typer.Exit(code=1)
        cluster_map.update(load_clusters_from_file(input_file))

    if ted_ids:
        parsed_ids = [x.strip() for x in ted_ids.replace(";", ",").split(",") if x.strip()]
        if parsed_ids:
            cluster_map.setdefault(cluster_id, []).extend(parsed_ids)

    if not cluster_map:
        console.print("[bold red]Error:[/bold red] Provide either --input <file> or --ted-ids <id1,id2,...>.")
        raise typer.Exit(code=1)

    pipeline = TedNamingPipeline(
        use_cache=not no_cache,
        concurrency=concurrency,
        use_llm=not no_llm,
        llm_provider=llm_provider,
        llm_model=llm_model,
        llm_base_url=llm_base_url,
    )

    with console.status("[bold cyan]Fetching TED, InterPro, CATH, and UniProt annotations...[/bold cyan]"):
        results = asyncio.run(
            pipeline.name_clusters(
                cluster_map=cluster_map,
                group_by=group_by,
                max_members_per_cluster=max_members,
            )
        )

    _render_cluster_results_to_console(results)

    if output:
        for out_path in output:
            save_output(results, out_path)
            console.print(f"[bold green]Saved output to:[/bold green] {out_path}")


@app.command("inspect")
def inspect_command(
    identifier: str = typer.Argument(
        ...,
        help="UniProt accession (e.g. 'P04637') or TED domain ID (e.g. 'AF-P04637-F1-model_v4_TED01').",
    ),
    no_llm: bool = typer.Option(
        False,
        "--no-llm",
        help="Disable LLM consensus refinement.",
    ),
    no_cache: bool = typer.Option(
        False,
        "--no-cache",
        help="Disable local SQLite API response caching.",
    ),
) -> None:
    """Inspect all TED domains on a UniProt protein or a specific TED domain with detailed boundary overlaps."""
    pipeline = TedNamingPipeline(
        use_cache=not no_cache,
        use_llm=not no_llm,
    )

    with console.status(f"[bold cyan]Inspecting {identifier} on ted.cathdb.info, InterPro, CATH, and UniProt...[/bold cyan]"):
        # Group each TED domain into its own cluster so each domain on the protein gets named individually
        initial_results = asyncio.run(
            pipeline.name_clusters(
                cluster_map={"inspect": [identifier]},
                group_by="single",
            )
        )

    if not initial_results or not initial_results[0].members:
        console.print(f"[bold yellow]No TED domains found for identifier:[/bold yellow] {identifier}")
        raise typer.Exit(code=1)

    members = initial_results[0].members
    # Name each domain individually to show its specific resolved name
    from ted_namer.aggregator import aggregate_cluster_evidence

    for m in members:
        single_res = aggregate_cluster_evidence(cluster_id=m.domain.short_ted_num, members=[m])
        prot_title = (
            f"{m.protein_context.recommended_name} ({m.domain.uniprot_acc}, {m.protein_context.organism})"
            if m.protein_context and m.protein_context.recommended_name
            else m.domain.uniprot_acc
        )
        cath_desc = m.domain.cath_label
        if m.cath_info:
            parts = [p for p in (m.cath_info.superfamily_name, m.cath_info.topology_name) if p]
            if parts:
                cath_desc = f"{m.domain.cath_label} ({' / '.join(parts)})"

        header_text = (
            f"[bold green]Proposed Name:[/bold green] [bold white]{single_res.consensus_name}[/bold white] "
            f"([yellow]{single_res.short_name}[/yellow], conf={single_res.confidence_score:.2f}, tier={single_res.specificity_tier.value})\n"
            f"[bold cyan]Parent Protein:[/bold cyan] {prot_title}\n"
            f"[bold cyan]Chopping:[/bold cyan] {m.domain.chopping} ({m.domain.nres_domain} res, {m.relative_position}, pLDDT={m.domain.plddt:.1f})\n"
            f"[bold cyan]CATH Classification:[/bold cyan] {cath_desc}\n"
            f"[dim]{single_res.rationale}[/dim]"
        )

        annot_table = Table(show_header=True, header_style="bold magenta", expand=True)
        annot_table.add_column("Database")
        annot_table.add_column("Accession", style="cyan")
        annot_table.add_column("Annotation Name", style="bold")
        annot_table.add_column("Type")
        annot_table.add_column("Chopping")
        annot_table.add_column("IoU", justify="right")
        annot_table.add_column("TED Cov", justify="right")
        annot_table.add_column("Match Category")

        for a in m.overlapping_annotations[:10]:
            annot_table.add_row(
                a.source_db,
                a.accession,
                a.name,
                a.annotation_type,
                a.chopping,
                f"{a.iou:.2f}",
                f"{a.ted_coverage:.2f}",
                a.overlap_category.value,
            )

        console.print(Panel(header_text, title=f"[bold]{m.domain.ted_id}[/bold]", border_style="green"))
        if m.overlapping_annotations:
            console.print(annot_table)
        else:
            console.print("[italic dim]No overlapping sequence annotations found; named via CATH/UniProt context.[/italic dim]")
        console.print()


if __name__ == "__main__":
    app()
