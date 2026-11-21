"""Compute protein-clustered & domain-level bootstrap 95% CIs and render paper/figure1_validation.png."""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter, defaultdict
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ted_namer.cath_classifier import _G3DSA_RE, classify_cath_with_rules
from ted_namer.models import TedDomainMember

ROOT = Path(__file__).resolve().parent.parent
JSONL_PATH = ROOT / "examples" / "random_100000_compact.jsonl"
SUMMARY_PATH = ROOT / "examples" / "benchmark_100000_summary.json"
FIG_PATH = ROOT / "paper" / "figure1_validation.png"

N_BOOT = 2000
SEED = 42

EXPECTED_PANEL_A = {
    "tier_1_specific_pct": 62.425,
    "tier_2_superfamily_pct": 24.537,
    "tier_1_plus_2_pct": 86.962,
    "tier_3_family_pct": 12.858,
    "tier_4_5_fallback_pct": 0.180,
}

EXPECTED_PANEL_B = {
    "direct_precision": [98.16, 98.57, 98.07, 99.78],
    "direct_claimed": [72130, 67115, 67110, 48913],
    "profile_cv_accuracy": [99.28, 98.63, 97.94, 92.29],
    "full_cascade_accuracy": [98.16, 95.22, 93.94, 91.10],
}


def _extract_g3dsa_parts(ann: dict) -> list[str] | None:
    """Extract 4-level CATH components from a CATH-Gene3D annotation accession."""
    match = _G3DSA_RE.match(str(ann.get("acc", "")).strip())
    if match:
        return list(match.groups())
    return None


def _build_fold_maps(
    official_entries: list[tuple[list[str], dict]],
    n_folds: int = 5,
) -> list[dict[str, tuple[list[str], int]]]:
    """Build blind hold-out 5-fold profile-to-CATH maps using the exact benchmark rules."""
    fold_maps: list[dict[str, tuple[list[str], int]]] = []
    for fold_idx in range(n_folds):
        acc_to_cath: dict[str, Counter[str]] = defaultdict(Counter)
        for i, (tp, rec) in enumerate(official_entries):
            if i % n_folds == fold_idx:
                continue
            cath_str = ".".join(tp)
            for ann in rec.get("annots", []):
                if ann.get("db") == "cathgene3d" or str(ann.get("acc", "")).upper().startswith("G3DSA:"):
                    continue
                if not (ann["iou"] >= 0.40 or (ann["tcov"] >= 0.65 and ann["acov"] >= 0.45)):
                    continue
                for key in (ann.get("acc"), ann.get("ipr")):
                    if key:
                        acc_to_cath[str(key)][cath_str] += 1

        fold_map: dict[str, tuple[list[str], int]] = {}
        for key, counts in acc_to_cath.items():
            tot = sum(counts.values())
            if tot == 0:
                continue
            top_cath, top_cath_cnt = counts.most_common(1)[0]
            h_purity = top_cath_cnt / tot
            cat_counts: Counter[str] = Counter()
            for cath_code, cnt in counts.items():
                cat_counts[".".join(cath_code.split(".")[:3])] += cnt
            top_cat, top_cat_cnt = cat_counts.most_common(1)[0]
            cat_purity = top_cat_cnt / tot
            if cat_purity < 0.85:
                continue
            k_depth = 3 if (tot >= 2 and cat_purity >= 0.90) else 1
            code = top_cath if h_purity >= 0.85 else f"{top_cat}.10"
            fold_map[key] = (code.split(".")[:4], k_depth)
        fold_maps.append(fold_map)
    return fold_maps


def _bootstrap_ratios(
    cluster_mat: np.ndarray,
    num_cols: list[int],
    den_cols: list[int],
    n_boot: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Run vector-indexed nonparametric cluster bootstrap for exact ratio estimators."""
    rng = np.random.default_rng(seed)
    n_clusters = cluster_mat.shape[0]
    k = len(num_cols)
    num_arr = np.asarray(num_cols, dtype=np.intp)
    den_arr = np.asarray(den_cols, dtype=np.intp)

    # Point estimate from total counts
    totals = cluster_mat.sum(axis=0, dtype=np.float64)
    point = 100.0 * totals[num_arr] / totals[den_arr]

    boot = np.empty((n_boot, k), dtype=np.float64)
    for b in range(n_boot):
        idx = rng.integers(0, n_clusters, size=n_clusters)
        sampled = cluster_mat[idx].sum(axis=0, dtype=np.float64)
        boot[b] = 100.0 * sampled[num_arr] / sampled[den_arr]
    return point, boot


def _summarize_ci(point: float, samples: np.ndarray) -> dict[str, float]:
    """Compute percentile 95% CI, bootstrap SE, and half-widths."""
    low, high = np.percentile(samples, [2.5, 97.5])
    se = float(np.std(samples, ddof=1))
    hw_sym = max(point - float(low), float(high) - point)
    return {
        "point_pct": round(float(point), 4),
        "ci95_low": round(float(low), 4),
        "ci95_high": round(float(high), 4),
        "halfwidth_sym": round(float(hw_sym), 4),
        "se_pct": round(se, 4),
    }


def load_and_bootstrap() -> dict:
    """Load the 100k dataset, verify exact point parity, and compute bootstrap 95% CIs."""
    t0 = time.perf_counter()
    if not JSONL_PATH.exists():
        raise FileNotFoundError(f"Missing input dataset: {JSONL_PATH}")

    records: list[dict] = []
    with JSONL_PATH.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                records.append(json.loads(line))

    if len(records) != 100_000:
        raise AssertionError(f"Expected 100,000 records, got {len(records)}")

    # Panel A counts per domain and per UniProt cluster:
    # cols: [total, t1, t2, t1_2, t3, t4_5]
    prot_a_idx: dict[str, int] = {}
    prot_a_rows: list[list[int]] = []
    dom_a_mat = np.zeros((len(records), 6), dtype=np.int32)

    official_entries: list[tuple[list[str], dict]] = []
    for i, rec in enumerate(records):
        tier = rec["specificity_tier"]
        t1 = int(tier == "specific_domain")
        t2 = int(tier == "superfamily")
        t3 = int(tier == "protein_family_context")
        t45 = int(tier in ("topology_fold", "novel_uncharacterized"))
        row = [1, t1, t2, t1 + t2, t3, t45]
        dom_a_mat[i] = row

        acc = rec["uniprot_acc"]
        if acc not in prot_a_idx:
            prot_a_idx[acc] = len(prot_a_rows)
            prot_a_rows.append([0, 0, 0, 0, 0, 0])
        p_row = prot_a_rows[prot_a_idx[acc]]
        for c in range(6):
            p_row[c] += row[c]

        cath_raw = (rec.get("cath_label") or "-").split(",")[0].strip()
        parts = [p.strip() for p in cath_raw.split(".") if p.strip() and p.strip() != "-"]
        if len(parts) >= 4:
            official_entries.append((parts[:4], rec))

    if len(official_entries) != 72_130:
        raise AssertionError(f"Expected 72,130 official 4-level CATH domains, got {len(official_entries)}")

    prot_a_mat = np.asarray(prot_a_rows, dtype=np.int32)
    print(f"Loaded 100,000 domains across {len(prot_a_mat):,} UniProt proteins ({time.perf_counter() - t0:.1f}s)")

    # Build 5-fold CV profile maps for Panel B
    fold_maps = _build_fold_maps(official_entries, n_folds=5)

    # Panel B matrix (18 columns):
    # [0]: off_tot (1)
    # [1]: prof_tot
    # [2..5]: cas_ok C..H
    # [6..9]: prof_ok C..H
    # [10..13]: direct_claimed C..H
    # [14..17]: direct_ok C..H
    prot_b_idx: dict[str, int] = {}
    prot_b_rows: list[list[int]] = []
    dom_b_mat = np.zeros((len(official_entries), 18), dtype=np.int32)

    for i, (tp, rec) in enumerate(official_entries):
        f_map = fold_maps[i % 5]
        best_g3d: tuple[float, list[str], float, float] | None = None
        best_prof: tuple[float, list[str], int, float] | None = None

        for ann in rec.get("annots", []):
            is_g3d = ann.get("db") == "cathgene3d" or str(ann.get("acc", "")).upper().startswith("G3DSA:")
            if is_g3d:
                if ann["iou"] >= 0.25 or ann["tcov"] >= 0.45:
                    g_parts = _extract_g3dsa_parts(ann)
                    if g_parts:
                        score = ann["iou"] + 0.5 * ann["tcov"]
                        if best_g3d is None or score > best_g3d[0]:
                            best_g3d = (score, g_parts, ann["iou"], ann["tcov"])
            else:
                if ann["iou"] >= 0.40 or (ann["tcov"] >= 0.65 and ann["acov"] >= 0.45):
                    for key in (ann.get("acc"), ann.get("ipr")):
                        if key and str(key) in f_map:
                            p_parts, p_depth = f_map[str(key)]
                            score = (
                                (ann["iou"] + 0.3 * ann["tcov"])
                                * (1.15 if ann.get("db") == "ssf" else 1.0)
                                * (1.1 if p_depth >= 3 else 1.0)
                            )
                            if best_prof is None or score > best_prof[0]:
                                best_prof = (score, p_parts, p_depth, ann["iou"])

        dom_obj = TedDomainMember(
            ted_id=rec["ted_id"],
            uniprot_acc=rec["uniprot_acc"],
            chopping=rec["chopping"],
            nres_domain=rec["nres_domain"],
            num_helix=rec.get("num_helix", 0),
            num_strand=rec.get("num_strand", 0),
            num_turn=rec.get("num_turn", 0),
            norm_rg=rec.get("norm_rg", 0.3),
            packing_density=rec.get("packing_density", 10.0),
            plddt=rec.get("plddt", 85.0),
            cath_label="-",
        )
        rule_parts, rule_depth, _ = classify_cath_with_rules(dom_obj, pdb_text=None)

        row = [0] * 18
        row[0] = 1
        if best_prof is not None:
            row[1] = 1
            p_parts = best_prof[1]
            for d in range(4):
                if p_parts[: d + 1] == tp[: d + 1]:
                    row[6 + d] = 1

        if best_g3d is not None:
            pred_a = best_g3d[1]
            k_depth = 4 if (best_g3d[2] >= 0.65 and best_g3d[3] >= 0.75) else (
                3 if (best_g3d[2] >= 0.40 or best_g3d[3] >= 0.65) else 1
            )
        elif best_prof is not None:
            pred_a = best_prof[1]
            k_depth = best_prof[2] if best_prof[3] >= 0.45 else 1
        else:
            pred_a = rule_parts
            k_depth = rule_depth

        for d in range(4):
            match_d = int(pred_a[: d + 1] == tp[: d + 1])
            row[2 + d] = match_d
            if d + 1 <= k_depth:
                row[10 + d] = 1
                row[14 + d] = match_d

        dom_b_mat[i] = row
        acc = rec["uniprot_acc"]
        if acc not in prot_b_idx:
            prot_b_idx[acc] = len(prot_b_rows)
            prot_b_rows.append([0] * 18)
        p_row = prot_b_rows[prot_b_idx[acc]]
        for c in range(18):
            p_row[c] += row[c]

    prot_b_mat = np.asarray(prot_b_rows, dtype=np.int32)
    print(f"Evaluated 72,130 official CATH domains across {len(prot_b_mat):,} UniProt proteins ({time.perf_counter() - t0:.1f}s)")

    # Run bootstrap (protein-clustered primary + domain-level comparison)
    a_num_cols = [1, 2, 3, 4, 5]
    a_den_cols = [0, 0, 0, 0, 0]
    a_pt, a_boot_clust = _bootstrap_ratios(prot_a_mat, a_num_cols, a_den_cols, N_BOOT, SEED)
    _, a_boot_dom = _bootstrap_ratios(dom_a_mat, a_num_cols, a_den_cols, N_BOOT, SEED + 1)

    b_num_cols = [14, 15, 16, 17, 6, 7, 8, 9, 2, 3, 4, 5, 1]
    b_den_cols = [10, 11, 12, 13, 1, 1, 1, 1, 0, 0, 0, 0, 0]
    b_pt, b_boot_clust = _bootstrap_ratios(prot_b_mat, b_num_cols, b_den_cols, N_BOOT, SEED + 2)
    _, b_boot_dom = _bootstrap_ratios(dom_b_mat, b_num_cols, b_den_cols, N_BOOT, SEED + 3)

    # Assert exact point-estimate parity with benchmark_100000_summary.json
    a_keys = [
        "tier_1_specific_pct",
        "tier_2_superfamily_pct",
        "tier_1_plus_2_pct",
        "tier_3_family_pct",
        "tier_4_5_fallback_pct",
    ]
    for idx, k in enumerate(a_keys):
        diff = abs(a_pt[idx] - EXPECTED_PANEL_A[k])
        if diff > 1e-2:
            raise AssertionError(f"Panel A parity mismatch on {k}: got {a_pt[idx]:.4f} vs expected {EXPECTED_PANEL_A[k]:.4f}")

    b_totals = dom_b_mat.sum(axis=0)
    claimed = [int(b_totals[10 + d]) for d in range(4)]
    if claimed != EXPECTED_PANEL_B["direct_claimed"]:
        raise AssertionError(f"Panel B claimed mismatch: got {claimed} vs expected {EXPECTED_PANEL_B['direct_claimed']}")

    for d in range(4):
        if abs(b_pt[d] - EXPECTED_PANEL_B["direct_precision"][d]) > 0.015:
            raise AssertionError(f"Direct precision parity mismatch at level {d+1}: {b_pt[d]:.4f} vs {EXPECTED_PANEL_B['direct_precision'][d]}")
        if abs(b_pt[4 + d] - EXPECTED_PANEL_B["profile_cv_accuracy"][d]) > 0.015:
            raise AssertionError(f"Profile CV parity mismatch at level {d+1}: {b_pt[4+d]:.4f} vs {EXPECTED_PANEL_B['profile_cv_accuracy'][d]}")
        if abs(b_pt[8 + d] - EXPECTED_PANEL_B["full_cascade_accuracy"][d]) > 0.015:
            raise AssertionError(f"Full cascade parity mismatch at level {d+1}: {b_pt[8+d]:.4f} vs {EXPECTED_PANEL_B['full_cascade_accuracy'][d]}")

    levels = ["C", "C.A", "C.A.T", "C.A.T.H"]
    results = {
        "method": "Nonparametric bootstrap (2.5th-97.5th percentiles)",
        "replicates": N_BOOT,
        "seed": SEED,
        "primary_estimator": "protein_clustered_by_uniprot_acc",
        "n_domains_all": int(len(records)),
        "n_uniprot_clusters_all": int(len(prot_a_mat)),
        "n_domains_official_cath": int(len(official_entries)),
        "n_uniprot_clusters_official_cath": int(len(prot_b_mat)),
        "panel_a_specificity_tiers": {
            "protein_clustered": {k: _summarize_ci(a_pt[idx], a_boot_clust[:, idx]) for idx, k in enumerate(a_keys)},
            "domain_iid_comparison": {k: _summarize_ci(a_pt[idx], a_boot_dom[:, idx]) for idx, k in enumerate(a_keys)},
        },
        "panel_b_cath_validation": {
            "confirmed_direct_precision": {
                "protein_clustered": {lvl: _summarize_ci(b_pt[d], b_boot_clust[:, d]) for d, lvl in enumerate(levels)},
                "domain_iid_comparison": {lvl: _summarize_ci(b_pt[d], b_boot_dom[:, d]) for d, lvl in enumerate(levels)},
            },
            "profile_to_cath_5fold_cv": {
                "coverage_pct": _summarize_ci(b_pt[12], b_boot_clust[:, 12]),
                "protein_clustered": {lvl: _summarize_ci(b_pt[4 + d], b_boot_clust[:, 4 + d]) for d, lvl in enumerate(levels)},
                "domain_iid_comparison": {lvl: _summarize_ci(b_pt[4 + d], b_boot_dom[:, 4 + d]) for d, lvl in enumerate(levels)},
            },
            "full_fallback_cascade": {
                "protein_clustered": {lvl: _summarize_ci(b_pt[8 + d], b_boot_clust[:, 8 + d]) for d, lvl in enumerate(levels)},
                "domain_iid_comparison": {lvl: _summarize_ci(b_pt[8 + d], b_boot_dom[:, 8 + d]) for d, lvl in enumerate(levels)},
            },
        },
    }
    print(f"Bootstrap completed in {time.perf_counter() - t0:.1f}s")
    return results


def _fmt2(val: float) -> str:
    """Format a float with standard half-up 2-decimal rounding."""
    return str(Decimal(str(val)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))


def render_figure(results: dict) -> None:
    """Render two-panel validation figure with 95% protein-clustered bootstrap error bars."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "axes.linewidth": 0.9,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
        }
    )

    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(13.6, 5.1),
        gridspec_kw={"width_ratios": [1.0, 1.08], "wspace": 0.23},
    )

    # --- Panel A: Naming Specificity Tiers ---
    a_stats = results["panel_a_specificity_tiers"]["protein_clustered"]
    tier_order = [
        ("tier_1_specific_pct", "Tier 1: Specific Domain\n(& Multi-Lobe Subdomain)", "#1f77b4"),
        ("tier_2_superfamily_pct", "Tier 2: Homologous\nSuperfamily (Curated)", "#2ca02c"),
        ("tier_3_family_pct", "Tier 3: Protein Family\n+ Positional Context", "#ff7f0e"),
        ("tier_4_5_fallback_pct", "Tier 4/5: Fold Fallback\n& Uncharacterized", "#7f7f7f"),
    ]
    labels_a = [item[1] for item in tier_order]
    colors_a = [item[2] for item in tier_order]
    vals_a = np.array([a_stats[item[0]]["point_pct"] for item in tier_order])
    lows_a = np.array([vals_a[i] - a_stats[item[0]]["ci95_low"] for i, item in enumerate(tier_order)])
    highs_a = np.array([a_stats[item[0]]["ci95_high"] - vals_a[i] for i, item in enumerate(tier_order)])
    hws_a = np.array([a_stats[item[0]]["halfwidth_sym"] for item in tier_order])

    y_pos = np.arange(len(tier_order))
    ax1.barh(
        y_pos,
        vals_a,
        height=0.58,
        color=colors_a,
        edgecolor="#222222",
        linewidth=0.85,
        zorder=3,
    )
    ax1.errorbar(
        vals_a,
        y_pos,
        xerr=[lows_a, highs_a],
        fmt="none",
        ecolor="#111111",
        elinewidth=1.4,
        capsize=3.5,
        capthick=1.1,
        zorder=4,
    )

    for i, (val, high_err, hw) in enumerate(zip(vals_a, highs_a, hws_a)):
        x_txt = val + high_err + 1.1
        ax1.text(
            x_txt,
            i,
            f"{_fmt2(val)} ± {_fmt2(hw)}%",
            va="center",
            ha="left",
            fontsize=9.2,
            fontweight="bold",
            color="#111111",
        )

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels_a, fontsize=9.5)
    ax1.invert_yaxis()
    ax1.set_xlim(0, 81)
    ax1.set_xlabel("Percentage of TED Domains (%)", fontsize=10, fontweight="bold")
    ax1.set_title("A. Naming Specificity Tiers (N = 100,000 Domains)", fontsize=11, fontweight="bold", pad=10)
    ax1.grid(axis="x", linestyle="--", alpha=0.45, color="#bbbbbb", zorder=0)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # --- Panel B: Blind Hold-Out CATH Accuracy by Level ---
    b_res = results["panel_b_cath_validation"]
    levels = ["C", "C.A", "C.A.T", "C.A.T.H"]
    x_labels = ["Class\n(C)", "Architecture\n(C.A)", "Topology\n(C.A.T)", "Superfamily\n(C.A.T.H)"]
    series_cfg = [
        ("confirmed_direct_precision", "Confirmed Direct Precision (N = 72,130)", "#1b9e77", -1),
        ("profile_to_cath_5fold_cv", "Profile-to-CATH 5-Fold CV", "#2b6cb0", 0),
        ("full_fallback_cascade", "Full Cascade (Mask TED CATH)", "#e53e3e", 1),
    ]

    x = np.arange(len(levels), dtype=float)
    width = 0.275

    ax2.axhline(
        90.0,
        color="#222222",
        linestyle=":",
        linewidth=1.1,
        label="90% Confirmed Calibration Threshold",
        zorder=2,
    )

    for key, legend_label, color, offset_idx in series_cfg:
        stats_map = b_res[key]["protein_clustered"]
        vals = np.array([stats_map[lvl]["point_pct"] for lvl in levels])
        lows = np.array([vals[j] - stats_map[lvl]["ci95_low"] for j, lvl in enumerate(levels)])
        highs = np.array([stats_map[lvl]["ci95_high"] - vals[j] for j, lvl in enumerate(levels)])
        hws = np.array([stats_map[lvl]["halfwidth_sym"] for lvl in levels])

        x_bars = x + offset_idx * width
        ax2.bar(
            x_bars,
            vals,
            width=width,
            color=color,
            edgecolor="#222222",
            linewidth=0.8,
            label=legend_label,
            zorder=3,
        )
        ax2.errorbar(
            x_bars,
            vals,
            yerr=[lows, highs],
            fmt="none",
            ecolor="#111111",
            elinewidth=1.25,
            capsize=2.6,
            capthick=1.0,
            zorder=4,
        )

        for xb, val, high_err, hw in zip(x_bars, vals, highs, hws):
            ax2.text(
                xb,
                val + high_err + 0.28,
                f"{_fmt2(val)}\n±{_fmt2(hw)}",
                ha="center",
                va="bottom",
                fontsize=6.7,
                fontweight="bold",
                linespacing=0.95,
                color="#111111",
            )

    ax2.set_xticks(x)
    ax2.set_xticklabels(x_labels, fontsize=9.5, fontweight="bold")
    ax2.set_ylim(80.0, 104.8)
    ax2.set_yticks([80, 85, 90, 95, 100])
    ax2.set_ylabel("Empirical Accuracy / Precision (%)", fontsize=10, fontweight="bold")
    ax2.set_title("B. Blind Hold-Out CATH Accuracy by Level", fontsize=11, fontweight="bold", pad=10)
    ax2.grid(axis="y", linestyle="--", alpha=0.45, color="#bbbbbb", zorder=0)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(
        loc="lower left",
        fontsize=8.2,
        frameon=True,
        facecolor="#ffffff",
        edgecolor="#cccccc",
        framealpha=0.95,
    )

    n_clust_a = results["n_uniprot_clusters_all"]
    n_clust_b = results["n_uniprot_clusters_official_cath"]
    fig.text(
        0.5,
        0.012,
        f"Error bars denote 95% protein-clustered nonparametric bootstrap confidence intervals (B = {N_BOOT:,} resamples by UniProt ID; "
        f"{n_clust_a:,} proteins in A, {n_clust_b:,} in B).",
        ha="center",
        va="bottom",
        fontsize=8.2,
        color="#444444",
    )

    fig.subplots_adjust(left=0.175, right=0.985, top=0.91, bottom=0.17)
    FIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PATH, dpi=300)
    plt.close(fig)
    print(f"Saved figure to {FIG_PATH}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute protein-clustered bootstrap 95% CIs and render paper/figure1_validation.png."
    )
    parser.add_argument(
        "--from-summary",
        action="store_true",
        help="Render the validation figure directly from existing bootstrap_95ci results in benchmark_100000_summary.json.",
    )
    args = parser.parse_args()

    if args.from_summary:
        summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
        results = summary.get("bootstrap_95ci")
        if not results:
            raise ValueError(f"Missing 'bootstrap_95ci' in {SUMMARY_PATH}; run without --from-summary first.")
        print(f"Loaded cached bootstrap_95ci from {SUMMARY_PATH}")
    else:
        results = load_and_bootstrap()
        summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
        summary["bootstrap_95ci"] = results
        SUMMARY_PATH.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"Updated {SUMMARY_PATH} with bootstrap_95ci")

    # Print formatted tables for inspection
    print("\n=== Panel A: Protein-Clustered vs Domain-IID Bootstrap 95% CI ===")
    a_pc = results["panel_a_specificity_tiers"]["protein_clustered"]
    a_dom = results["panel_a_specificity_tiers"]["domain_iid_comparison"]
    for k, s in a_pc.items():
        d = a_dom[k]
        print(
            f"  {k:24s}: {s['point_pct']:6.3f}% | clust CI=[{s['ci95_low']:6.3f}, {s['ci95_high']:6.3f}] (±{s['halfwidth_sym']:.3f}%) "
            f"| dom CI=[{d['ci95_low']:6.3f}, {d['ci95_high']:6.3f}] (±{d['halfwidth_sym']:.3f}%)"
        )

    print("\n=== Panel B: Protein-Clustered vs Domain-IID Bootstrap 95% CI ===")
    for group in ("confirmed_direct_precision", "profile_to_cath_5fold_cv", "full_fallback_cascade"):
        print(f"  [{group}]")
        pc_map = results["panel_b_cath_validation"][group]["protein_clustered"]
        dom_map = results["panel_b_cath_validation"][group]["domain_iid_comparison"]
        for lvl in ("C", "C.A", "C.A.T", "C.A.T.H"):
            s = pc_map[lvl]
            d = dom_map[lvl]
            print(
                f"    {lvl:7s}: {s['point_pct']:6.2f}% | clust CI=[{s['ci95_low']:6.2f}, {s['ci95_high']:6.2f}] (±{s['halfwidth_sym']:.2f}%) "
                f"| dom CI=[{d['ci95_low']:6.2f}, {d['ci95_high']:6.2f}] (±{d['halfwidth_sym']:.2f}%)"
            )
    cov = results["panel_b_cath_validation"]["profile_to_cath_5fold_cv"]["coverage_pct"]
    print(f"  Profile-to-CATH coverage: {cov['point_pct']:.2f}% [{cov['ci95_low']:.2f}, {cov['ci95_high']:.2f}] (±{cov['halfwidth_sym']:.2f}%)")

    render_figure(results)


if __name__ == "__main__":
    main()
