# Accuracy, Literature Validity & Cross-Validation Benchmark (`ted-namer` — `v0.2.0`, 100,000 TED Domains)

This document reports the ground-truth classification accuracy, naming hygiene audit, and decile-by-decile statistical stability of `ted-namer` (`v0.2.0`) evaluated across **100,000 unique TED structural domains** (100 batches of 1,000 domains / 10 deciles of 10,000 domains across all kingdoms of life), reproducible via [`scripts/benchmark_100000_random.py`](../scripts/benchmark_100000_random.py) and [`scripts/reevaluate_100k.py`](../scripts/reevaluate_100k.py).

---

## 1. Blind Hold-Out CATH Classification Accuracy ($N = 72,130$)

To measure true predictive accuracy against experimental/Foldseek structural assignments, we took all **72,130 domains** in our 100,000-domain benchmark that carry an official 4-level TED CATH assignment (`ted_official`), masked out their official CATH label (`cath_label = "-"`), and evaluated three progressively stricter hold-out scenarios:

| Validation Scenario ($N = 72,130$) | Class ($C$) Accuracy | Architecture ($C.A$) Accuracy | Topology ($C.A.T$) Accuracy | Superfamily ($C.A.T.H$) Accuracy |
| :--- | :---: | :---: | :---: | :---: |
| **Scenario A: Full Fallback Cascade**<br>*(Mask TED `cath_label`; allow `G3DSA` $\to$ 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **98.2%** (`98.15%`) | **95.2%** (`95.22%`) | **93.9%** (`93.94%`) | **91.1%** (`91.10%`) |
| **Scenario B: No Structural CATH DBs**<br>*(Mask both TED `cath_label` AND `G3DSA`; allow 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **97.0%** (`96.99%`) | **90.4%** (`90.35%`) | **87.7%** (`87.73%`) | **81.3%** (`81.33%`) |
| **Scenario B1: `Profile-to-CATH` Bridge Alone**<br>*(5-fold out-of-sample Pfam/SUPERFAMILY/CDD/IPR $\to$ CATH, covering 87.9% of domains)* | **99.3%** (`99.28%`) | **98.6%** (`98.63%`) | **97.9%** (`97.93%`) | **92.3%** (`92.29%`) |
| **Scenario C: Pure 3D/SSE Rules Only**<br>*(Mask all sequence/profile annotations; geometry & SSE counts only)* | **85.6%**<br>*(96.7% pure $\alpha/\beta$)* | **41.5%** | **25.4%** | — |

### Precision of Confirmed Direct (Non-Wildcard `?`-Free) CATH Levels ($N = 72,130$)

Aligned with TED (*Lau et al., Science 2024*) and CATH v4.4 (*Waman et al., NAR 2025*), `ted-namer` avoids conflating computational estimates with CATH's `T`-level (Topology) code. A hierarchy level is emitted as a confirmed integer **only when empirical out-of-sample precision exceeds 90%**; deeper unconfirmed levels use explicit **`?` wildcards** (e.g., `2.60.40.?` or `3.?.?.?`), while numeric estimates are preserved separately in `cath_fold_hint` alongside `novel_fold_candidate = True`:

| CATH Level Emitted as Confirmed Integer (Non-`?`) | Claimed Count (out of 72,130) | Empirical Precision | 95% Protein-Clustered Bootstrap CI |
| :--- | :---: | :---: | :---: |
| **Level 1: Class ($C$)** | 72,130 (100.0%) | **98.15%** | [98.05%, 98.25%] ($\pm 0.11\%$) |
| **Level 2: Architecture ($C.A$)** | 67,115 (93.0%) | **98.57%** | [98.48%, 98.66%] ($\pm 0.09\%$) |
| **Level 3: Topology ($C.A.T$)** | 67,110 (93.0%) | **98.07%** | [97.96%, 98.18%] ($\pm 0.11\%$) |
| **Level 4: Homologous Superfamily ($C.A.T.H$)** | 48,913 (67.8%) | **99.78%** | [99.73%, 99.82%] ($\pm 0.05\%$) |

> **Key Takeaway**: Across all 100,000 domains, `72,130` (`72.13%`) have full 4-level official CATH codes, `17,130` (`17.13%`) are mapped via partial TED `T`-level (`9,692`), `Profile-to-CATH` (`4,320`), or `G3DSA` (`3,118`), and `10,740` (`10.74%`) lack homology and are explicitly flagged as `novel_fold_candidate = True` with `?` wildcard notation (`C.?.?.?`).

### Bootstrap Uncertainty Quantification ($B = 2,000$, Protein-Clustered)

Because multi-domain UniProt proteins contribute multiple TED domains (up to 6 domains per protein across `41,138` proteins in the 100,000-domain benchmark and `34,443` proteins among the `72,130` official CATH domains), treating domains as independent i.i.d. draws underestimates sampling variance. We therefore compute **nonparametric protein-clustered bootstrap 95% confidence intervals** ($B = 2,000$ resamples of `uniprot_acc` clusters with replacement, `seed = 42`, alongside domain-level i.i.d. bootstrap for comparison) via [`scripts/generate_validation_figure.py`](../scripts/generate_validation_figure.py), persisted in [`examples/benchmark_100000_summary.json`](../examples/benchmark_100000_summary.json) (`bootstrap_95ci`):

| Metric (Panel A — $N = 100,000$ Domains) | Point Estimate | Protein-Clustered 95% CI (41,138 Proteins) | Domain-IID 95% CI Comparison |
| :--- | :---: | :---: | :---: |
| **Tier 1: Specific Domain (& Subdomain)** | 62.43% (`62.425%`) | [62.04%, 62.80%] ($\pm 0.38\%$) | [62.13%, 62.72%] ($\pm 0.30\%$) |
| **Tier 2: Homologous Superfamily (Curated)** | 24.54% (`24.537%`) | [24.22%, 24.85%] ($\pm 0.32\%$) | [24.28%, 24.80%] ($\pm 0.26\%$) |
| **Tier 1 + Tier 2 Combined Resolution** | 86.96% (`86.962%`) | [86.70%, 87.23%] ($\pm 0.27\%$) | [86.75%, 87.17%] ($\pm 0.22\%$) |
| **Tier 3: Protein Family + Context** | 12.86% (`12.858%`) | [12.59%, 13.12%] ($\pm 0.27\%$) | [12.65%, 13.07%] ($\pm 0.21\%$) |
| **Tier 4/5: Fold Fallback & Uncharacterized** | 0.18% (`0.180%`) | [0.15%, 0.22%] ($\pm 0.04\%$) | [0.15%, 0.21%] ($\pm 0.03\%$) |

| Metric (Panel B — $N = 72,130$ Official CATH Domains) | Level | Point Estimate | Protein-Clustered 95% CI (34,443 Proteins) | Domain-IID 95% CI |
| :--- | :---: | :---: | :---: | :---: |
| **Profile-to-CATH 5-Fold CV** (`87.94%` coverage) | $C$ / $C.A$ / $C.A.T$ / $C.A.T.H$ | 99.28% / 98.63% / 97.93% / 92.29% | [99.21, 99.34] / [98.54, 98.73] / [97.82, 98.06] / [92.08, 92.51]% | $\pm 0.07$ / $\pm 0.09$ / $\pm 0.11$ / $\pm 0.20\%$ |
| **Full Fallback Cascade** (Mask TED CATH) | $C$ / $C.A$ / $C.A.T$ / $C.A.T.H$ | 98.15% / 95.22% / 93.94% / 91.10% | [98.05, 98.25] / [95.06, 95.39] / [93.75, 94.12] / [90.86, 91.32]% | $\pm 0.10$ / $\pm 0.16$ / $\pm 0.17$ / $\pm 0.20\%$ |

---

## 2. Literature Validity & Naming Hygiene Audit ($N = 100,000$)

Summary of cumulative naming hygiene improvements across the 100,000-domain benchmark up through **`v0.2.0`**:

| Validity / Hygiene Criterion (across 100,000 Domains) | Initial Pipeline | `v0.2.0` Final | Resolution |
| :--- | :---: | :---: | :--- |
| **`T`-prefix collisions with CATH `T`-level** | `26,187` (26.2%) | **`0` (0.0%)** | Replaced 100% with hierarchical `?` wildcards (`26,187`) + `cath_fold_hint` + `novel_fold_candidate` (`10,740`). |
| **Bare positional names (`N-terminal domain`, `NTD domain`)** | `149` | **`0` (0.0%)** | **100% eliminated** — enriched with parent protein/family context (e.g. `Argonaute N-terminal domain`). |
| **Bracketed gene / EC tags (`[bioU]`, etc.)** | `703` | **`0` (0.0%)** | **100% eliminated** by tag stripping. |
| **Adjacent double-domain suffixes (`domain-like domain`)** | `8` | **`0` (0.0%)** | **100% eliminated** by suffix deduplication. |
| **Top legacy PDB titles (`Phosphorylase Kinase domain 1`, `Transferase(Phosphotransferase) domain 1`)** | `2,146` (`1,152 + 994`) | **`0` (0.0%)** | **100% eliminated** via curated unintegrated H-names & scoping legacy PDB penalty cleanly. |
| **False `0.65x` penalty on curated Tier 1 Pfam/InterPro domains** | `820` | **`0` (0.0%)** | **100% eliminated** (e.g. `RNA polymerase Rpb1 domain 1` retains unpenalized confidence). |
| **Plain lowercase-starting names (`protein kinase ...`)** | `5,054` | **`744`** *(biochemical only)* | Plain words capitalized; `p53`, `tRNA`, `uS5`, `von Willebrand`, `beta-` preserved. |
| **Strict legacy PDB noise (`; domain \d+`, `. Chain`, `Vaccinia`)** | — | **`36` (0.036%)** | **99.96% clean** across 100,000 domains. |
| **Sentence-like descriptions / verbose >90-char names** | `2,454` / `2,230` | **`81` / `64`** | **96.7% / 97.1% reduction** via concise InterPro/short titles + `0.40x` penalty. |

---

## 3. Statistical Stability Across 10 Deciles of 10,000 Domains ($N = 100,000$)

Evaluated across **10 consecutive deciles of 10,000 unique TED domains each** (`v0.2.0`, stored in [`examples/random_100000_results.tsv.gz`](../examples/random_100000_results.tsv.gz)):

| Decile (10,000 Domains Each) | Mean Conf | Tier 1: Domain | Tier 2: Superfamily | **Tier 1 + Tier 2** | Tier 3: Family Context | Tier 4 / 5: Fold / Unchar | CATH Official | Multi-Lobe Subdomains |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **D1 (`1–10k`)** | 0.804 | 62.54% | 23.24% | **85.78%** | 14.20% | 0.02% | 69.79% | 415 |
| **D2 (`10k–20k`)** | 0.818 | 63.29% | 24.98% | **88.27%** | 10.58% | 1.15% | 73.20% | 471 |
| **D3 (`20k–30k`)** | 0.800 | 63.31% | 23.87% | **87.18%** | 12.35% | 0.47% | 70.52% | 629 |
| **D4 (`30k–40k`)** | 0.796 | 62.92% | 24.35% | **87.27%** | 12.70% | 0.03% | 72.82% | 643 |
| **D5 (`40k–50k`)** | 0.796 | 62.51% | 24.55% | **87.06%** | 12.92% | 0.02% | 72.05% | 610 |
| **D6 (`50k–60k`)** | 0.793 | 62.22% | 24.91% | **87.13%** | 12.83% | 0.04% | 72.09% | 670 |
| **D7 (`60k–70k`)** | 0.793 | 61.78% | 24.93% | **86.71%** | 13.28% | 0.01% | 73.19% | 609 |
| **D8 (`70k–80k`)** | 0.793 | 61.44% | 25.35% | **86.79%** | 13.18% | 0.03% | 72.50% | 661 |
| **D9 (`80k–90k`)** | 0.797 | 62.43% | 24.55% | **86.98%** | 12.99% | 0.03% | 73.21% | 592 |
| **D10 (`90k–100k`)** | 0.791 | 61.81% | 24.64% | **86.45%** | 13.55% | 0.00% | 71.93% | 634 |
| **TOTAL (`100,000`)** | **0.798** | **62.43%** | **24.54%** | **86.96%** | **12.86%** | **0.18%** | **72.13%** | **5,934** |

---

## 4. Leave-One-Database-Out Cross-Validation & Structural Quality Flags ($N = 100,000$)

### A. Leave-One-Database-Out Cross-Validation ($N = 36,843$)
For all **36,843 domains** across the 100,000-domain dataset where `Pfam` contributed to a Tier 1 (`specific_domain`) assignment alongside at least one other domain database (`CDD`, `SMART`, `InterPro`, etc.), masking out `Pfam` completely yielded:
- **Exact InterPro ID / Normalized Name Concordance**: **91.23%** (`33,611 / 36,843`)
- **Semantic / Functional Keyword Concordance**: **96.45%** (`35,534 / 36,843`)

### B. Structural Quality Flags (Decoupled from Naming Confidence)
| `structure_quality_flags` | Domain Count across 100,000 | Share |
| :--- | :---: | :---: |
| **`ok`** (pLDDT $\ge 70$, globular `norm_rg` $\le 0.50$) | `96,190` | **96.19%** |
| **`low_plddt`** (`plddt < 70` or member `< 60`) | `2,078` | **2.08%** |
| **`non_globular`** (`norm_rg > 0.50`) | `1,850` | **1.85%** |

### C. Calibrated Confidence Score Distribution ($N = 100,000$)
| Specificity Tier | Share of 100,000 Domains | Mean Calibrated Confidence |
| :--- | :---: | :---: |
| **Tier 1: `specific_domain`** (incl. 5,934 multi-lobe subdomains) | **62.43%** (62,425) | **0.889** |
| **Tier 2: `superfamily`** (Curated Integrated IPR + H-overrides) | **24.54%** (24,537) | **0.769** |
| **Tier 3: `protein_family_context`** (Family + positional context) | **12.86%** (12,858) | **0.422** |
| **Tier 4: `topology_fold`** (Confirmed-`T` CATH Topology fallback) | **0.07%** (70) | **0.209** |
| **Tier 5: `novel_uncharacterized`** | **0.11%** (110) | **0.151** |
| **Overall (100,000 TED Domains)** | **100.00%** | **0.798** |
