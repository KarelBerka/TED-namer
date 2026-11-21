# Accuracy, Literature Validity & Cross-Validation Benchmark (`ted-namer` — 100,000 TED Domains)

This document reports the ground-truth classification accuracy, naming hygiene audit, and decile-by-decile statistical stability of `ted-namer` evaluated across **100,000 unique TED structural domains** (100 batches of 1,000 domains / 10 deciles of 10,000 domains across all kingdoms of life), reproducible via [`scripts/benchmark_100000_random.py`](../scripts/benchmark_100000_random.py) and [`scripts/reevaluate_100k.py`](../scripts/reevaluate_100k.py).

---

## 1. Blind Hold-Out CATH Classification Accuracy ($N = 72,130$)

To measure true predictive accuracy against experimental/Foldseek structural assignments, we took all **72,130 domains** in our 100,000-domain benchmark that carry an official 4-level TED CATH assignment (`ted_official`), masked out their official CATH label (`cath_label = "-"`), and evaluated three progressively stricter hold-out scenarios:

| Validation Scenario ($N = 72,130$) | Class ($C$) Accuracy | Architecture ($C.A$) Accuracy | Topology ($C.A.T$) Accuracy | Superfamily ($C.A.T.H$) Accuracy |
| :--- | :---: | :---: | :---: | :---: |
| **Scenario A: Full Fallback Cascade**<br>*(Mask TED `cath_label`; allow `G3DSA` $\to$ 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **98.2%** (`98.16%`) | **95.2%** (`95.22%`) | **93.9%** (`93.94%`) | **91.1%** (`91.10%`) |
| **Scenario B: No Structural CATH DBs**<br>*(Mask both TED `cath_label` AND `G3DSA`; allow 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **97.0%** (`96.99%`) | **90.4%** (`90.35%`) | **87.7%** (`87.73%`) | **81.3%** (`81.33%`) |
| **Scenario B1: `Profile-to-CATH` Bridge Alone**<br>*(5-fold out-of-sample Pfam/SUPERFAMILY/CDD/IPR $\to$ CATH, covering 87.9% of domains)* | **99.3%** (`99.28%`) | **98.6%** (`98.63%`) | **97.9%** (`97.94%`) | **92.3%** (`92.29%`) |
| **Scenario C: Pure 3D/SSE Rules Only**<br>*(Mask all sequence/profile annotations; geometry & SSE counts only)* | **85.6%**<br>*(96.7% pure $\alpha/\beta$)* | **41.5%** | **25.4%** | — |

### Precision of Confirmed Direct (Non-Wildcard `?`-Free) CATH Levels ($N = 72,130$)

Aligned with TED (*Lau et al., Science 2024*) and CATH v4.4 (*Waman et al., NAR 2025*), `ted-namer` avoids conflating computational estimates with CATH's `T`-level (Topology) code. A hierarchy level is emitted as a confirmed integer **only when empirical out-of-sample precision exceeds 90%**; deeper unconfirmed levels use explicit **`?` wildcards** (e.g., `2.60.40.?` or `3.?.?.?`), while numeric estimates are preserved separately in `cath_fold_hint` alongside `novel_fold_candidate = True`:

| CATH Level Emitted as Confirmed Integer (Non-`?`) | Claimed Count (out of 72,130) | Empirical Precision |
| :--- | :---: | :---: |
| **Level 1: Class ($C$)** | 72,130 (100.0%) | **98.16%** |
| **Level 2: Architecture ($C.A$)** | 67,115 (93.0%) | **98.57%** |
| **Level 3: Topology ($C.A.T$)** | 67,110 (93.0%) | **98.07%** |
| **Level 4: Homologous Superfamily ($C.A.T.H$)** | 48,913 (67.8%) | **99.78%** |

> **Key Takeaway**: Across all 100,000 domains, `72,130` (`72.1%`) have full 4-level official CATH codes, `17,130` (`17.1%`) are mapped via partial TED `T`-level (`9,692`), `Profile-to-CATH` (`4,320`), or `G3DSA` (`3,118`), and `10,740` (`10.74%`) lack homology and are explicitly flagged as `novel_fold_candidate = True` with `?` wildcard notation (`C.?.?.?`).

---

## 2. Before vs. After Literature & Naming Hygiene Audit ($N = 100,000$)

Comparing the refreshed 100,000-domain outputs against the pre-audit pipeline highlights the exact impact of enforcing InterPro `annotation_type` rules, curated unintegrated CATH H-overrides, naming sanitization, and strict `?` wildcard notation:

| Validity / Hygiene Criterion (across 100,000 Domains) | Before Audit | After Audit & Upgrade | Change / Resolution |
| :--- | :---: | :---: | :--- |
| **Maturing `T`-prefix collisions with CATH `T`-level** | `26,187` (26.2%) | **`0` (0.0%)** | Replaced 100% with hierarchical `?` wildcards (`26,187`) + `cath_fold_hint` + `novel_fold_candidate` (`10,740`). |
| **Bare positional names (`N-terminal domain`, `NTD domain`)** | `149` | **`0` (0.0%)** | **100% eliminated** — enriched with parent protein/family context (e.g. `argonaute N-terminal domain`). |
| **Bracketed gene / EC tags (`[bioU]`, etc.)** | `703` | **`0` (0.0%)** | **100% eliminated** by tag stripping. |
| **Adjacent double-domain suffixes (`domain-like domain`)** | `8` | **`0` (0.0%)** | **100% eliminated** by suffix deduplication. |
| **Top legacy PDB titles (`Phosphorylase Kinase domain 1`, `Transferase(Phosphotransferase) domain 1`)** | `2,146` (`1,152 + 994`) | **`0` (0.0%)** | **100% eliminated** via curated unintegrated H-names & removing inflated `cath_h` boost. |
| **Strict legacy PDB noise (`; domain \d+`, `Chain`, `Vaccinia`)** | — | **`36` (0.036%)** | **99.96% clean** across 100,000 domains. |
| **Sentence-like descriptions / verbose >90-char names** | `2,454` / `2,230` | **`82` / `65`** | **96.7% / 97.1% reduction** by preferring concise InterPro/short titles + `0.40x` penalty. |

---

## 3. Statistical Stability Across 10 Deciles of 10,000 Domains ($N = 100,000$)

Evaluated across **10 consecutive deciles of 10,000 unique TED domains each** under the strict updated rules:

| Decile (10,000 Domains Each) | Mean Conf | Tier 1: Domain | Tier 2: Superfamily | **Tier 1 + Tier 2** | Tier 3: Family Context | Tier 4 / 5: Fold / Unchar | CATH Official | Multi-Lobe Subdomains |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **D1 (`1–10k`)** | 0.800 | 62.55% | 23.11% | **85.66%** | 14.32% | 0.02% | 69.79% | 423 |
| **D2 (`10k–20k`)** | 0.813 | 63.15% | 24.95% | **88.10%** | 10.75% | 1.15% | 73.20% | 478 |
| **D3 (`20k–30k`)** | 0.794 | 63.11% | 23.90% | **87.01%** | 12.52% | 0.47% | 70.52% | 640 |
| **D4 (`30k–40k`)** | 0.792 | 62.64% | 24.50% | **87.14%** | 12.83% | 0.03% | 72.82% | 639 |
| **D5 (`40k–50k`)** | 0.790 | 62.00% | 24.98% | **86.98%** | 13.00% | 0.02% | 72.05% | 603 |
| **D6 (`50k–60k`)** | 0.787 | 62.00% | 24.99% | **86.99%** | 12.97% | 0.04% | 72.09% | 667 |
| **D7 (`60k–70k`)** | 0.787 | 61.45% | 25.16% | **86.61%** | 13.38% | 0.01% | 73.19% | 601 |
| **D8 (`70k–80k`)** | 0.788 | 61.09% | 25.54% | **86.63%** | 13.35% | 0.02% | 72.50% | 651 |
| **D9 (`80k–90k`)** | 0.792 | 62.02% | 24.90% | **86.92%** | 13.05% | 0.03% | 73.21% | 578 |
| **D10 (`90k–100k`)** | 0.786 | 61.39% | 24.89% | **86.28%** | 13.72% | 0.00% | 71.93% | 621 |
| **TOTAL (`100,000`)** | **0.793** | **62.14%** | **24.69%** | **86.83%** | **12.99%** | **0.18%** | **72.13%** | **5,901** |

---

## 4. Leave-One-Database-Out Cross-Validation & Structural Quality Flags ($N = 100,000$)

### A. Leave-One-Database-Out Cross-Validation ($N = 36,511$)
For all **36,511 domains** across the 100,000-domain dataset where `Pfam` contributed to a Tier 1 (`specific_domain`) assignment alongside at least one other domain database (`CDD`, `SMART`, `InterPro`, etc.), masking out `Pfam` completely yielded:
- **Exact InterPro ID / Normalized Name Concordance**: **91.06%** (`33,247 / 36,511`) — *up from 90.05% before hygiene fixes!*
- **Semantic / Functional Keyword Concordance**: **96.39%** (`35,192 / 36,511`)

### B. Structural Quality Flags (Decoupled from Naming Confidence)
| `structure_quality_flags` | Domain Count across 100,000 | Share |
| :--- | :---: | :---: |
| **`ok`** (pLDDT $\ge 70$, globular `norm_rg` $\le 0.50$) | `96,190` | **96.19%** |
| **`low_plddt`** (`plddt < 70` or member `< 60`) | `2,078` | **2.08%** |
| **`non_globular`** (`norm_rg > 0.50`) | `1,850` | **1.85%** |

### C. Calibrated Confidence Score Distribution ($N = 100,000$)
| Specificity Tier | Share of 100,000 Domains | Mean Calibrated Confidence |
| :--- | :---: | :---: |
| **Tier 1: `specific_domain`** (incl. 5,901 multi-lobe subdomains) | **62.14%** (62,140) | **0.884** |
| **Tier 2: `superfamily`** (Curated Integrated IPR + H-overrides) | **24.69%** (24,692) | **0.765** |
| **Tier 3: `protein_family_context`** (Family + positional context) | **12.99%** (12,989) | **0.420** |
| **Tier 4: `topology_fold`** (Confirmed-`T` CATH Topology fallback) | **0.07%** (69) | **0.203** |
| **Tier 5: `novel_uncharacterized`** | **0.11%** (110) | **0.151** |
| **Overall (100,000 TED Domains)** | **100.00%** | **0.793** |
