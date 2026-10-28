# Accuracy & Cross-Validation Benchmark (`ted-namer` — 100,000 TED Domains)

This document reports the ground-truth accuracy and decile-by-decile statistical stability of `ted-namer` evaluated across **100,000 unique TED structural domains** (100 batches of 1,000 domains / 10 deciles of 10,000 domains across all kingdoms of life), reproducible via [`scripts/benchmark_100000_random.py`](../scripts/benchmark_100000_random.py).

---

## 1. Blind Hold-Out CATH Classification Accuracy ($N = 72,130$)

To measure true predictive accuracy rather than internal confidence heuristics, we took all **72,130 domains** in our 100,000-domain benchmark that possess an official 4-level TED CATH assignment (`ted_official`), masked out their official CATH label (`cath_label = "-"`), and evaluated three progressively stricter hold-out scenarios:

| Validation Scenario ($N = 72,130$) | Class ($C$) Accuracy | Architecture ($C.A$) Accuracy | Topology ($C.A.T$) Accuracy | Superfamily ($C.A.T.H$) Accuracy |
| :--- | :---: | :---: | :---: | :---: |
| **Scenario A: Full Fallback Cascade**<br>*(Mask TED `cath_label`; allow `G3DSA` $\to$ 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **98.2%** (`98.16%`) | **95.2%** (`95.22%`) | **93.9%** (`93.94%`) | **91.1%** (`91.10%`) |
| **Scenario B: No Structural CATH DBs**<br>*(Mask both TED `cath_label` AND `G3DSA`; allow 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **97.0%** (`96.99%`) | **90.4%** (`90.35%`) | **87.7%** (`87.73%`) | **81.3%** (`81.33%`) |
| **Scenario B1: `Profile-to-CATH` Bridge Alone**<br>*(5-fold out-of-sample Pfam/SUPERFAMILY/CDD/IPR $\to$ CATH, covering 87.9% of domains)* | **99.3%** (`99.28%`) | **98.6%** (`98.63%`) | **97.9%** (`97.94%`) | **92.3%** (`92.29%`) |
| **Scenario C: Pure 3D/SSE Rules Only**<br>*(Mask all sequence/profile annotations; geometry & SSE counts only)* | **85.6%**<br>*(96.7% pure $\alpha/\beta$)* | **41.5%** | **25.4%** | — |

### Precision of Non-`T`-Prefixed (Direct) CATH Levels ($N = 72,130$)

In `ted-namer`, a CATH level is written as a plain number only when its empirical precision exceeds **90%**; otherwise it carries the `T` estimation prefix (e.g., `3.T40.T50.T10`). Across the 72,130 held-out official CATH domains in Scenario A:

| CATH Level Emitted Without `T` Prefix | Claimed Count (out of 72,130) | Empirical Precision |
| :--- | :---: | :---: |
| **Level 1: Class ($C$)** | 72,130 (100.0%) | **98.16%** |
| **Level 2: Architecture ($C.A$)** | 67,115 (93.0%) | **98.57%** |
| **Level 3: Topology ($C.A.T$)** | 67,110 (93.0%) | **98.07%** |
| **Level 4: Homologous Superfamily ($C.A.T.H$)** | 48,913 (67.8%) | **99.78%** |

> **Key Takeaway**: Because pure 3D/SSE geometric rules without sequence homology achieve **85.6%** accuracy on Class $C$ (and **96.7%** on pure $\alpha$ or pure $\beta$ domains) but only **41.5%** on Architecture ($C.A$) and **25.4%** on Topology ($C.A.T$), pure 3D/SSE predictions strictly emit only Class $C$ as a plain number and prefix estimated Architecture, Topology, and Superfamily levels with `T` (e.g., `3.T40.T50.T10`).

---

## 2. Statistical Stability Across 10 Deciles of 10,000 Domains ($N = 100,000$)

To verify statistical convergence and absence of drift or overfitting, the 100,000 TED domains were evaluated across **10 consecutive deciles of 10,000 unique domains each**:

| Decile (10,000 Domains Each) | Mean Conf | Tier 1: Domain | Tier 2: Superfamily | **Tier 1 + Tier 2** | Tier 3: Family Context | Tier 4: Fold | CATH Official | Multi-Lobe Subdomains |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **D1 (`1–10k`)** | 0.828 | 64.41% | 23.87% | **88.28%** | 11.69% | 0.03% | 69.79% | 430 |
| **D2 (`10k–20k`)** | 0.843 | 64.71% | 26.08% | **90.79%** | 8.06% | 1.15% | 73.20% | 457 |
| **D3 (`20k–30k`)** | 0.825 | 62.08% | 26.82% | **88.90%** | 10.63% | 0.47% | 70.52% | 632 |
| **D4 (`30k–40k`)** | 0.821 | 62.92% | 26.91% | **89.83%** | 10.14% | 0.03% | 72.82% | 655 |
| **D5 (`40k–50k`)** | 0.820 | 62.31% | 27.29% | **89.60%** | 10.38% | 0.02% | 72.05% | 629 |
| **D6 (`50k–60k`)** | 0.819 | 62.26% | 27.99% | **90.25%** | 9.70% | 0.05% | 72.09% | 686 |
| **D7 (`60k–70k`)** | 0.820 | 61.48% | 28.15% | **89.63%** | 10.36% | 0.01% | 73.19% | 631 |
| **D8 (`70k–80k`)** | 0.819 | 61.49% | 28.35% | **89.84%** | 10.13% | 0.03% | 72.50% | 663 |
| **D9 (`80k–90k`)** | 0.822 | 62.00% | 28.10% | **90.10%** | 9.87% | 0.03% | 73.21% | 603 |
| **D10 (`90k–100k`)** | 0.818 | 61.94% | 27.42% | **89.36%** | 10.64% | 0.00% | 71.93% | 636 |
| **TOTAL (`100,000`)** | **0.823** | **62.56%** | **27.10%** | **89.66%** | **10.16%** | **0.18%** | **72.13%** | **6,022** |

---

## 3. Domain Naming Accuracy & Cross-Database Validation ($N = 100,000$)

### A. Leave-One-Database-Out Cross-Validation ($N = 34,648$)
For all **34,648 domains** across the 100,000-domain dataset where `Pfam` contributed to a Tier 1 (`specific_domain`) assignment alongside at least one other domain database (`CDD`, `SMART`, `PROSITE`, `InterPro`), we masked out `Pfam` completely and tested whether the remaining databases independently recovered the same domain designation:

- **Exact InterPro ID / Normalized Name Concordance**: **90.05%** (`31,201 / 34,648`)
- **Semantic / Functional Keyword Concordance**: **95.69%** (`33,156 / 34,648`)

### B. Elimination of Legacy PDB-Header Naming Artifacts
By always prioritizing the curated integrated InterPro superfamily name (`integrated_ipr`) over raw `G3DSA` headers and stripping PDB chain/domain suffixes:
- Legacy PDB-header artifacts across all 100,000 domains were restricted to **2 / 100,000 (`0.002%`)** (**>99.9% reduction**).

### C. Resolution of Multi-Lobe Functional Domains
- With subdomain lobe detection ($\text{Cov}_{\text{TED}} \ge 0.80$, $0.20 \le \text{Cov}_{\text{Annot}} < 0.50$), **6,022 multi-lobe domains** (`6.02%` of all 100,000 TED domains) are accurately named with their positional subdomain identity (`N-terminal`, `central`, or `C-terminal subdomain`).

---

## 4. Calibrated Confidence Score Distribution ($N = 100,000$)

| Specificity Tier | Share of 100,000 Domains | Mean Calibrated Confidence |
| :--- | :---: | :---: |
| **Tier 1: `specific_domain`** | **62.56%** (62,560) | **0.899** |
| **Tier 2: `superfamily`** | **27.10%** (27,098) | **0.809** |
| **Tier 3: `protein_family_context`** | **10.16%** (10,160) | **0.410** |
| **Tier 4: `topology_fold`** | **0.18%** (182) | **0.235** |
| **Overall (100,000 TED Domains)** | **100.00%** | **0.823** |
