# Accuracy & Cross-Validation Benchmark (`ted-namer`)

This document reports the ground-truth accuracy and batch-to-batch stability of `ted-namer` evaluated across **20,000 randomly sampled TED structural domains** (20 independent batches of 1,000 domains across all kingdoms of life: Batches 1–10 and out-of-sample replication Batches 11–20), reproducible via [`scripts/benchmark_20000_random.py`](../scripts/benchmark_20000_random.py) and [`scripts/evaluate_accuracy.py`](../scripts/evaluate_accuracy.py).

---

## 1. Blind Hold-Out CATH Classification Accuracy ($N = 14,299$)

To measure true predictive accuracy rather than internal confidence heuristics, we took all **14,299 domains** in our 20,000-domain benchmark that possess an official 4-level TED CATH assignment (`ted_official`), masked out their official CATH label (`cath_label = "-"`), and evaluated three progressively stricter hold-out scenarios:

| Validation Scenario | Class ($C$) Accuracy | Architecture ($C.A$) Accuracy | Topology ($C.A.T$) Accuracy | Superfamily ($C.A.T.H$) Accuracy |
| :--- | :---: | :---: | :---: | :---: |
| **Scenario A: Full Fallback Cascade**<br>*(Mask TED `cath_label`; allow `G3DSA` $\to$ 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **98.1%** | **95.1%** | **93.8%** | **91.1%** |
| **Scenario B: No Structural CATH DBs**<br>*(Mask both TED `cath_label` AND `G3DSA`; allow 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **96.2%** | **87.9%** | **84.8%** | **78.9%** |
| **Scenario B1: `Profile-to-CATH` Bridge Alone**<br>*(5-fold out-of-sample Pfam/SUPERFAMILY/CDD/IPR $\to$ CATH, covering 85.4% of domains)* | **99.0%** | **97.7%** | **96.8%** | **92.2%** |
| **Scenario C: Pure 3D/SSE Rules Only**<br>*(Mask all sequence/profile annotations; geometry & SSE counts only)* | **85.9%**<br>*(95.8% pure $\alpha/\beta$)* | **41.2%** | **24.6%** | — |

### Precision of Non-`T`-Prefixed (Direct) CATH Levels

In `ted-namer`, a CATH level is written as a plain number only when its empirical precision exceeds **90%**; otherwise it carries the `T` estimation prefix (e.g., `3.T40.T50.T10`). Across the 14,299 held-out official CATH domains in Scenario A:

| CATH Level Emitted Without `T` Prefix | Claimed Count (out of 14,299) | Empirical Precision |
| :--- | :---: | :---: |
| **Level 1: Class ($C$)** | 14,299 (100.0%) | **98.1%** |
| **Level 2: Architecture ($C.A$)** | 13,295 (93.0%) | **98.4%** |
| **Level 3: Topology ($C.A.T$)** | 13,293 (93.0%) | **98.2%** |
| **Level 4: Homologous Superfamily ($C.A.T.H$)** | 10,274 (71.9%) | **99.8%** |

> **Key Takeaway**: Because pure 3D/SSE geometric rules without sequence homology achieve **85.9%** accuracy on Class $C$ (and **95.8%** on pure $\alpha$ or pure $\beta$ domains) but only **41.2%** on Architecture ($C.A$) and **24.6%** on Topology ($C.A.T$), pure 3D/SSE predictions strictly emit only Class $C$ as a plain number and prefix estimated Architecture, Topology, and Superfamily levels with `T` (e.g., `3.T40.T50.T10`).

---

## 2. Out-of-Sample Replication Across 20,000 TED Domains (Batches 1–10 vs. Batches 11–20)

To verify that all methodology improvements generalize without overfitting to the initial 10,000 domains (Batches 1–10), we sampled and evaluated an additional **10 independent batches of 1,000 strictly disjoint TED domains** (Batches 11–20):

| Metric | Batches 1–10 (Initial 10,000) | Batches 11–20 (Replication 10,000) | **Combined Total (20,000 Domains)** |
| :--- | :---: | :---: | :---: |
| **Domain Count** | 10,000 | 10,000 | **20,000** |
| **Mean Calibrated Confidence** | 0.827 | 0.843 | **0.835** |
| **Tier 1: `specific_domain` (incl. subdomains)** | 64.4% | 64.7% | **64.6%** (12,912) |
| **Tier 2: `superfamily` (Integrated IPR priority)** | 23.9% | 26.1% | **25.0%** (4,995) |
| **Tier 1 + Tier 2 (Domain or Superfamily)** | **88.3%** | **90.8%** | **89.5%** (17,907) |
| **Tier 3: `protein_family_context`** | 11.7% | 8.1% | **9.9%** (1,975) |
| **Tier 4: `topology_fold`** | 0.03% | 1.15% | **0.59%** (118) |
| **CATH: Official 4-Level (`ted_official`)** | 69.8% | 73.2% | **71.5%** |
| **CATH: TED 3-Level (`ted_t_plus_*`)** | 11.6% | 10.1% | **10.8%** |
| **CATH: `Profile-to-CATH` + `G3DSA` Overlap** | 8.7% | 6.7% | **7.7%** |
| **CATH: Pure 3D/SSE Rules (`rule_*`)** | 9.9% | 10.0% | **10.0%** |

---

## 3. Domain Naming Accuracy & Cross-Database Validation ($N = 20,000$)

### A. Leave-One-Database-Out Cross-Validation ($N = 6,842$)
For all 6,842 domains across the 20,000-domain dataset where `Pfam` contributed to a Tier 1 (`specific_domain`) assignment alongside at least one other domain database (`CDD`, `SMART`, `PROSITE`, `InterPro`), we masked out `Pfam` completely and tested whether the remaining databases independently recovered the same domain designation:

- **Exact InterPro ID / Normalized Name Concordance**: **84.4%** (`5,776 / 6,842`)
- **Semantic / Functional Keyword Concordance**: **94.1%** (`6,436 / 6,842`)

### B. Elimination of Legacy PDB-Header Naming Artifacts
In raw CATH/Gene3D (`G3DSA`), several superfamilies carry historical 1990s PDB structure headers as their raw names (for example, `G3DSA:3.40.50.150` is named `"Vaccinia Virus protein VP39"` instead of its curated InterPro name `"S-adenosyl-L-methionine-dependent methyltransferase superfamily"`, and `G3DSA:3.30.200.20` is named `"Phosphorylase Kinase; domain 1"`).

By always prioritizing the curated integrated InterPro superfamily name (`integrated_ipr`) over raw `G3DSA` headers and stripping PDB chain/domain suffixes:
- Legacy PDB-header artifacts across all 20,000 domains were restricted to **2 / 20,000 (0.01%)** (**>99% reduction**).

### C. Resolution of Multi-Lobe Functional Domains
TED frequently chops a single canonical sequence domain (such as a 260-residue `Protein kinase domain` `PF00069` or `Cytochrome P450` `PF00067`) into 2–3 structural lobes/subdomains.
- With subdomain lobe detection ($\text{Cov}_{\text{TED}} \ge 0.80$, $0.20 \le \text{Cov}_{\text{Annot}} < 0.50$), **887 multi-lobe domains** (`4.4%` of all 20,000 TED domains) are accurately named with their positional subdomain identity (e.g., `Protein kinase N-terminal subdomain`, `Cytochrome P450 central subdomain`).

---

## 4. Calibrated Confidence Score Distribution ($N = 20,000$)

| Specificity Tier | Share of 20,000 Domains | Mean Calibrated Confidence |
| :--- | :---: | :---: |
| **Tier 1: `specific_domain`** | **64.6%** (12,912) | **0.911** |
| **Tier 2: `superfamily`** | **25.0%** (4,995) | **0.816** |
| **Tier 3: `protein_family_context`** | **9.9%** (1,975) | **0.428** |
| **Tier 4: `topology_fold`** | **0.59%** (118) | **0.235** |
| **Overall (20,000 TED Domains)** | **100.0%** | **0.835** |
