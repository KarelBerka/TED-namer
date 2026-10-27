# Accuracy & Cross-Validation Benchmark (`ted-namer`)

This document reports the ground-truth accuracy of `ted-namer` evaluated via **blind hold-out testing** and **5-fold cross-validation** across **10,000 diverse TED structural domains** (`examples/random_10000_results.tsv`), reproducible via [`scripts/evaluate_accuracy.py`](../scripts/evaluate_accuracy.py).

---

## 1. Blind Hold-Out CATH Classification Accuracy ($N = 6,979$)

To measure true predictive accuracy rather than internal confidence heuristics, we took all **6,979 domains** in our 10,000-domain benchmark that possess an official 4-level TED CATH assignment (`ted_official`), masked out their official CATH label (`cath_label = "-"`), and evaluated three progressively stricter hold-out scenarios:

| Validation Scenario | Class ($C$) Accuracy | Architecture ($C.A$) Accuracy | Topology ($C.A.T$) Accuracy | Superfamily ($C.A.T.H$) Accuracy |
| :--- | :---: | :---: | :---: | :---: |
| **Scenario A: Full Fallback Cascade**<br>*(Mask TED `cath_label`; allow `G3DSA` $\to$ 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **97.9%** | **93.9%** | **92.3%** | **89.4%** |
| **Scenario B: No Structural CATH DBs**<br>*(Mask both TED `cath_label` AND `G3DSA`; allow 5-fold CV `Profile-to-CATH` $\to$ `3D/SSE rules`)* | **96.0%** | **85.6%** | **81.7%** | **74.8%** |
| **Scenario B1: `Profile-to-CATH` Bridge Alone**<br>*(5-fold out-of-sample Pfam/SUPERFAMILY/CDD/IPR $\to$ CATH, covering 81.6% of domains)* | **99.0%** | **97.3%** | **96.4%** | **91.4%** |
| **Scenario C: Pure 3D/SSE Rules Only**<br>*(Mask all sequence/profile annotations; geometry & SSE counts only)* | **85.2%**<br>*(95.6% pure $\alpha/\beta$)* | **40.6%** | **24.4%** | — |

### Precision of Non-`T`-Prefixed (Direct) CATH Levels

In `ted-namer`, a CATH level is written as a plain number only when its empirical precision exceeds **90%**; otherwise it carries the `T` estimation prefix (e.g., `3.T40.T50.T10`). Across the 6,979 held-out domains in Scenario A:

| CATH Level Emitted Without `T` Prefix | Claimed Count (out of 6,979) | Empirical Precision |
| :--- | :---: | :---: |
| **Level 1: Class ($C$)** | 6,979 (100.0%) | **97.9%** |
| **Level 2: Architecture ($C.A$)** | 6,363 (91.2%) | **98.2%** |
| **Level 3: Topology ($C.A.T$)** | 6,362 (91.2%) | **97.9%** |
| **Level 4: Homologous Superfamily ($C.A.T.H$)** | 4,875 (69.9%) | **99.8%** |

> **Key Takeaway**: Because pure 3D/SSE geometric rules without sequence homology achieve **85.2%** accuracy on Class $C$ (and **95.6%** on pure $\alpha$ or pure $\beta$ domains) but only **40.6%** on Architecture ($C.A$) and **24.4%** on Topology ($C.A.T$), pure 3D/SSE predictions strictly emit only Class $C$ as a plain number and prefix estimated Architecture, Topology, and Superfamily levels with `T` (e.g., `3.T40.T50.T10`).

---

## 2. Domain Naming Accuracy & Cross-Database Validation ($N = 10,000$)

### A. Leave-One-Database-Out Cross-Validation ($N = 3,397$)
For all 3,397 domains where `Pfam` contributed to a Tier 1 (`specific_domain`) assignment alongside at least one other domain database (`CDD`, `SMART`, `PROSITE`, `InterPro`), we masked out `Pfam` completely and tested whether the remaining databases independently recovered the same domain designation:

- **Exact InterPro ID / Normalized Name Concordance**: **87.6%** (`2,975 / 3,397`)
- **Semantic / Functional Keyword Concordance**: **95.3%** (`3,239 / 3,397`)

### B. Elimination of Legacy PDB-Header Naming Artifacts
In raw CATH/Gene3D (`G3DSA`), several superfamilies carry historical 1990s PDB structure headers as their raw names (for example, `G3DSA:3.40.50.150` is named `"Vaccinia Virus protein VP39"` instead of its curated InterPro name `"S-adenosyl-L-methionine-dependent methyltransferase superfamily"`, and `G3DSA:3.30.200.20` is named `"Phosphorylase Kinase; domain 1"`).

By always prioritizing the curated integrated InterPro superfamily name (`integrated_ipr`) over raw `G3DSA` headers and stripping PDB chain/domain suffixes:
- Legacy PDB-header artifacts across the 10,000-domain benchmark dropped from **208 domains to 2** (**99.0% reduction**).

### C. Resolution of Multi-Lobe Functional Domains
TED frequently chops a single canonical sequence domain (such as a 260-residue `Protein kinase domain` `PF00069` or `Cytochrome P450` `PF00067`) into 2–3 structural lobes/subdomains. Because each lobe covers only ~25–35% of the full Pfam profile ($\text{IoU} < 0.35$), naive IoU thresholding previously demoted the specific domain to Tier 3 (`protein_family_context`), allowing a generic Tier 2 fold name to win.
- With subdomain lobe detection ($\text{Cov}_{\text{TED}} \ge 0.80$, $0.20 \le \text{Cov}_{\text{Annot}} < 0.50$), **430 multi-lobe domains** in the 10,000-domain dataset are now accurately named with their positional subdomain identity (e.g., `Protein kinase N-terminal subdomain`, `Cytochrome P450 central subdomain`).

---

## 3. Calibrated Confidence Score Distribution ($N = 10,000$)

After removing the artificial single-domain frequency floor (`freq = 1/1`) and replacing hardcoded whole-protein fallback overlaps with the true domain-to-protein length ratio ($L_{\text{domain}} / L_{\text{protein}}$), the confidence score accurately reflects annotation specificity and boundary agreement:

| Specificity Tier | Share of 10,000 Domains | Mean Calibrated Confidence |
| :--- | :---: | :---: |
| **Tier 1: `specific_domain`** | **64.4%** | **0.905** |
| **Tier 2: `superfamily`** | **23.9%** | **0.817** |
| **Tier 3: `protein_family_context`** | **11.7%** | **0.426** |
| **Tier 4: `topology_fold`** | **0.03%** (3 / 10,000) | **0.235** |
| **Overall (10,000 TED Domains)** | **100.0%** | **0.827** |
