---
title: 'ted-namer: Automated biological naming and calibrated CATH v4.4 classification for structural domains in The Encyclopedia of Domains (TED)'
tags:
  - Python
  - bioinformatics
  - structural biology
  - protein domains
  - AlphaFold
  - TED
  - CATH
  - InterPro
authors:
  - name: Karel Berka
    orcid: 0000-0001-9472-2589
    affiliation: 1
affiliations:
  - name: Department of Physical Chemistry, Faculty of Science, Palacký University Olomouc, 17. listopadu 12, 771 46 Olomouc, Czech Republic
    index: 1
    ror: 04qxnmv42
date: 21 November 2026
bibliography: paper.bib
---

# Summary

**The Encyclopedia of Domains (TED)** catalogs over 365 million structural domains segmented across the AlphaFold Protein Structure Database (AFDB) [@Jumper2021AlphaFold; @Varadi2024AFDB; @Lau2024TED] and underpins **CATH v4.4**, which integrated ~90 million predicted domains and introduced Class 6 (`Special`) alongside new architectures and superfamilies [@Waman2025CATHv44]. While TED provides continuous and discontinuous boundaries, consensus quality (`high`, `medium`, `low`), and partial CATH fold labels, individual domains and clusters carry only numeric tags (e.g., `AF-P04637-F1-model_v4_TED01`).

`ted-namer` is an open-source Python library and CLI tool (`ted-namer inspect` and `ted-namer name`) that assigns standardized human-readable domain names, specificity tiers, calibrated confidence scores, structural quality flags, and provenance-aware CATH v4.4 classifications to TED domains and multi-member structural clusters. Across a 100,000-domain benchmark, `ted-namer` resolves specific functional domain or homologous superfamily names for 87.0% of domains (99.8% including positional family context) with 100% CATH Class coverage.

# Statement of Need

Translating numeric TED identifiers into valid biological names presents five challenges that defeat naive sequence lookups:

1. **1D vs. 3D Boundary Discrepancies and Multi-Lobe Splitting:** Structural choppers such as Merizo [@Lau2023Merizo] and Chainsaw [@Wells2024Chainsaw] produce discontinuous segments (`1-146_284-319`) and frequently split canonical sequence domains (e.g., protein kinases) into separate N- and C-terminal lobes that fail naive Jaccard thresholds ($\text{IoU} \ge 0.35$).
2. **Legacy `CATH-Gene3D` PDB Headers and Score Inflation:** Unintegrated `CATH-Gene3D` entries often retain 1990s PDB headers (`"Phosphorylase Kinase; domain 1"`, `"Transferase(Phosphotransferase) domain 1"`), which outscore curated Tier 1 domain hits when `cath_h` receives fixed high overlap weights.
3. **`annotation_type` Conflation and Text Noise:** Whole-protein family databases (`HAMAP`, `NCBIfam` `family`) [@Blum2025InterPro], short signatures (`PROSITE conserved_site`), UniProt `region` prose, bracketed gene tags (`[bioU]`), and bare positional strings (`"N-terminal domain"`) require strict type- and syntax-aware filtering.
4. **Cross-Database Synonym Redundancy:** Synonymous `Pfam` [@PaysanLafosse2025Pfam], `CDD` [@Wang2023CDD], and `InterPro` hits must be unified rather than double-counted.
5. **Hypothetical Novel Folds:** In TED and CATH v4.4, `T` denotes a confirmed *Topology* match, whereas ~18--20% of TED domains lack homology (`cath_label = "-"`) and represent hypothetical novel folds whose geometric hints should not collide with confirmed `T`-level codes.

# State of the Field

Domain resources address parts of this workflow independently: `Merizo` [@Lau2023Merizo] and `Chainsaw` [@Wells2024Chainsaw] predict structural boundaries and **TED** [@Lau2024TED] assigns Foldseek/Foldclass `cath_label` codes without biological names; `InterPro` [@Blum2025InterPro], `Pfam` [@PaysanLafosse2025Pfam], and `CDD` [@Wang2023CDD] curate 1D functional signatures without reconciling 3D multi-lobe choppings or cluster consensus; and **CATH v4.4** [@Waman2025CATHv44] (`CATH-AlphaFlow` [@Waman2024CATH]) leaves unmapped AlphaFold folds unclassified. `ted-namer` unifies 3D boundary geometry, CATH v4.4 resolution, and curated InterPro/UniProt [@UniProt2025] nomenclature.

# Methodology and Software Architecture

For single domains or homologous clusters, `TedNamingPipeline` executes six cached asynchronous stages:

1. **Multi-Source Retrieval:** Queries TED (`chopping`, SSEs, `norm_rg`, pLDDT, `consensus_level`), EBI InterPro, UniProtKB (batch `/accessions`), CATH v4.4 (fallback `v4_3_0`, Classes 1--4 and 6), and TED $\text{C}\alpha$ coordinates when sequence profiles fail.
2. **Residue-Set Boundary Overlap:** Expands discontinuous segments into residue sets $R_{\text{TED}}, R_{\text{Annot}}$ and computes exact $\text{IoU}$, $\text{Cov}_{\text{TED}}$, and $\text{Cov}_{\text{Annot}}$.
3. **5-Step CATH v4.4 Cascade:** Assigns `C.A.T.H` via official TED codes, `G3DSA` overlap, a 19,304-entry `Profile-to-CATH` bridge ($\ge 85\%$ fold purity), and 3D $\text{C}\alpha$/SSE rules [@Levitt1977SSE] (helix crossing angles, $\beta$-strand contact graphs, SVD $\beta$-sheet plane layer distances). Plain integers appear only where empirical precision exceeds 90%; unconfirmed levels use **`?` wildcards** (`2.60.40.?`, `3.?.?.?`), with geometric hints saved in `cath_fold_hint` and unmapped folds flagged as `novel_fold_candidate`.
4. **Strict `annotation_type` 4-Tier Hierarchy:** Assigns **Tier 1 (`specific_domain`, $W = 1.00$)** to curated `domain`/`repeat` hits, strict 1:1 `Pfam family` matches ($\text{IoU} \ge 0.65, \text{Cov}_{\text{Annot}} \ge 0.70, \text{Cov}_{\text{TED}} \ge 0.65$), and multi-lobe structural subdomains ($\text{Cov}_{\text{TED}} \ge 0.80, 0.20 \le \text{Cov}_{\text{Annot}} < 0.50$, suffixed as `N-terminal`/`central`/`C-terminal subdomain`); **Tier 2 (`superfamily`, $W = 0.85$)** using integrated InterPro names and curated unintegrated H-overrides weighted by real `G3DSA` overlap; **Tier 3 (`protein_family_context`, $W = 0.65$)** for whole-protein families (`HAMAP`, `NCBIfam family`, `PANTHER`) qualified by domain position and $L_{\text{domain}}/L_{\text{protein}}$; and **Tier 4 (`topology_fold`, $W = 0.45$)** strictly for confirmed topologies.
5. **Naming Hygiene, Synonym Unification & Quality Decoupling:** Groups synonyms under `IPR:{id}:{tier}` or `CATH_H:{code}`, strips `[gene]` tags and duplicate `domain` suffixes, enriches bare positional strings with parent family context, capitalizes titles while preserving biochemical prefixes (`p53`, `tRNA`, `von Willebrand`, `beta-`), penalizes residual prose ($0.40\times$), unintegrated PDB titles ($0.65\times$), and low-information strings (`DUF`, `Unknown domain`, $0.25\times$), and reports `structure_quality_flags` independently of `confidence_score`.
6. **Outputs:** Exports JSON, TSV/CSV, and Markdown via CLI (`inspect`, `name`) and Python API.

# Empirical Validation and Accuracy

We evaluated `ted-namer` (`v0.2.0`) across **100,000 randomly sampled TED domains** (`41,138` UniProt proteins; 10 deciles of 10,000 domains) and performed blind hold-out 5-fold cross-validation across all $N = 72,130$ official 4-level CATH domains (`34,443` UniProt proteins; \autoref{fig:validation}). Uncertainty is quantified via protein-clustered nonparametric bootstrap 95% confidence intervals ($B = 2,000$ resamples grouped by `uniprot_acc`).

![Empirical validation of `ted-namer` (`v0.2.0`) across 100,000 TED domains. (A) Distribution of naming specificity tiers across 41,138 UniProt proteins, resolving 86.96% (±0.27%) to specific domains/subdomains (Tier 1) or curated superfamilies (Tier 2). (B) Blind hold-out CATH accuracy across hierarchy levels ($N = 72,130$ domains across 34,443 proteins). Error bars denote 95% protein-clustered bootstrap confidence intervals ($B = 2,000$).\label{fig:validation}](figure1_validation.png)

Across 100,000 domains, **62.43% ± 0.38%** (`62,425`, decile range 61.4--63.3%, including `5,934` multi-lobe subdomains) received Tier 1 names (mean confidence $0.889$), **24.54% ± 0.32%** (`24,537`) received Tier 2 superfamily names ($0.769$), **12.86% ± 0.27%** (`12,858`) received Tier 3 family-context names ($0.422$), and **0.18% ± 0.04%** (`180`) fell to Tier 4/5 fallbacks. Enforcing hygiene rules eliminated bare positional names (`149` $\to$ `0`), bracketed gene tags (`703` $\to$ `0`), adjacent duplicate `domain` strings (`8` $\to$ `0`), and `T`-prefix collisions (`26,187` $\to$ `0`, with `10,740` = `10.74%` flagged `novel_fold_candidate`), reduced legacy lowercase starts from `5,054` to `744` biochemical prefixes, and restricted strict legacy PDB artifacts to **0.036%** (`36 / 100,000`). Leave-one-database-out testing (masking `Pfam`, $N = 36,843$) achieved **91.23%** exact normalized name concordance and **96.45%** semantic concordance.

| Validation Benchmark ($N = 72,130$) | Class ($C$) | Architecture ($C.A$) | Topology ($C.A.T$) | Superfamily ($C.A.T.H$) |
| :--- | :---: | :---: | :---: | :---: |
| **Confirmed Direct (Non-`?`) Precision** | **98.15 ± 0.11%** | **98.57 ± 0.09%** | **98.07 ± 0.11%** | **99.78 ± 0.05%** |
| **`Profile-to-CATH` Bridge (5-Fold CV)** | **99.28 ± 0.07%** | **98.63 ± 0.10%** | **97.93 ± 0.12%** | **92.29 ± 0.22%** |
| **Full Fallback Cascade** (Mask TED CATH) | **98.15 ± 0.11%** | **95.22 ± 0.17%** | **93.94 ± 0.19%** | **91.10 ± 0.23%** |

# Availability

`ted-namer` (`v0.2.0`) is MIT-licensed at [https://github.com/KarelBerka/TED-namer](https://github.com/KarelBerka/TED-namer), including the gzipped 100,000-domain benchmark (`examples/random_100000_results.tsv.gz`) and reproducible evaluation scripts (`scripts/reevaluate_100k.py`).

# Acknowledgements

Supported by the Ministry of Education, Youth and Sports of the Czech Republic (ELIXIR-CZ, LM2023055). AI-assisted coding tools (Google Antigravity) assisted with software implementation and documentation; all methodology, empirical experiments, and final code were verified by the author.

# References
