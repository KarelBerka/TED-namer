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

The expansion of the AlphaFold Protein Structure Database (AFDB) to over 214 million predicted structures [@Jumper2021AlphaFold; @Varadi2024AFDB] enabled the creation of **The Encyclopedia of Domains (TED)**, which catalogs over 365 million structural domains segmented across the protein universe [@Lau2024TED], and directly drove **CATH v4.4**, expanding CATH with roughly 90 million predicted TED domains and introducing Class 6 (`Special`) alongside hundreds of new architectures and superfamilies [@Waman2025CATHv44]. While TED provides high-accuracy continuous and discontinuous structural domain boundaries, consensus confidence (`high`, `medium`, `low`), and partial CATH structural assignments (`H`-level homologous superfamilies vs. `T`-level topologies), individual TED domains and clusters are identified purely by positional numeric tags (e.g., `AF-P04637-F1-model_v4_TED02`).

`ted-namer` is an open-source, asynchronous Python package and command-line tool that assigns standardized, human-readable biological names, specificity tiers, calibrated naming confidence scores, structural quality flags, and calibrated CATH v4.4 classifications (`C.A.T.H` with explicit `?` wildcards for unconfirmed deeper levels and `novel_fold_candidate` provenance) to individual TED domains and multi-member structural clusters. By combining exact residue-set boundary overlap across EBI InterPro [@Blum2025InterPro], Pfam [@PaysanLafosse2025Pfam], CDD [@Wang2023CDD], UniProtKB [@UniProt2025], and CATH v4.4 [@Waman2025CATHv44] with strict InterPro `annotation_type` semantics, curated unintegrated superfamily overrides, a 19,304-entry sequence-profile-to-CATH bridge, and direct 3D $\text{C}\alpha$ coordinate geometry hints, `ted-namer` resolves specific functional domain or homologous superfamily names for 86.8% of TED domains (99.8% when including positional family context) across a 100,000-domain benchmark.

# Statement of Need

Translating numeric TED identifiers into biologically valid domain names presents six persistent challenges that defeat naive sequence-database lookups:

1. **1D Sequence vs. 3D Structural Boundary Discrepancies:** Sequence profile Hidden Markov Models (HMMs) and 3D choppers such as Merizo [@Lau2023Merizo] and Chainsaw [@Wells2024Chainsaw] rarely share exact start/end endpoints, and many TED domains comprise discontinuous chain segments (e.g., `1-146_284-319`).
2. **Multi-Lobe Domain Splitting:** 3D choppers routinely split canonical sequence domains (such as eukaryotic protein kinases or Cytochrome P450 cores) into N-terminal and C-terminal structural lobes, causing standard Jaccard thresholds ($\text{IoU} \ge 0.35$) to miss both lobes because each covers only 25--45% of the full profile.
3. **Unintegrated Legacy CATH-Gene3D PDB Headers:** Many high-frequency `CATH-Gene3D` superfamilies remain unintegrated in InterPro or carry 1990s Protein Data Bank (PDB) crystal headers (such as `"Phosphorylase Kinase; domain 1"` or `"Transferase(Phosphotransferase) domain 1"`). Worse, naive scoring of `CATH-Gene3D` can let these Tier 2 labels outrank Tier 1 curated domain or subdomain names.
4. **`annotation_type` Conflation & Sentence/Positional Noise:** InterPro member databases such as `HAMAP` (100% `family`) and `NCBIfam` (98.1% `family`) annotate whole multi-domain proteins, while `PROSITE` (`conserved_site`) and UniProt features (`region`, bare `"N-terminal domain"`, or prose descriptions) frequently lack full domain scope or include bracketed gene symbols (e.g., `[bioU]`).
5. **Cross-Database Synonym Redundancy:** A single domain often matches `Pfam`, `CDD`, `SMART`, `InterPro`, `SUPERFAMILY`, and `CATH-Gene3D` under synonymous names that must be unified under their parent InterPro accession rather than double-counted.
6. **Hypothetical Novel Folds & `T`-Level Terminology Collision:** In TED [@Lau2024TED] and CATH v4.4 [@Waman2025CATHv44], `T` denotes a confirmed 3-digit *Topology* match, while domains lacking any structural/sequence homology (`cath_label = "-"`) represent hypothetical novel folds whose computationally predicted architectures/topologies should not be conflated with confirmed CATH `T`-level assignments.

`ted-namer` resolves all six issues reproducibly in one async, SQLite-cached CLI and Python library.

# State of the Field

Existing tools address parts of the domain annotation workflow in isolation:

- **3D Domain Choppers** (`Merizo` [@Lau2023Merizo], `Chainsaw` [@Wells2024Chainsaw]) segment protein structures into residue ranges, and **TED** [@Lau2024TED] aggregates consensus boundaries and Foldseek/Foldclass hits into `cath_label` codes, but does not synthesize human-readable domain or sub-lobe names from InterPro/UniProt.
- **Sequence Classification Resources** (`InterPro` [@Blum2025InterPro], `Pfam` [@PaysanLafosse2025Pfam], `CDD` [@Wang2023CDD]) annotate 1D sequence signatures, but do not reconcile discontinuous 3D structural choppings, multi-lobe splitting, legacy `G3DSA` titles, or cluster-wide consensus naming.
- **Structural Classification Resources** (`CATH v4.4` [@Waman2025CATHv44], `CATH-AlphaFlow` [@Waman2024CATH]) classify domains structurally into hierarchies (`C.A.T.H` including Class 6 `Special`), yet leave roughly 18--20% of TED domains unassigned without biological names or protein-family context fallbacks.

`ted-namer` bridges 3D structural segmentation and curated sequence nomenclature into one unified, provenance-aware pipeline.

# Methodology and Software Architecture

Given one or more TED domain IDs, [`TedNamingPipeline`](https://github.com/KarelBerka/TED-namer) executes six stages:

1. **Asynchronous Multi-Source Retrieval:** Concurrently queries TED (`chopping`, SSE counts, radius of gyration, pLDDT, `consensus_level`, raw CATH), EBI InterPro member databases (using `accessions` batch prefetching for UniProtKB), and the **CATH v4.4 REST API** (with automatic fallback to `v4_3_0`) covering Classes 1--4 and **Class 6 (`Special`)**, plus optional $\text{C}\alpha$ coordinates for unclassified domains.
2. **Residue-Set Boundary Overlap Engine:** Expands continuous and discontinuous segments into exact residue sets $R_{\text{TED}}$ and $R_{\text{Annot}}$, computing Jaccard $\text{IoU}$, $\text{Cov}_{\text{TED}}$, and $\text{Cov}_{\text{Annot}}$.
3. **Calibrated 5-Step CATH Classification Cascade & Novel Fold Flagging:** Assigns a four-slot hierarchical code (`C.A.T.H`) via: (i) official 4-level TED CATH, (ii) official 3-level TED topology completed at the $H$-level via high-overlap `G3DSA` or formatted `C.A.T.?`, (iii) high-IoU `CATH-Gene3D` (`G3DSA`) overlap, (iv) a 19,304-entry high-purity ($\ge 85\%$) `Profile-to-CATH` bridge (`Pfam`, `SUPERFAMILY`, `CDD`, `SMART`, `InterPro`), and (v) direct 3D $\text{C}\alpha$ coordinate geometry (Levitt--Greer SSE extraction [@Levitt1977SSE], helix crossing angles, strand contact graphs, and SVD $\beta$-sheet layer analysis). Plain integer levels are emitted **only when empirical out-of-sample precision exceeds 90%** (`known_depth`); deeper unconfirmed levels use explicit **`?` wildcards** (e.g., `2.60.40.?` for `T`-level matches or `3.?.?.?` for Class 3 novel fold candidates), while full geometric predictions are exported separately in `cath_fold_hint` alongside the boolean `novel_fold_candidate` flag.
4. **4-Tier Specificity Classification With Strict `annotation_type` Semantics:** Categorizes hits into **Tier 1 (`specific_domain`, $W = 1.00$)** strictly for curated `domain`/`repeat` entries (including `NCBIfam` `domain` hits), strict 1:1 `Pfam` `family` matches ($\text{IoU} \ge 0.65$, $\text{Cov}_{\text{Annot}} \ge 0.70$, $\text{Cov}_{\text{TED}} \ge 0.65$), and multi-lobe structural subdomains ($\text{Cov}_{\text{TED}} \ge 0.80$, $0.20 \le \text{Cov}_{\text{Annot}} < 0.50$, suffixed as `N-terminal`, `central`, or `C-terminal subdomain`); **Tier 2 (`superfamily`, $W = 0.85$)** honoring curated integrated InterPro titles and curated overrides for unintegrated legacy CATH H-names without artificial `cath_h` overlap inflation; **Tier 3 (`protein_family_context`, $W = 0.65$)** for whole-protein families (`HAMAP`, `NCBIfam` family, `PANTHER`, `InterPro family`, `conserved_site` motifs, and parent UniProt `SIMILARITY` comments qualified by relative position and scaled by $L_{\text{domain}} / L_{\text{protein}}$); and **Tier 4 (`topology_fold`, $W = 0.45$)** strictly for confirmed 3-digit CATH topologies (`?`-free).
5. **Naming Hygiene, Synonym Unification & Structural Quality Separation:** Unifies cross-database hits under `IPR:{id}:{tier}` or `CATH_H:{code}`, strips bracketed gene/EC tags (`[bioU]`) and duplicate `domain` suffixes, enriches bare positional labels (`"N-terminal domain"`) with parent family context (`"argonaute N-terminal domain"`), penalizes residual sentence-like prose ($0.40\times$), residual legacy PDB headers ($0.65\times$), and low-information strings (`DUF`, `Unknown domain`, $0.25\times$), and decouples structural quality warnings (`structure_quality_flags`: `ok`, `low_plddt`, `low_consensus`, `non_globular`) from annotation `confidence_score`.
6. **Consensus Synthesis:** Outputs JSON, TSV/CSV, and Markdown reports via the `ted-namer` CLI (`inspect`, `name`, `batch`) or Python API.

# Empirical Validation and Accuracy

We evaluated the complete pipeline across **100,000 randomly sampled TED domains** (10 deciles of 10,000 unique domains spanning all kingdoms of life) and performed blind hold-out 5-fold cross-validation on all $N = 72,130$ domains carrying official 4-level CATH assignments (\autoref{fig:validation}).

![Empirical validation of `ted-namer` across 100,000 randomly sampled TED domains (10 deciles of 10,000 domains). (A) Distribution of assigned naming specificity tiers under strict `annotation_type` semantics, with 86.8% resolved to a specific domain/subdomain (Tier 1) or curated homologous superfamily (Tier 2). (B) Blind hold-out CATH classification accuracy across hierarchy levels ($N = 72,130$) for confirmed direct (non-wildcard) levels, 5-fold cross-validation of the Profile-to-CATH bridge, and the full fallback cascade.\label{fig:validation}](figure1_validation.png)

Across all 100,000 domains, **62.14%** (`62,140`, decile range 61.1--63.2%, including `5,901` multi-lobe subdomains) received a Tier 1 specific domain or subdomain name (mean confidence $0.884$), **24.69%** (`24,692`) received a Tier 2 homologous superfamily name (mean confidence $0.765$), **12.99%** (`12,989`) received a Tier 3 positional protein family name (mean confidence $0.420$), and **0.18%** (`179`) required a Tier 4 confirmed-topology fold (`69`) or uncharacterized (`110`) fallback. Enforcing naming hygiene and curated unintegrated H-level overrides completely eliminated bare positional names (`149` $\to$ `0`), bracketed gene/EC tags (`703` $\to$ `0`), adjacent double-domain strings (`8` $\to$ `0`), and `T`-prefix notation collisions (`26,187` $\to$ `0`, cleanly represented as `?` wildcards with `10,740` = `10.74%` flagged as `novel_fold_candidate`), while reducing sentence-like descriptions by **96.7%** (`2,454` $\to$ `82`) and legacy PDB header noise to **0.036%** (`36 / 100,000`). In leave-one-database-out testing ($N = 36,511$ domains where `Pfam` was held out), `ted-namer` achieved **91.1%** exact normalized name concordance and **96.4%** semantic keyword concordance.

| Validation Benchmark ($N = 72,130$) | Class ($C$) | Architecture ($C.A$) | Topology ($C.A.T$) | Superfamily ($C.A.T.H$) |
| :--- | :---: | :---: | :---: | :---: |
| **Confirmed Direct (Non-`?`) Precision** | **98.2%** | **98.6%** | **98.1%** | **99.8%** |
| **`Profile-to-CATH` Bridge (5-Fold CV)** | **99.3%** | **98.6%** | **97.9%** | **92.3%** |
| **Full Fallback Cascade** (Mask TED CATH) | **98.2%** | **95.2%** | **93.9%** | **91.1%** |

# Availability

`ted-namer` is available under the MIT License at [https://github.com/KarelBerka/TED-namer](https://github.com/KarelBerka/TED-namer), including the refreshed 100,000-domain benchmark dataset (`examples/random_100000_results.tsv`), summary statistics (`examples/benchmark_100000_summary.json`), and reproducible validation scripts (`scripts/benchmark_100000_random.py` and `scripts/reevaluate_100k.py`).

# Acknowledgements

This work was supported by the Ministry of Education, Youth and Sports of the Czech Republic (ELIXIR-CZ, LM2023055). AI-assisted coding tools (Google Antigravity) were used during software implementation, benchmarking, and documentation drafting; all methodology, literature validation criteria, empirical experiments, and final code were reviewed and verified by the author.

# References
