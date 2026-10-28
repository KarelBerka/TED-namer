---
title: 'ted-namer: Automated biological naming and calibrated CATH classification for structural domains in The Encyclopedia of Domains (TED)'
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
date: 27 October 2026
bibliography: paper.bib
---

# Summary

The expansion of the AlphaFold Protein Structure Database (AFDB) to over 214 million predicted structures [@Jumper2021AlphaFold; @Varadi2024AFDB] enabled the creation of **The Encyclopedia of Domains (TED)**, which catalogs more than 365 million structural domains segmented across the protein universe [@Lau2024TED]. While TED provides high-accuracy continuous and discontinuous structural domain boundaries alongside partial CATH fold assignments, individual TED domains and structural clusters are identified solely by positional numeric tags (e.g., `AF-P04637-F1-model_v4_TED02`).

`ted-namer` is an open-source, asynchronous Python package and command-line tool that assigns standardized, human-readable biological names, specificity tiers, calibrated confidence scores, and four-level CATH structural classifications (`C.A.T.H`) to individual TED domains and multi-member structural clusters. By integrating exact residue-set boundary overlap across EBI InterPro [@PaysanLafosse2023InterPro], Pfam [@Mistry2021Pfam], CDD [@Wang2023CDD], UniProtKB [@UniProt2023], and CATH [@Sillitoe2021CATH; @Waman2024CATH] with a 19,304-entry sequence-profile-to-CATH bridge and direct 3D $\text{C}\alpha$ coordinate geometry rules, `ted-namer` resolves specific domain or homologous superfamily names for 89.7% of TED domains across a 100,000-domain benchmark and achieves 100% CATH classification coverage with calibrated uncertainty notation.

# Statement of Need

Structural biologists, evolutionary bioinformaticians, and database curators routinely inspect TED domain decompositions or cluster millions of AlphaFold domains by structural similarity. However, translating numeric TED identifiers into biologically meaningful domain names presents five persistent challenges that prevent naive coordinate matching:

1. **1D Sequence vs. 3D Structural Boundary Discrepancies:** Sequence profile Hidden Markov Models (HMMs) and 3D structural choppers such as Merizo [@Lau2023Merizo] and Chainsaw [@Wells2024Chainsaw] rarely share identical start and end coordinates. Furthermore, 16% of TED domains consist of discontinuous peptide segments (e.g., `1-146_284-319`) that fold into a single 3D core.
2. **Multi-Lobe Domain Splitting:** 3D choppers frequently divide canonical sequence domains (such as eukaryotic protein kinases or ABC transporter cores) into separate N-terminal and C-terminal structural lobes. Standard Jaccard overlap thresholds reject both lobes because each covers only 30--45% of the parent sequence domain.
3. **Missing CATH H-Level Names and Legacy PDB Header Artifacts:** In the CATH REST API, four-digit Homologous Superfamily (`H-level`) names are frequently null or populated with raw 1990s Protein Data Bank (PDB) header strings from the first crystal structure solved in that fold (e.g., naming bacterial methyltransferases *"Vaccinia Virus protein VP39"*).
4. **Cross-Database Synonym Redundancy:** A single structural domain often matches multiple member databases (`Pfam`, `CDD`, `SMART`, `PROSITE`, `SUPERFAMILY`, `CATH-Gene3D`) under synonymous titles that must be unified rather than double-counted.
5. **Unclassified Structural Domains (`cath_label = "-"`):** Approximately 18--20% of TED domains lack an official CATH assignment, requiring automated fold inference with clear provenance separating confirmed classifications from computational estimates.

`ted-namer` addresses these challenges in a single reproducible pipeline with SQLite disk caching, rich CLI reporting, and programmatic Python APIs.

# State of the Field

Existing tools address parts of the domain annotation workflow in isolation:

- **3D Domain Choppers** (`Merizo` [@Lau2023Merizo], `Chainsaw` [@Wells2024Chainsaw]) segment protein structures into residue ranges and **TED** [@Lau2024TED] aggregates their consensus boundaries, but neither generates human-readable functional domain or sub-lobe names.
- **Sequence Annotation Platforms** (`InterProScan` [@PaysanLafosse2023InterPro], `CDD` [@Wang2023CDD], `Pfam` [@Mistry2021Pfam]) annotate 1D protein sequences with profile HMMs, but do not reconcile discontinuous 3D structural choppings, resolve split structural lobes, or synthesize consensus names across multi-protein structural clusters.
- **Structural Classification Resources** (`CATH-AlphaFlow` [@Waman2024CATH]) assign structural fold codes, but leave ~20% of TED domains unassigned and do not bridge curated InterPro superfamily nomenclature with sequence-domain and protein-family positional fallbacks.

`ted-namer` bridges 3D structural segmentation and 1D functional curation into a unified, boundary-aware nomenclature engine.

# Methodology and Software Architecture

Given one or more TED domain IDs, [`TedNamingPipeline`](https://github.com/KarelBerka/TED-namer) executes six stages:

1. **Asynchronous Multi-Source Retrieval:** Concurrently fetches domain choppings, secondary structure element (SSE) counts, radius of gyration, pLDDT, and raw CATH labels from the TED API; all overlapping member database hits from EBI InterPro; recommended protein names, Swiss-Prot status, `SIMILARITY` comments, and `DOMAIN` features from UniProtKB; and $\text{C}\alpha$ coordinates for unclassified domains.
2. **Residue-Set Boundary Overlap Engine:** Expands continuous and discontinuous domain segments into exact integer residue sets $R_{\text{TED}}$ and $R_{\text{Annot}}$, computing Jaccard Intersection-over-Union ($\text{IoU} = |R_{\text{TED}} \cap R_{\text{Annot}}| / |R_{\text{TED}} \cup R_{\text{Annot}}|$), $\text{Cov}_{\text{TED}}$, and $\text{Cov}_{\text{Annot}}$.
3. **Calibrated 5-Step CATH Classification Cascade:** Assigns a four-level `C.A.T.H` code via: (i) official 4-level TED CATH, (ii) official 3-level TED topology completed at the $H$-level via high-overlap `G3DSA`, (iii) high-IoU `CATH-Gene3D` (`G3DSA`) overlap, (iv) a bundled lookup table of **19,304 high-purity ($\ge 85\%$ fold purity) sequence profiles** (`Pfam`, `SUPERFAMILY`, `CDD`, `SMART`, `InterPro`) mapped to CATH, and (v) direct 3D $\text{C}\alpha$ coordinate geometry (Levitt--Greer SSE extraction [@Levitt1977SSE], inter-helical crossing angles, $\beta$-strand contact graph closure, and Singular Value Decomposition of $\beta$-sheet planes for 3-layer Rossmann vs. 2-layer plait folds). Hierarchy levels are emitted as plain numbers **only when empirical precision exceeds 90%** (`known_depth`); lower-confidence heuristic levels are prefixed with `T` (e.g., `3.40.50.T300` or `3.T40.T50.T10`).
4. **4-Tier Specificity Classification & Multi-Lobe Subdomain Resolution:** Categorizes annotations into **Tier 1 (`specific_domain`, $W = 1.00$)** for direct domain matches or split structural lobes ($\text{Cov}_{\text{TED}} \ge 0.80$, $0.20 \le \text{Cov}_{\text{Annot}} < 0.50$, automatically suffixed as `N-terminal`, `central`, or `C-terminal subdomain`); **Tier 2 (`superfamily`, $W = 0.85$)** prioritizing curated integrated InterPro superfamily titles over raw `G3DSA` PDB headers; **Tier 3 (`protein_family_context`, $W = 0.65$)** combining parent protein family with spatial position scaled by $L_{\text{domain}} / L_{\text{protein}}$; and **Tier 4 (`topology_fold`, $W = 0.45$)** for fold-only fallbacks.
5. **InterPro Synonym Unification & Calibrated Scoring:** Unifies synonymous cross-database hits under `IPR:{id}:{tier}`, selects the cleanest title via curation-weighted voting, applies a $0.25\times$ penalty to low-information strings (`DUF`, `Uncharacterized`), and computes a calibrated confidence score $S(c) \in [0, 1]$ tailored to single-domain ($N = 1$) or multi-member cluster ($N > 1$) queries.
6. **Consensus Synthesis:** Outputs structured JSON, TSV, and Markdown reports via the `ted-namer` CLI (`inspect`, `name`, `batch`) or Python API.

# Empirical Validation and Accuracy

We evaluated `ted-namer` on **100,000 randomly sampled TED domains** across 100 batches of 1,000 domains (10 consecutive deciles of 10,000 unique domains spanning all kingdoms of life) and performed blind hold-out cross-validation on the $N = 72,130$ domains carrying official 4-level CATH assignments (\autoref{fig:validation}).

![Empirical validation of `ted-namer` across 100,000 randomly sampled TED domains (100 batches of 1,000 domains). (A) Distribution of assigned biological naming specificity tiers, with 89.7% of domains resolved to a specific functional domain/subdomain (Tier 1) or homologous superfamily (Tier 2). (B) Blind hold-out CATH classification accuracy across hierarchy levels ($N = 72,130$) for non-`T` direct predictions, 5-fold cross-validation of the Profile-to-CATH bridge, and the full fallback cascade.\label{fig:validation}](figure1_validation.png)

Across all 100,000 domains, **62.6%** (62,560; decile range 61.5--64.7%, including 6,022 multi-lobe subdomains) received a Tier 1 specific domain or subdomain name (mean confidence $0.899$), **27.1%** (27,098; decile range 23.9--28.4%) received a Tier 2 homologous superfamily name (mean confidence $0.809$), **10.2%** (10,160; decile range 8.1--11.7%) received a Tier 3 positional protein family name (mean confidence $0.410$), and **0.18%** (182, predominantly unreviewed TrEMBL proteins annotated as *"Uncharacterized protein"*) received a Tier 4 fold-only fallback (mean confidence $0.235$). Combined Tier 1 + Tier 2 resolution remained tightly bounded between **88.3% and 90.8%** across all ten 10,000-domain deciles (**89.7%** overall). Prioritizing integrated InterPro superfamily names over raw `G3DSA` headers restricted legacy PDB title artifacts to **0.002%** (2 / 100,000 domains, a $>99.9\%$ reduction). In leave-one-database-out testing ($N = 34,648$ domains where `Pfam` was held out), `ted-namer` achieved **90.1%** exact normalized name concordance (`31,201 / 34,648`) and **95.7%** semantic keyword concordance (`33,156 / 34,648`).

| Validation Benchmark ($N = 72,130$) | Class ($C$) | Architecture ($C.A$) | Topology ($C.A.T$) | Superfamily ($C.A.T.H$) |
| :--- | :---: | :---: | :---: | :---: |
| **Non-`T` Direct Precision** | **98.2%** | **98.6%** | **98.1%** | **99.8%** |
| **`Profile-to-CATH` Bridge (5-Fold CV)** | **99.3%** | **98.6%** | **97.9%** | **92.3%** |
| **Full Fallback Cascade** (Mask TED CATH) | **98.2%** | **95.2%** | **93.9%** | **91.1%** |

# Availability

`ted-namer` is available under the MIT License at [https://github.com/KarelBerka/TED-namer](https://github.com/KarelBerka/TED-namer), including the full 100,000-domain benchmark dataset (`examples/random_100000_results.tsv`) and reproducible accuracy evaluation suite (`scripts/benchmark_100000_random.py`).

# Acknowledgements

This work was supported by the Ministry of Education, Youth and Sports of the Czech Republic (ELIXIR-CZ, LM2023055). AI-assisted coding tools (Google Antigravity) were used during software implementation, benchmarking, and documentation drafting; all methodology, empirical validation experiments, and final code were reviewed and verified by the author.

# References
