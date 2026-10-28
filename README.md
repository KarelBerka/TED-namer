# `ted-namer`: Structural Domain & Cluster Namer for TED

[![CI](https://github.com/KarelBerka/TED-namer/actions/workflows/ci.yml/badge.svg)](https://github.com/KarelBerka/TED-namer/actions/workflows/ci.yml)
[![JOSS Paper Draft](https://github.com/KarelBerka/TED-namer/actions/workflows/joss-paper.yml/badge.svg)](https://github.com/KarelBerka/TED-namer/actions/workflows/joss-paper.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Code style: uv](https://img.shields.io/badge/managed%20by-uv-orange.svg)](https://github.com/astral-sh/uv)

Generates standardized, biologically meaningful human-readable names and CATH structural classifications for domains and domain clusters from **TED (The Encyclopedia of Domains, [ted.cathdb.info](https://ted.cathdb.info))**.

Detailed documentation of the entire naming methodology is available in [**`docs/NAMING_PROCESS.md`**](docs/NAMING_PROCESS.md).

---

## Key Features

- **Boundary-Aware Overlap Engine**: Computes exact residue-level Intersection-over-Union (IoU), TED coverage, and annotation coverage across both continuous (`109-287`) and discontinuous (`15-65_140-210`) TED segments against **InterPro, Pfam, CDD, SMART, PROSITE, PANTHER, SUPERFAMILY, CATH-Gene3D, and UniProt features**.
- **CATH H-Level & T-Level Resolution**: Prioritizes curated 4-digit Homologous Superfamily names ($H$-level, e.g., *"p53/RUNT-type transcription factor DNA-binding domain"*) resolved via InterPro CATH-Gene3D (`G3DSA`) and CATH FunFams, rather than settling for generic 3-digit Topology names ($T$-level).
- **Direct 3D Structural Rules & CATH Classification Engine**:
  - Class ($C$): Determined from secondary structure content ($\alpha$, $\beta$).
  - Architecture ($A$) & Topology ($T$): Computed directly from AlphaFold C-alpha 3D coordinates (helix crossing angles $\theta_{ij}$, beta-strand adjacency & barrel closure, signed distance to beta-sheet planes for Rossmann 3-layer vs 2-layer).
  - Standardized **`T`-prefix notation** (e.g. `3.40.50.T10`, `1.10.T10.T10`) distinguishing experimentally curated levels from computational estimations.
- **4-Tier Specificity Hierarchy**:
  - **Tier 1 (Specific Domain)**: Functional domain profiles (CDD, Pfam, InterPro domain, UniProt sequence feature).
  - **Tier 2 (Superfamily)**: Broad homologous superfamily (CATH H-level, SUPERFAMILY, PANTHER).
  - **Tier 3 (Protein Family Context)**: Positional family context (e.g., *BioU C-terminal domain*) for uncharacterized or supra-domain regions.
  - **Tier 4 (Topology / Fold)**: Structural fold (e.g., *Immunoglobulin-like fold*) as a fallback when sequence homology is absent.
- **Novel Fold & Unannotated Fallbacks**: Clean handling of unreviewed TrEMBL proteins and novel structural folds (`cath_label = "-"`) using parent UniProt `SIMILARITY` comments, GO terms, and relative domain position (`N-terminal`, `central`, `C-terminal`).
- **High-Performance Async Pipeline**: SQLite-backed caching of API calls and PDB models with concurrent batch processing reaching >20 domains/second.

---

## Installation

### Via `uv` (Recommended)

```bash
# Clone repository
git clone https://github.com/KarelBerka/TED-namer.git
cd TED-namer

# Install dependencies and sync virtual environment
uv sync --extra dev
```

### Via `pip`

```bash
pip install .
```

---

## Quickstart

### 1. Inspect All TED Domains for a UniProt Protein

Retrieve all AlphaFold-predicted TED domains, residue boundaries, and consensus CATH labels for a UniProt accession (e.g., Human p53, `P04637`):

```bash
uv run ted-namer inspect P04637
```

Output:
```
TED Domains for UniProt P04637:
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┳━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┓
┃ TED ID                     ┃ Range     ┃ Length ┃ pLDDT ┃ CATH Label  ┃ SSE (H/E/C)    ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━╇━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━┩
│ AF-P04637-F1-model_v4_TED01 │ 94-293    │    200 │ 92.4  │ 2.60.40.720 │ 12% / 38% / 50%│
│ AF-P04637-F1-model_v4_TED02 │ 324-356   │     33 │ 84.1  │ 1.20.1280.10│ 85% / 0% / 15% │
└────────────────────────────┴───────────┴────────┴───────┴─────────────┴────────────────┘
```

### 2. Name a Single TED Domain or Homologous Cluster

Name a cluster of orthologous/paralogous TED domains (e.g., DNA-binding domains of p53, p63, and p73):

```bash
uv run ted-namer name \
  --ted-ids AF-P04637-F1-model_v4_TED01,AF-Q9H3D4-F1-model_v4_TED01,AF-O15350-F1-model_v4_TED01 \
  --cluster-id p53_family_DBD \
  -o results.json -o results.tsv -o report.md
```

### 3. Batch Naming from a Table

Provide a TSV or CSV with columns `cluster_id` and `ted_id`:

```bash
uv run ted-namer batch \
  --input examples/random_1000_queries.tsv \
  --concurrency 25 \
  --output batch_results.tsv \
  --report batch_report.md
```

---

## Python API Usage

```python
import asyncio
from ted_namer.pipeline import TedNamingPipeline

async def main():
    pipeline = TedNamingPipeline(use_cache=True, concurrency=20)

    # Name domains grouped by cluster ID
    clusters = {
        "p53_DBD": ["AF-P04637-F1-model_v4_TED01"],
        "p53_tetramer": ["AF-P04637-F1-model_v4_TED02"],
    }
    results = await pipeline.name_clusters(clusters)

    for res in results:
        print(f"Cluster: {res.cluster_id}")
        print(f"  Proposed Name:  {res.proposed_name}")
        print(f"  Short Label:    {res.short_label}")
        print(f"  Tier:           {res.specificity_tier.value}")
        print(f"  Confidence:     {res.confidence_score:.3f}")
        for m in res.members:
            print(f"  Member: {m.domain.ted_id} -> CATH: {m.domain.cath_classification_display}")

if __name__ == "__main__":
    asyncio.run(main())
```

---

## The 4-Tier Specificity Hierarchy

| Tier | Name | Typical Sources | Example |
| :---: | :--- | :--- | :--- |
| **1** | **Specific Domain** | CDD, Pfam, InterPro domain, UniProt feature | *p53 DNA-binding domain* |
| **2** | **Superfamily** | CATH H-level (Gene3D), SUPERFAMILY, PANTHER | *p53/RUNT transcription factor superfamily* |
| **3** | **Protein Family Context** | Parent UniProt `SIMILARITY` family + position | *BioU C-terminal domain* |
| **4** | **Topology / Fold** | CATH T-level, 3D structural rules | *Immunoglobulin-like fold* |

---

## CATH Classification & 3D Structural Rules

When official CATH classifications are incomplete, `ted-namer` applies 3D coordinate rules directly to AlphaFold structures:
- **Class $C$**: From $\% \alpha$ and $\% \beta$ secondary structure content.
- **Architecture $A$**:
  - Orthogonal bundle (`1.10`) vs. Up-down bundle (`1.20`) via helix crossing angles $\theta_{ij}$.
  - Beta-sandwich (`2.60`) vs. Beta-barrel (`2.40`) via strand adjacency graphs.
  - Rossmann 3-layer $\alpha/\beta/\alpha$ (`3.40`) vs. 2-layer $\alpha/\beta$ (`3.30`) via signed distance to sheet planes.
- **Estimated Levels**: Prefixed with `T` (e.g., `3.40.50.T10`, `1.10.T10.T10`).

### Examples of Domains with Unknown / Inferred CATH (`cath_label = "-"`)

In TED, ~18–20% of structural domains lack an official CATH classification. `ted-namer` resolves their 3D fold geometry and biological names seamlessly:

| Domain ID | Protein & Taxon | TED Raw CATH | Assigned CATH Code | Detection Mechanism | Assigned Name & Tier |
| :--- | :--- | :---: | :---: | :--- | :--- |
| `AF-A4F4K3-F1-model_v4_TED02` | Cytochrome P450 CYP1A1 | `-` | **`1.10.630.T10`** | Sequence overlap with `G3DSA:1.10.630` + multi-lobe subdomain | **Cytochrome P450 central subdomain** *(Tier 1, conf: 0.772)* |
| `AF-Q8U3D2-F1-model_v4_TED01` | Argonaute ago (*Pyrococcus*) | `-` | **`3.90.70.T180`** | `Profile-to-CATH` bridge (`known_depth = 3`) | **N-terminal domain** *(Tier 1, conf: 0.820)* |
| `AF-O05442-F1-model_v4_TED01` | CpnT (*Bacillus*) | `-` | **`1.T25.T40.T10`** | 3D helix curvature ($\alpha$-solenoid rule) | **NTD domain** *(Tier 1, conf: 0.787)* |
| `AF-P0ABG4-F1-model_v4_TED02` | FtsW (*E. coli*) | `-` | **`3.T30.T70.T10`** | 3D single-sided helix-sheet packing | **Cell cycle protein domain** *(Tier 1, conf: 0.605)* |
| `AF-A0A1G9Z131-F1-model_v4_TED02` | BioU (*Corynebacterium*) | `-` | **`3.T40.T50.T720`** | Low-IoU `Profile-to-CATH` / 3D Rossmann rule | **(S)-8-amino-7-oxononanoate synthase BioU [bioU] domain** *(Tier 1, conf: 0.597)* |
| `AF-P0A710-F1-model_v4_TED01` | YciB (*E. coli*) | `1.20.58` | **`1.20.58.T10`** | TED 3-level topology + estimated $H$ | **Intracellular septation protein A domain** *(Tier 1, conf: 0.656)* |

*Full case studies and geometric derivations are detailed in [docs/NAMING_PROCESS.md](docs/NAMING_PROCESS.md).*

---

## Blind Cross-Validation & Ground-Truth Accuracy (100,000 TED Domains)

Detailed ground-truth validation across **100,000 unique TED domains** (100 batches of 1,000 domains / 10 deciles of 10,000 domains) is documented in [**`docs/ACCURACY_VALIDATION.md`**](docs/ACCURACY_VALIDATION.md) and reproducible via [`scripts/benchmark_100000_random.py`](scripts/benchmark_100000_random.py):

### 1. Blind Hold-Out CATH Classification Accuracy ($N = 72,130$ Official CATH Domains)

| Hold-Out Scenario | Class ($C$) | Architecture ($C.A$) | Topology ($C.A.T$) | Superfamily ($C.A.T.H$) |
| :--- | :---: | :---: | :---: | :---: |
| **Full Fallback Cascade** *(Mask TED `cath_label`)* | **98.2%** | **95.2%** | **93.9%** | **91.1%** |
| **`Profile-to-CATH` Bridge (5-Fold CV)** *(Mask TED CATH & `G3DSA`; 87.9% cov)* | **99.3%** | **98.6%** | **97.9%** | **92.3%** |
| **Pure 3D/SSE Rules Only** *(Zero sequence/profile homology)* | **85.6%** *(96.7% pure $\alpha/\beta$)* | **41.5%** *(marked `T`)* | **25.4%** *(marked `T`)* | — |

- **Precision of Direct (Non-`T`-Prefixed) Levels ($N = 72,130$)**: **98.2%** for $C$, **98.6%** for $C.A$, **98.1%** for $C.A.T$, and **99.8%** for $C.A.T.H$.
- **Leave-One-Database-Out Naming Concordance ($N = 34,648$)**: **90.1%** exact InterPro/normalized name concordance (`31,201 / 34,648`) and **95.7%** semantic keyword concordance (`33,156 / 34,648`) when holding out `Pfam`.

---

## Large-Scale Benchmark & Statistical Stability (100,000 TED Domains)

Tested across **100,000 unique TED domains** (100 batches of 1,000 domains / 10 deciles of 10,000 domains across all kingdoms of life):

| Metric | Range Across 10 Deciles (`10x 10k`) | **Total Across 100,000 TED Domains** |
| :--- | :---: | :---: |
| **Domain Count** | 10,000 per decile | **100,000** |
| **Mean Calibrated Confidence Score** | 0.818 – 0.843 | **0.823** *(Tier 1: 0.899, Tier 2: 0.809, Tier 3: 0.410)* |
| **Tier 1: Specific Domain & Subdomain** | 61.5% – 64.7% | **62.6%** (62,560 / 100,000; incl. 6,022 subdomains) |
| **Tier 2: Superfamily (Curated IPR)** | 23.9% – 28.4% | **27.1%** (27,098 / 100,000) |
| **Tier 1 + Tier 2 (Domain or Superfamily)** | **88.3% – 90.8%** | **89.7%** (89,658 / 100,000) |
| **Tier 3: Protein Family Context** | 8.1% – 11.7% | **10.2%** (10,160 / 100,000) |
| **Tier 4: Topology / Fold** | 0.00% – 1.15% | **0.18%** (182 / 100,000) |
| **CATH: Official (4-level)** | 69.8% – 73.2% | **72.1%** (72,130 / 100,000) |
| **CATH: TED T-level (`ted_t_plus_*`)** | 8.7% – 11.6% | **9.7%** (9,692 / 100,000) |
| **CATH: Profile-to-CATH & Gene3D Overlap** | 6.7% – 8.7% | **7.4%** *(4.3% Profile-to-CATH + 3.1% Gene3D)* |
| **CATH: Pure 3D/SSE Rules (`C.T...`)** | 9.6% – 11.4% | **10.7%** (10,740 / 100,000) |
| **Low pLDDT (< 60)** | 0.4% – 1.2% | **0.7%** (686 / 100,000) |

Full benchmark reports and datasets are available in the [`examples/`](examples/) directory.

---

## Running Tests

```bash
uv run --extra dev pytest
```

---

## Contributing & Community Guidelines

Contributions, bug reports, and feature requests are welcome!
- See [**`CONTRIBUTING.md`**](CONTRIBUTING.md) for development setup, testing instructions, and pull request guidelines.
- Participation is governed by our [**`CODE_OF_CONDUCT.md`**](CODE_OF_CONDUCT.md).
- Report issues or edge cases at [https://github.com/KarelBerka/TED-namer/issues](https://github.com/KarelBerka/TED-namer/issues).

---

## Citation & Manuscript

If you use `ted-namer` in your research, please cite our manuscript in [`paper/paper.md`](paper/paper.md) (or via [`CITATION.cff`](CITATION.cff)):

```bibtex
@software{Berka2026tednamer,
  author    = {Berka, Karel},
  title     = {{ted-namer}: Automated biological naming and calibrated {CATH} classification for structural domains in {The Encyclopedia of Domains (TED)}},
  year      = {2026},
  url       = {https://github.com/KarelBerka/TED-namer},
  version   = {0.1.0}
}
```

- **Repository**: [https://github.com/KarelBerka/TED-namer](https://github.com/KarelBerka/TED-namer)
- **Author**: Karel Berka ([ORCID: 0000-0001-9472-2589](https://orcid.org/0000-0001-9472-2589)), Palacký University Olomouc
- **License**: [MIT](LICENSE)
