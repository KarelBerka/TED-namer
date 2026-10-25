# `ted-namer`: Structural Domain & Cluster Namer for TED

[![CI](https://github.com/KarelBerka/TED-namer/actions/workflows/ci.yml/badge.svg)](https://github.com/KarelBerka/TED-namer/actions/workflows/ci.yml)
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

In TED, over 20% of structural domains lack an official CATH classification. `ted-namer` resolves their 3D fold geometry and biological names seamlessly:

| Domain ID | Protein & Taxon | TED Raw CATH | Assigned CATH Code | Detection Mechanism | Assigned Name & Tier |
| :--- | :--- | :---: | :---: | :--- | :--- |
| `AF-A0A1G9Z131-F1-model_v4_TED02` | BioU (*Corynebacterium*) | `-` | **`3.40.50.T10`** | 3D SVD sheet plane + helix distances | **BioU C-terminal domain** *(Tier 3)* |
| `AF-Q8U3D2-F1-model_v4_TED01` | Argonaute ago (*Pyrococcus*) | `-` | **`3.20.20.T10`** | 3D cylindrical closure of 8 $\beta/\alpha$ | **N-terminal domain** *(Tier 1)* |
| `AF-O05442-F1-model_v4_TED01` | CpnT (*Bacillus*) | `-` | **`1.25.40.T10`** | 3D helix curvature (solenoid) | **NTD domain** *(Tier 1)* |
| `AF-P0ABG4-F1-model_v4_TED02` | FtsW (*E. coli*) | `-` | **`3.30.70.T10`** | 3D single-sided helix-sheet packing | **SEDS central domain** *(Tier 3)* |
| `AF-A4F4K3-F1-model_v4_TED02` | Cytochrome P450 CYP1A1 | `-` | **`1.T10.T630.T10`** | Sequence overlap with `G3DSA:1.10.630` | **Cytochrome P450** *(Tier 2)* |
| `AF-P0A710-F1-model_v4_TED01` | YciB (*E. coli*) | `1.20.58` | **`1.20.58.T10`** | TED partial topology + estimated $H$ | **YciB N-terminal domain** *(Tier 3)* |

*Full case studies and geometric derivations are detailed in [docs/NAMING_PROCESS.md](docs/NAMING_PROCESS.md).*

---

## Large-Scale Benchmark & Stability (10,000 TED Domains)

Tested across **10,000 diverse TED domains** (10 independent batches of 1,000 domains each) spanning 27 taxonomic kingdoms, model organisms, and unreviewed TrEMBL sets:

| Metric | Batches 1–5 (1–5k) | Batches 6–10 (6–10k) | **TOTAL 10,000 TED Domains** |
| :--- | :---: | :---: | :---: |
| **Domain Count** | 5,000 | 5,000 | **10,000** |
| **Mean Confidence Score** | 0.844 | 0.861 | **0.852** |
| **Tier 1: Specific Domain** | 53.9% | 57.1% | **55.5%** |
| **Tier 2: Superfamily** | 24.5% | 25.2% | **24.9%** |
| **Tier 3: Protein Family Context** | 21.5% | 17.6% | **19.6%** |
| **Tier 4: Topology / Fold** | 0.1% | 0.0% | **0.0%** (3 / 10,000) |
| **CATH: Official (4-level)** | 67.0% | 72.5% | **69.8%** |
| **CATH: Direct 3D/SSE Rules** | 17.3% | 13.7% | **15.5%** |
| **CATH: TED T-level (`.T10`)** | 12.5% | 10.6% | **11.6%** |
| **CATH: Gene3D Overlap** | 3.2% | 3.2% | **3.2%** |
| **Difficult Cases (Confidence < 0.50)** | 0.1% | 0.0% | **0.1%** (7 / 10,000) |
| **Low pLDDT (< 60)** | 0.8% | 1.0% | **0.9%** |

Full benchmark reports and datasets are available in the [`examples/`](examples/) directory.

---

## Running Tests

```bash
uv run --extra dev pytest
```

---

## Repository & Citation

- **Repository**: [https://github.com/KarelBerka/TED-namer](https://github.com/KarelBerka/TED-namer)
- **Author**: Karel Berka
- **License**: MIT
