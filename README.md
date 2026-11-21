# `ted-namer`: Structural Domain & Cluster Namer for TED

[![CI](https://github.com/KarelBerka/TED-namer/actions/workflows/ci.yml/badge.svg)](https://github.com/KarelBerka/TED-namer/actions/workflows/ci.yml)
[![JOSS Paper Draft](https://github.com/KarelBerka/TED-namer/actions/workflows/joss-paper.yml/badge.svg)](https://github.com/KarelBerka/TED-namer/actions/workflows/joss-paper.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Code style: uv](https://img.shields.io/badge/managed%20by-uv-orange.svg)](https://github.com/astral-sh/uv)

Generates standardized, biologically meaningful human-readable names and calibrated **CATH v4.4** structural classifications for domains and domain clusters from **TED (The Encyclopedia of Domains, [ted.cathdb.info](https://ted.cathdb.info))**.

Detailed methodology and literature validation criteria are documented in [**`docs/NAMING_PROCESS.md`**](docs/NAMING_PROCESS.md) and [**`docs/ACCURACY_VALIDATION.md`**](docs/ACCURACY_VALIDATION.md).

---

## Key Features

- **Boundary-Aware Overlap Engine**: Computes exact residue-level Intersection-over-Union (IoU), TED coverage, and annotation coverage across both continuous (`109-287`) and discontinuous (`15-65_140-210`) TED segments against **InterPro, Pfam, CDD, SMART, Profile, HAMAP, NCBIfam, PANTHER, SUPERFAMILY, CATH-Gene3D, and UniProt features**.
- **Strict `annotation_type` & Multi-Lobe Subdomain Resolution**: Honors InterPro entry types (distinguishing true `domain`/`repeat` matches and strict 1:1 `Pfam family` domains from whole-protein `HAMAP`/`NCBIfam` families and `PROSITE` sites) and automatically resolves split structural lobes (e.g., *Protein kinase N-terminal subdomain*).
- **CATH v4.4 & `?`-Wildcard Classification Engine**: Queries the **CATH v4.4 REST API** (with `v4_3_0` fallback) covering Classes 1–4 and **Class 6 (`Special`)**, prioritizes curated integrated InterPro superfamily names and curated unintegrated H-overrides over legacy PDB titles, and avoids colliding with TED's `T`-level nomenclature by using **`?` wildcards** for unconfirmed levels (`2.60.40.?`, `3.?.?.?`) while exporting geometric/profile hints in `cath_fold_hint` and flagging unmapped folds via `novel_fold_candidate`.
- **Naming Hygiene & Structural Quality Separation**: Sanitizes sentence-like descriptions, strips `[gene]` tags and duplicate `domain` suffixes, enriches bare positional labels (*"N-terminal domain"* $\to$ *"Argonaute N-terminal domain"*), properly capitalizes consensus titles while preserving biochemical prefixes (`p53`, `tRNA`, `von Willebrand`, `beta-`), and reports TED structural quality (`structure_quality_flags`: `ok`, `low_plddt`, `low_consensus`, `non_globular`) independently of naming confidence.
- **High-Performance Async Pipeline**: WAL-backed SQLite API caching, batch UniProt `/accessions` prefetching (100 proteins/request), and concurrent cluster processing.

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

Retrieve all AlphaFold-predicted TED domains, residue boundaries, structural quality flags, and consensus CATH classifications for a UniProt accession (e.g., Human p53, `P04637`):

```bash
uv run ted-namer inspect P04637 --no-llm
```

Output:
```text
╭──────────────────────── AF-P04637-F1-model_v4_TED01 ────────────────────────╮
│ Proposed Name: P53 DNA-binding domain (P53_DNA-bd, conf=0.99,               │
│ tier=specific_domain)                                                       │
│ Parent Protein: Cellular tumor antigen p53 (P04637, Homo sapiens)           │
│ Chopping: 109-287 (179 res, central, pLDDT=95.1)                            │
│ CATH Classification: 2.60.40.720 (p53/RUNT-type transcription factor,       │
│ DNA-binding domain superfamily / Immunoglobulin-like)                       │
│ Novel Fold Candidate: no | Structure Quality: ok                            │
│ Selected 'P53 DNA-binding domain' (specific_domain) supported by 1/1        │
│ cluster members (100%) across [cdd, interpro, pfam] (IPR011615, PF00870,    │
│ cd08367) with mean residue IoU=0.99 and TED coverage=1.00.                  │
╰─────────────────────────────────────────────────────────────────────────────╯
```

### 2. Name a Single TED Domain or Homologous Cluster

Name a cluster of orthologous/paralogous TED domains (e.g., DNA-binding domains of p53, p63, and p73):

```bash
uv run ted-namer name \
  --ted-ids AF-P04637-F1-model_v4_TED01,AF-Q9H3D4-F1-model_v4_TED01,AF-O15350-F1-model_v4_TED01 \
  --cluster-id p53_family_DBD \
  --no-llm \
  -o results.json -o results.tsv -o report.md
```

### 3. Batch Naming from a Table

Provide a TSV or CSV with columns `cluster_id` and `ted_id`:

```bash
uv run ted-namer name \
  --input examples/sample_clusters.tsv \
  --concurrency 25 \
  --no-llm \
  -o batch_results.tsv -o batch_report.md
```

---

## Python API Usage

```python
import asyncio
from ted_namer.pipeline import TedNamingPipeline

async def main():
    pipeline = TedNamingPipeline(use_cache=True, concurrency=20, use_llm=False)

    clusters = {
        "p53_DBD": ["AF-P04637-F1-model_v4_TED01"],
        "p63_SAM": ["AF-Q9H3D4-F1-model_v4_TED02"],
    }
    results = await pipeline.name_clusters(clusters)

    for res in results:
        print(f"Cluster: {res.cluster_id}")
        print(f"  Consensus Name:  {res.consensus_name}")
        print(f"  Short Name:      {res.short_name}")
        print(f"  Tier:            {res.specificity_tier.value}")
        print(f"  Confidence:      {res.confidence_score:.3f}")
        print(f"  CATH / Hint:     {res.dominant_cath_label} (hint={res.cath_fold_hint}, novel={res.novel_fold_candidate})")
        print(f"  Quality Flags:   {', '.join(res.structure_quality_flags)}")
        for m in res.members:
            print(f"  Member: {m.domain.ted_id} -> {m.domain.inferred_cath_label}")

if __name__ == "__main__":
    asyncio.run(main())
```

---

## The 4-Tier Specificity Hierarchy

| Tier | Name | Typical Sources | Example |
| :---: | :--- | :--- | :--- |
| **1** | **Specific Domain / Subdomain** | CDD, Pfam (`domain`/`repeat`/1:1 `family`), SMART, InterPro domain, multi-lobe subdomains | *p53 DNA-binding domain*, *Protein kinase N-terminal subdomain* |
| **2** | **Homologous Superfamily** | CATH H-level (curated IPR / unintegrated H-overrides), SUPERFAMILY | *Protein kinase-like N-terminal lobe superfamily* |
| **3** | **Protein Family Context** | `HAMAP`, `NCBIfam family`, `PANTHER`, UniProt `SIMILARITY` family + position | *BioU C-terminal domain*, *SEDS central domain* |
| **4** | **Topology / Fold** | Confirmed-`T` CATH v4.4 Topology fallback (when uncharacterized) | *Immunoglobulin-like fold* |

---

## CATH v4.4 Classification, `?` Wildcards & Novel Fold Candidates

To avoid colliding with CATH/TED's **Topology (`T`) level**, unconfirmed hierarchy levels use explicit **`?` wildcards**, while underlying `Profile-to-CATH` and 3D $\text{C}\alpha$ geometric estimates (Class $C$ from SSEs, helix crossing angles, $\beta$-strand contact graphs, SVD $\beta$-sheet layers) are stored in `cath_fold_hint` alongside `novel_fold_candidate`:
- **Confirmed Direct Depth**: Plain integers only where empirical accuracy exceeds **90%** (`2.60.40.720`, `2.60.40.?`, `4.10.?.?`, `3.?.?.?`).
- **CATH v4.4 & Class 6**: Queries `v4_4_0` API (`v4_3_0` fallback) and resolves Classes **1–4** and **Class 6 (`Special`)**.

### Examples of Domains Resolved Across Literature Edge Cases (`N = 100,000`)

| Domain ID | TED Raw CATH | Assigned CATH (`fold_hint`, `novel`) | Resolution Mechanism | Assigned Name & Tier |
| :--- | :---: | :---: | :--- | :--- |
| `AF-Q8IWB3-F1-model_v4_TED01` | `1.10.510.10` | **`1.10.510.10`** *(`-`, `false`)* | Legacy `Transferase(Phosphotransferase) domain 1` supplanted by Pfam | **Protein kinase domain** *(Tier 1, conf: 0.722)* |
| `AF-Q6DBX4-F1-model_v4_TED01` | `1.10.510.10` | **`1.10.510.10`** *(`-`, `false`)* | Uninflated `cath_h` + curated `IPR011009` beats legacy PDB header | **Protein kinase-like domain superfamily** *(Tier 2, conf: 0.729)* |
| `AF-A4F4K3-F1-model_v4_TED02` | `-` | **`1.10.630.?`** *(`1.10.630.10`, `false`)* | `g3dsa_overlap` (`known_depth=3`) + multi-lobe subdomain | **Cytochrome P450 central subdomain** *(Tier 1, conf: 0.772)* |
| `AF-Q8U3D2-F1-model_v4_TED01` | `-` | **`3.90.70.?`** *(`3.90.70.180`, `false`)* | `Profile-to-CATH` bridge + bare `"N-terminal domain"` enrichment | **Argonaute N-terminal domain** *(Tier 1, conf: 0.820)* |
| `AF-A0A1G9Z131-F1-model_v4_TED02` | `-` | **`3.?.?.?`** *(`3.40.50.720`, `false`)* | `HAMAP` whole-protein `family` routed to Tier 3 + `[bioU]` stripped | **BioU C-terminal domain** *(Tier 3, conf: 0.451)* |
| `AF-O05442-F1-model_v4_TED01` | `-` | **`1.?.?.?`** *(`1.25.40.10`, `true`)* | Bare `NTD domain` enriched + unassigned `rule_sse` fold candidate | **Outer membrane channel protein CpnT N-terminal domain** *(Tier 3, conf: 0.511)* |

*Full case studies and geometric derivations are detailed in [docs/NAMING_PROCESS.md](docs/NAMING_PROCESS.md).*

---

## Blind Cross-Validation & Ground-Truth Accuracy (100,000 TED Domains)

Detailed ground-truth validation across **100,000 unique TED domains** (100 batches of 1,000 domains / 10 deciles of 10,000 domains) is documented in [**`docs/ACCURACY_VALIDATION.md`**](docs/ACCURACY_VALIDATION.md) and reproducible via [`scripts/benchmark_100000_random.py`](scripts/benchmark_100000_random.py) and [`scripts/reevaluate_100k.py`](scripts/reevaluate_100k.py):

### 1. Blind Hold-Out CATH Classification Accuracy ($N = 72,130$ Official CATH Domains)

| Hold-Out Scenario | Class ($C$) | Architecture ($C.A$) | Topology ($C.A.T$) | Superfamily ($C.A.T.H$) |
| :--- | :---: | :---: | :---: | :---: |
| **Full Fallback Cascade** *(Mask TED `cath_label`)* | **98.2%** | **95.2%** | **93.9%** | **91.1%** |
| **`Profile-to-CATH` Bridge (5-Fold CV)** *(Mask TED CATH & `G3DSA`; 87.9% cov)* | **99.3%** | **98.6%** | **97.9%** | **92.3%** |
| **Pure 3D/SSE Rules Only** *(Zero sequence/profile homology)* | **85.6%** *(96.7% pure $\alpha/\beta$)* | **41.5%** *(in `cath_fold_hint`)* | **25.4%** *(in `cath_fold_hint`)* | — |

- **Precision of Confirmed Direct (Non-`?`) Levels ($N = 72,130$)**: **98.2%** for $C$, **98.6%** for $C.A$, **98.1%** for $C.A.T$, and **99.8%** for $C.A.T.H$.
- **Leave-One-Database-Out Naming Concordance ($N = 36,843$)**: **91.2%** exact InterPro/normalized name concordance (`33,611 / 36,843`) and **96.4%** semantic keyword concordance (`35,534 / 36,843`) when holding out `Pfam`.

---

## Large-Scale Benchmark & Naming Hygiene Audit (100,000 TED Domains)

Evaluated across **100,000 unique TED domains** (10 consecutive deciles of 10,000 domains across all kingdoms of life):

| Metric | Range Across 10 Deciles (`10x 10k`) | **Total Across 100,000 TED Domains** |
| :--- | :---: | :---: |
| **Domain Count** | 10,000 per decile | **100,000** |
| **Mean Calibrated Confidence Score** | 0.791 – 0.818 | **0.798** *(Tier 1: 0.889, Tier 2: 0.769, Tier 3: 0.422)* |
| **Tier 1: Specific Domain & Subdomain** | 61.4% – 63.3% | **62.43%** (62,425 / 100,000; incl. 5,934 subdomains) |
| **Tier 2: Superfamily (Curated IPR & H-overrides)** | 23.2% – 25.4% | **24.54%** (24,537 / 100,000) |
| **Tier 1 + Tier 2 (Domain or Superfamily)** | **85.8% – 88.3%** | **86.96%** (86,962 / 100,000) |
| **Tier 3: Protein Family Context** | 10.6% – 14.2% | **12.86%** (12,858 / 100,000) |
| **Tier 4 / 5: Confirmed Fold / Uncharacterized** | 0.00% – 1.15% | **0.18%** (70 fold + 110 uncharacterized) |
| **CATH: Official (4-level, Classes 1–4 & 6)** | 69.8% – 73.2% | **72.13%** (72,130 / 100,000) |
| **CATH: Partial TED (`C.A.T.?`) & Homology Bridge** | 15.6% – 19.7% | **17.13%** *(9.7% TED `T` + 4.3% Profile + 3.1% `G3DSA`)* |
| **CATH: Novel Fold Candidates (`C.?.?.?` / `4.10.?.?`)** | 9.6% – 11.4% | **10.74%** (10,740 / 100,000 flagged `novel_fold_candidate`) |
| **Hygiene Artifacts Eliminated (`bare_pos` / `[gene]` / `T`-prefix)** | — | **0** remaining (strict PDB noise **0.036%**) |
| **Structure Quality (`ok` / `low_plddt` / `non_globular`)** | — | **96.2%** / **2.1%** / **1.9%** |

Full benchmark outputs ([`examples/random_100000_results.tsv.gz`](examples/random_100000_results.tsv.gz), [`examples/benchmark_100000_summary.json`](examples/benchmark_100000_summary.json)) and worked examples ([`examples/sample_clusters.tsv`](examples/sample_clusters.tsv), [`examples/sample_results.tsv`](examples/sample_results.tsv), [`examples/sample_report.md`](examples/sample_report.md)) are included in [`examples/`](examples/).

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
  version   = {0.2.0}
}
```

- **Repository**: [https://github.com/KarelBerka/TED-namer](https://github.com/KarelBerka/TED-namer)
- **Author**: Karel Berka ([ORCID: 0000-0001-9472-2589](https://orcid.org/0000-0001-9472-2589)), Palacký University Olomouc
- **License**: [MIT](LICENSE)
