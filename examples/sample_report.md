# TED Structural Domain Cluster Naming Report

| Cluster ID | Members | Consensus Domain Name | Short Label | Tier | CATH Label | Confidence |
| :--- | :---: | :--- | :--- | :--- | :---: | :---: |
| `p53_DBD_cluster` | 3 | **P53 DNA-binding domain** | `P53_DNA-bd` | `specific_domain` | `2.60.40.720` | 0.99 |
| `p63_p73_SAM_cluster` | 2 | **Sterile alpha motif domain** | `Sterile_alpha_motif_dom` | `specific_domain` | `1.10.150.50` | 0.99 |
| `ErbB_Receptor_L_cluster` | 2 | **Receptor L-domain** | `Receptor_L-dom` | `superfamily` | `3.80.20.20` | 0.89 |

## Cluster `p53_DBD_cluster`: P53 DNA-binding domain (`P53_DNA-bd`)

- **Confidence Score**: `0.990` (`rule_based`)
- **Specificity Tier**: `specific_domain`
- **Dominant CATH Classification**: `2.60.40.720` — Superfamily: *p53/RUNT-type transcription factor, DNA-binding domain superfamily*, Topology: *Immunoglobulin-like*
- **Consensus Relative Position**: `central`
- **Rationale**: Selected 'P53 DNA-binding domain' (specific_domain) supported by 3/3 cluster members (100%) across [cdd, interpro, pfam] (IPR011615, PF00870, cd08367) with mean residue IoU=0.98 and TED coverage=0.99. Structural fold classification: CATH 2.60.40.720 (Immunoglobulin-like).

### Top Ranked Candidate Annotations

| Rank | Candidate Name | Tier | Freq | Mean IoU | TED Cov | Databases | Accessions | Score |
| :---: | :--- | :--- | :---: | :---: | :---: | :--- | :--- | :---: |
| 1 | P53 DNA-binding domain | `specific_domain` | 100% (3/3) | 0.98 | 0.99 | cdd, interpro, pfam | `IPR011615, PF00870, cd08367` | **1.000** |
| 2 | p53/RUNT-type transcription factor DNA-binding domain superfamily | `superfamily` | 100% (3/3) | 0.89 | 1.00 | cath_h, cathgene3d, interpro | `CATH:2.60.40.720, G3DSA:2.60.40.720, IPR012346` | **0.879** |
| 3 | p53-like transcription factor DNA-binding domain superfamily | `superfamily` | 100% (3/3) | 0.92 | 0.98 | interpro, ssf | `IPR008967, SSF49417` | **0.871** |
| 4 | DNA-binding domain | `specific_domain` | 33% (1/3) | 0.99 | 0.99 | uniprot_feature | `UP_FEAT:O15350:131-310` | **0.722** |
| 5 | Interaction with AXIN1 domain | `specific_domain` | 33% (1/3) | 0.93 | 0.96 | uniprot_feature | `UP_FEAT:P04637:116-292` | **0.695** |
| 6 | Required for interaction with ZNF385A domain | `specific_domain` | 33% (1/3) | 0.89 | 1.00 | uniprot_feature | `UP_FEAT:P04637:100-300` | **0.682** |

### Cluster Member TED Domains

| TED Domain ID | UniProt | Protein Name | Chopping | Position | pLDDT | CATH Label | Top Overlapping Hit (IoU) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| [AF-P04637-F1-model_v4_TED01](https://ted.cathdb.info/uniprot/P04637) | `P04637` | Cellular tumor antigen p53 | `109-287` | `central` | 95.1 | `2.60.40.720` | P53 DNA-binding domain (`cd08367`, IoU=0.99) |
| [AF-Q9H3D4-F1-model_v4_TED01](https://ted.cathdb.info/uniprot/Q9H3D4) | `Q9H3D4` | Tumor protein 63 | `177-361` | `central` | 93.9 | `2.60.40.720` | P53 DNA-binding domain (`cd08367`, IoU=0.98) |
| [AF-O15350-F1-model_v4_TED01](https://ted.cathdb.info/uniprot/O15350) | `O15350` | Tumor protein p73 | `129-310` | `N-terminal` | 94.4 | `2.60.40.720` | DNA-binding (`UP_FEAT:O15350:131-310`, IoU=0.99) |

## Cluster `p63_p73_SAM_cluster`: Sterile alpha motif domain (`Sterile_alpha_motif_dom`)

- **Confidence Score**: `0.990` (`rule_based`)
- **Specificity Tier**: `specific_domain`
- **Dominant CATH Classification**: `1.10.150.50` — Superfamily: *Transcription Factor, Ets-1*, Topology: *DNA polymerase; domain 1*
- **Consensus Relative Position**: `C-terminal`
- **Rationale**: Selected 'Sterile alpha motif domain' (specific_domain) supported by 2/2 cluster members (100%) across [interpro, pfam, smart] (IPR001660, PF07647, SM00454) with mean residue IoU=0.88 and TED coverage=0.98. Structural fold classification: CATH 1.10.150.50 (DNA polymerase; domain 1).

### Top Ranked Candidate Annotations

| Rank | Candidate Name | Tier | Freq | Mean IoU | TED Cov | Databases | Accessions | Score |
| :---: | :--- | :--- | :---: | :---: | :---: | :--- | :--- | :---: |
| 1 | Sterile alpha motif domain | `specific_domain` | 100% (2/2) | 0.88 | 0.98 | interpro, pfam, smart | `IPR001660, PF07647, SM00454` | **1.000** |
| 2 | SAM domain | `specific_domain` | 100% (2/2) | 0.88 | 0.98 | uniprot_feature | `UP_FEAT:O15350:485-551, UP_FEAT:Q9H3D4:541-607` | **0.986** |
| 3 | Transcription Factor Ets-1 | `superfamily` | 100% (2/2) | 0.91 | 1.00 | cath_h, cathgene3d, interpro, ssf | `CATH:1.10.150.50, G3DSA:1.10.150.50, IPR013761, SSF47769` | **0.903** |
| 4 | SAM domain of tumor-p63 proteins | `specific_domain` | 50% (1/2) | 0.91 | 1.00 | cdd, interpro | `IPR037611, cd09572` | **0.783** |
| 5 | SAM domain of tumor-p73 proteins | `specific_domain` | 50% (1/2) | 0.90 | 0.97 | cdd, interpro | `IPR037612, cd09571` | **0.773** |
| 6 | p53 C-terminal domain | `protein_family_context` | 100% (2/2) | 0.60 | 0.80 | uniprot_family | `O15350, Q9H3D4` | **0.548** |

### Cluster Member TED Domains

| TED Domain ID | UniProt | Protein Name | Chopping | Position | pLDDT | CATH Label | Top Overlapping Hit (IoU) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| [AF-Q9H3D4-F1-model_v4_TED02](https://ted.cathdb.info/uniprot/Q9H3D4) | `Q9H3D4` | Tumor protein 63 | `548-606` | `C-terminal` | 90.7 | `1.10.150.50` | SAM domain of tumor-p63 proteins (`cd09572`, IoU=0.91) |
| [AF-O15350-F1-model_v4_TED02](https://ted.cathdb.info/uniprot/O15350) | `O15350` | Tumor protein p73 | `492-553` | `C-terminal` | 85.9 | `1.10.150.50` | SAM/Pointed domain (`SSF47769`, IoU=0.93) |

## Cluster `ErbB_Receptor_L_cluster`: Receptor L-domain (`Receptor_L-dom`)

- **Confidence Score**: `0.894` (`rule_based`)
- **Specificity Tier**: `superfamily`
- **Dominant CATH Classification**: `3.80.20.20` — Superfamily: *Receptor L-domain*, Topology: *24 nucleotide stem-loop, u2 snrnp hairpin iv. U2 a'; Chain A*
- **Consensus Relative Position**: `N-terminal`
- **Rationale**: Selected 'Receptor L-domain' (superfamily) supported by 2/2 cluster members (100%) across [cath_h, cathgene3d, interpro] (CATH:3.80.20.20, G3DSA:3.80.20.20, IPR036941) with mean residue IoU=0.94 and TED coverage=0.99.

### Top Ranked Candidate Annotations

| Rank | Candidate Name | Tier | Freq | Mean IoU | TED Cov | Databases | Accessions | Score |
| :---: | :--- | :--- | :---: | :---: | :---: | :--- | :--- | :---: |
| 1 | Receptor L-domain | `superfamily` | 100% (2/2) | 0.94 | 0.99 | cath_h, cathgene3d, interpro | `CATH:3.80.20.20, G3DSA:3.80.20.20, IPR036941` | **0.894** |
| 2 | L domain-like | `superfamily` | 100% (2/2) | 0.90 | 0.93 | ssf | `SSF52058` | **0.839** |
| 3 | Receptor L-domain | `specific_domain` | 100% (2/2) | 0.60 | 0.61 | interpro, pfam | `IPR000494, PF01030` | **0.829** |
| 4 | Approximate domain | `specific_domain` | 50% (1/2) | 0.49 | 0.75 | uniprot_feature | `UP_FEAT:P00533:75-300` | **0.556** |
| 5 | protein kinase N-terminal domain | `protein_family_context` | 100% (2/2) | 0.60 | 0.80 | uniprot_family | `P00533, P04626` | **0.548** |
| 6 | Tyrosine protein kinase EGF/ERB/XmrK receptor N-terminal domain | `protein_family_context` | 100% (2/2) | 0.15 | 1.00 | interpro | `IPR016245` | **0.446** |

### Cluster Member TED Domains

| TED Domain ID | UniProt | Protein Name | Chopping | Position | pLDDT | CATH Label | Top Overlapping Hit (IoU) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| [AF-P00533-F1-model_v4_TED01](https://ted.cathdb.info/uniprot/P00533) | `P00533` | Epidermal growth factor receptor | `31-207` | `N-terminal` | 92.0 | `3.80.20.20` | L domain-like (`SSF52058`, IoU=0.99) |
| [AF-P04626-F1-model_v4_TED01](https://ted.cathdb.info/uniprot/P04626) | `P04626` | Receptor tyrosine-protein kinase erbB-2 | `20-123_133-215` | `N-terminal` | 91.1 | `3.80.20.20` | Receptor L-domain (`G3DSA:3.80.20.20`, IoU=0.93) |
