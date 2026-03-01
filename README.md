# 🧬 Population-Specific Variant Analysis of Chromosome 22  
### 1000 Genomes Project | EUR vs SAS | Functional Impact Profiling

---

## 🚀 Project Overview

This project performs a **population-specific variant analysis** on **Chromosome 22** using data from the 1000 Genomes Project.

Two populations were analyzed:

- 🇪🇺 European (EUR)  
- 🇸🇦 South Asian (SAS)  

Two filtering strategies were applied:

- 🔒 Strict threshold  
- 🔓 Loose threshold  

Variants were functionally annotated using **SnpEff**, and population-specific differences were evaluated across:

- SNP and INDEL burden  
- Functional impact categories (HIGH / MODERATE / LOW / MODIFIER)  
- Effect type distribution  
- Missense and loss-of-function variants  
- High-impact gene identification  

---

# 📊 Key Results

## 🔹 Variant Burden

| Population | Threshold | SNPs | INDELs |
|------------|------------|------|--------|
| EUR | Loose | 1511 | 692 |
| SAS | Loose | 4699 | 1125 |
| SAS | Strict | 4125 | 1024 |

➡ The South Asian population consistently shows **2–3× higher variant counts** across thresholds.

---

## 🔹 Functional Impact Distribution

### EUR (Loose)
- 13,898 MODIFIER  
- 70 LOW  
- 54 MODERATE  
- 6 HIGH  

### SAS (Loose)
- 30,731 MODIFIER  
- 187 LOW  
- 112 MODERATE  
- 13 HIGH  

🔥 SAS exhibits:
- ~2× more MODERATE variants  
- >2× more HIGH-impact variants  
- Greater overall functional diversity  

---

## 🔹 Effect Type Distribution

The majority of variants in both populations were:

- intron_variant  
- upstream_gene_variant  
- downstream_gene_variant  
- intergenic_region  

Example (Loose threshold):

- EUR intron variants: 8,096  
- SAS intron variants: 17,534  

➡ Most population-specific variants reside in **non-coding regulatory regions**.

---

## 🔹 Coding Variant Insights

### Missense Variants
- EUR (Loose): 50  
- SAS (Loose): 94  

### Severe Functional Variants (SAS Loose)
- 1 stop_gained  
- 1 frameshift_variant  
- 2 start_lost  

➡ Indicates broader protein-altering variation in SAS.

---

# 🧬 High-Impact Genes

## Shared Across Populations
- PI4KAP1  
- TMPRSS6  

## Unique to South Asian Population
- ADORA2A-AS1  
- CRYBB2P1  
- POM121L9P  
- PRR5-ARHGAP8  
- CTA-299D3.8  
- VPREB1  

❗ No HIGH-impact genes were unique to EUR.

---

# 🧠 Biological Relevance of Identified Genes

### TMPRSS6
- Regulates iron homeostasis  
- Controls hepcidin expression  
- Associated with iron-refractory anemia  

### VPREB1
- Essential component of pre-B cell receptor  
- Critical for B-cell development  
- Immune-related gene  

### ADORA2A-AS1
- Antisense RNA regulating ADORA2A  
- Involved in inflammatory and immune signaling  

### PRR5-ARHGAP8
- Involved in cytoskeleton and Rho GTPase signaling  
- Influences cell migration  

➡ Findings suggest ancestry-specific divergence particularly in **immune and regulatory pathways**.

---

# 🛠️ Pipeline Workflow

```bash
1. Extract chr22 from 1000G VCF
2. Filter population-specific variants
3. Apply strict & loose thresholds
4. Annotate using SnpEff
5. Extract:
   - SNP / INDEL counts
   - Impact category counts
   - Effect type distribution
   - HIGH-impact genes
6. Perform EUR vs SAS comparison


population-specific-variant-analysis-chr22/
│
├── analysis_pipeline/
├── data
├── Figures
├── result files/
│   ├── strict_threshold/
│   └── loose_threshold/   
├── results and conclusions/
└── README.md
