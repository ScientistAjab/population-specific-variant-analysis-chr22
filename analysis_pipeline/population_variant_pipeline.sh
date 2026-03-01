#!/usr/bin/env bash
set -euo pipefail

echo "Starting population-specific analysis (INFO-based)..."

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${PROJECT_DIR}/data"
RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/analysis_pipeline/logs"

mkdir -p "${RESULTS_DIR}/comparison" \
         "${RESULTS_DIR}/comparison/loose_threshold" \
         "${LOG_DIR}"

RAW_VCF="${DATA_DIR}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
GENOME_BUILD="GRCh37.75"

############################################
# FUNCTION: Clean + Annotate + Summarize
############################################

process_population () {

    POP=$1
    INPUT_VCF=$2
    OUTDIR=$3

    echo "Processing ${POP}..."

    # Remove structural variants
    bcftools view -v snps,indels "$INPUT_VCF" -Oz \
    -o "${OUTDIR}/${POP}_clean.vcf.gz"
    tabix -p vcf "${OUTDIR}/${POP}_clean.vcf.gz"

    # SNP / INDEL counts
    bcftools view -v snps "${OUTDIR}/${POP}_clean.vcf.gz" | grep -v "^#" | wc -l \
    > "${OUTDIR}/${POP}_snp_count.txt"

    bcftools view -v indels "${OUTDIR}/${POP}_clean.vcf.gz" | grep -v "^#" | wc -l \
    > "${OUTDIR}/${POP}_indel_count.txt"

    # Annotation
    snpEff -Xmx8g ${GENOME_BUILD} \
    -stats "${OUTDIR}/${POP}_summary.html" \
    "${OUTDIR}/${POP}_clean.vcf.gz" \
    | bgzip > "${OUTDIR}/${POP}_annotated.vcf.gz"

    tabix -p vcf "${OUTDIR}/${POP}_annotated.vcf.gz"

    # Impact counts (proper ANN parsing)
    bcftools query -f '%INFO/ANN\n' \
    "${OUTDIR}/${POP}_annotated.vcf.gz" \
    | tr ',' '\n' \
    | awk -F'|' '{print $3}' \
    | sort | uniq -c | sort -nr \
    > "${OUTDIR}/${POP}_impact_counts.txt"

    # Variant effect type counts (missense etc.)
    bcftools query -f '%INFO/ANN\n' \
    "${OUTDIR}/${POP}_annotated.vcf.gz" \
    | tr ',' '\n' \
    | awk -F'|' '{print $2}' \
    | sort | uniq -c | sort -nr \
    > "${OUTDIR}/${POP}_effect_type_counts.txt"

    # HIGH impact gene names
    bcftools query -f '%INFO/ANN\n' \
    "${OUTDIR}/${POP}_annotated.vcf.gz" \
    | tr ',' '\n' \
    | awk -F'|' '$3=="HIGH" {print $4}' \
    | sort -u \
    > "${OUTDIR}/${POP}_HIGH_impact_genes.txt"

}

############################################
# 1️⃣ STRICT THRESHOLD
############################################

echo "Running STRICT threshold filtering..."

bcftools view -i 'INFO/EUR_AF>0.3 && INFO/SAS_AF<0.05' \
"${RAW_VCF}" -Oz \
-o "${RESULTS_DIR}/comparison/eur_strict.vcf.gz"

bcftools view -i 'INFO/EUR_AF<0.01 && INFO/SAS_AF>0.05' \
"${RAW_VCF}" -Oz \
-o "${RESULTS_DIR}/comparison/sas_strict.vcf.gz"

tabix -p vcf "${RESULTS_DIR}/comparison/eur_strict.vcf.gz"
tabix -p vcf "${RESULTS_DIR}/comparison/sas_strict.vcf.gz"

process_population "eur_strict" \
"${RESULTS_DIR}/comparison/eur_strict.vcf.gz" \
"${RESULTS_DIR}/comparison"

process_population "sas_strict" \
"${RESULTS_DIR}/comparison/sas_strict.vcf.gz" \
"${RESULTS_DIR}/comparison"

############################################
# 2️⃣ LOOSE THRESHOLD
############################################

echo "Running LOOSE threshold filtering..."

bcftools view -i 'INFO/EUR_AF>0.2 && INFO/SAS_AF<0.1' \
"${RAW_VCF}" -Oz \
-o "${RESULTS_DIR}/comparison/loose_threshold/eur_loose.vcf.gz"

bcftools view -i 'INFO/EUR_AF<0.05 && INFO/SAS_AF>0.1' \
"${RAW_VCF}" -Oz \
-o "${RESULTS_DIR}/comparison/loose_threshold/sas_loose.vcf.gz"

tabix -p vcf "${RESULTS_DIR}/comparison/loose_threshold/eur_loose.vcf.gz"
tabix -p vcf "${RESULTS_DIR}/comparison/loose_threshold/sas_loose.vcf.gz"

process_population "eur_loose" \
"${RESULTS_DIR}/comparison/loose_threshold/eur_loose.vcf.gz" \
"${RESULTS_DIR}/comparison/loose_threshold"

process_population "sas_loose" \
"${RESULTS_DIR}/comparison/loose_threshold/sas_loose.vcf.gz" \
"${RESULTS_DIR}/comparison/loose_threshold"

echo "Population-specific analysis complete."

echo "Generating comparative high-impact gene analysis..."

STRICT_DIR="${RESULTS_DIR}/comparison"

EUR_GENES="${STRICT_DIR}/eur_strict_HIGH_impact_genes.txt"
SAS_GENES="${STRICT_DIR}/sas_strict_HIGH_impact_genes.txt"

# Ensure sorted
sort -u "$EUR_GENES" > "${STRICT_DIR}/eur_sorted.txt"
sort -u "$SAS_GENES" > "${STRICT_DIR}/sas_sorted.txt"

# Genes common in both
comm -12 "${STRICT_DIR}/eur_sorted.txt" "${STRICT_DIR}/sas_sorted.txt" \
> "${STRICT_DIR}/HIGH_genes_shared.txt"

# Unique to EUR (common in EUR, rare in SAS)
comm -23 "${STRICT_DIR}/eur_sorted.txt" "${STRICT_DIR}/sas_sorted.txt" \
> "${STRICT_DIR}/HIGH_genes_unique_to_EUR.txt"

# Unique to SAS (common in SAS, rare in EUR)
comm -13 "${STRICT_DIR}/eur_sorted.txt" "${STRICT_DIR}/sas_sorted.txt" \
> "${STRICT_DIR}/HIGH_genes_unique_to_SAS.txt"

# Count summary
echo "HIGH impact gene summary (STRICT):" \
> "${STRICT_DIR}/HIGH_gene_summary.txt"

echo "EUR-specific HIGH genes:" \
>> "${STRICT_DIR}/HIGH_gene_summary.txt"
wc -l "${STRICT_DIR}/HIGH_genes_unique_to_EUR.txt" \
>> "${STRICT_DIR}/HIGH_gene_summary.txt"

echo "SAS-specific HIGH genes:" \
>> "${STRICT_DIR}/HIGH_gene_summary.txt"
wc -l "${STRICT_DIR}/HIGH_genes_unique_to_SAS.txt" \
>> "${STRICT_DIR}/HIGH_gene_summary.txt"

echo "Shared HIGH genes:" \
>> "${STRICT_DIR}/HIGH_gene_summary.txt"
wc -l "${STRICT_DIR}/HIGH_genes_shared.txt" \
>> "${STRICT_DIR}/HIGH_gene_summary.txt"

############################################
# 🔬 LOOSE Threshold Comparative Analysis
############################################

echo "Generating comparative high-impact gene analysis (LOOSE)..."

LOOSE_DIR="${RESULTS_DIR}/comparison/loose_threshold"

EUR_GENES="${LOOSE_DIR}/eur_loose_HIGH_impact_genes.txt"
SAS_GENES="${LOOSE_DIR}/sas_loose_HIGH_impact_genes.txt"

sort -u "$EUR_GENES" > "${LOOSE_DIR}/eur_loose_sorted.txt"
sort -u "$SAS_GENES" > "${LOOSE_DIR}/sas_loose_sorted.txt"

# Shared
comm -12 "${LOOSE_DIR}/eur_loose_sorted.txt" \
         "${LOOSE_DIR}/sas_loose_sorted.txt" \
> "${LOOSE_DIR}/HIGH_genes_shared_loose.txt"

# Unique to EUR
comm -23 "${LOOSE_DIR}/eur_loose_sorted.txt" \
         "${LOOSE_DIR}/sas_loose_sorted.txt" \
> "${LOOSE_DIR}/HIGH_genes_unique_to_EUR_loose.txt"

# Unique to SAS
comm -13 "${LOOSE_DIR}/eur_loose_sorted.txt" \
         "${LOOSE_DIR}/sas_loose_sorted.txt" \
> "${LOOSE_DIR}/HIGH_genes_unique_to_SAS_loose.txt"

# Summary
echo "HIGH impact gene summary (LOOSE):" \
> "${LOOSE_DIR}/HIGH_gene_summary_loose.txt"

echo "EUR-specific HIGH genes:" \
>> "${LOOSE_DIR}/HIGH_gene_summary_loose.txt"
wc -l "${LOOSE_DIR}/HIGH_genes_unique_to_EUR_loose.txt" \
>> "${LOOSE_DIR}/HIGH_gene_summary_loose.txt"

echo "SAS-specific HIGH genes:" \
>> "${LOOSE_DIR}/HIGH_gene_summary_loose.txt"
wc -l "${LOOSE_DIR}/HIGH_genes_unique_to_SAS_loose.txt" \
>> "${LOOSE_DIR}/HIGH_gene_summary_loose.txt"

echo "Shared HIGH genes:" \
>> "${LOOSE_DIR}/HIGH_gene_summary_loose.txt"
wc -l "${LOOSE_DIR}/HIGH_genes_shared_loose.txt" \
>> "${LOOSE_DIR}/HIGH_gene_summary_loose.txt"

























