#!/bin/bash

# =============================================================================
# Main script: generate aligned count files in batch and integrate tumor ASE data
# =============================================================================

set -e

BASE_DIR="/path/to/data_root/ASE"
WORK_DIR="$BASE_DIR"

echo "========================================================"
echo "ASE data processing pipeline - full version"
echo "========================================================"
echo "Working directory: $WORK_DIR"
echo "Start time: $(date)"
echo ""

echo "========================================================"
echo "Step 1: generating aligned files in batch"
echo "========================================================"
echo ""

cd "$WORK_DIR"

SAMPLES=()
for file in "${BASE_DIR}"/counts_for_R_pipeline_T[0-9]*_with_GQ_geneBody.csv; do
    if [ -f "$file" ]; then
        basename=$(basename "$file")
        sample=$(echo "$basename" | sed 's/counts_for_R_pipeline_//; s/_with_GQ_geneBody\.csv//')
        SAMPLES+=("$sample")
    fi
done

echo "Found ${#SAMPLES[@]} tumor samples to process"
echo ""

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "Error: no tumor-sample counts files found!"
    exit 1
fi

SUCCESS=0
FAILED=0

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE"

    OUTPUT_FILE1="${BASE_DIR}/counts_for_R_pipeline_${SAMPLE}_with_GQ_geneBody_aligned.csv"
    OUTPUT_FILE2="${BASE_DIR}/counts_for_R_pipeline_${SAMPLE}_with_GQ_geneBody.aligned.csv"

    if [ -f "$OUTPUT_FILE1" ] || [ -f "$OUTPUT_FILE2" ]; then
        echo "  ✓ aligned file already exists, skipping"
        SUCCESS=$((SUCCESS + 1))
        continue
    fi

    python3 $HOME/bin/make_aligned_csv.py \
        --sample "$SAMPLE" \
        --base-dir "$BASE_DIR" \
        --out-dir "$BASE_DIR" && {
        echo "  ✓ success"
        SUCCESS=$((SUCCESS + 1))
    } || {
        echo "  ✗ failed"
        FAILED=$((FAILED + 1))
    }
done

echo ""
echo "Step 1 done: $SUCCESS succeeded, $FAILED failed"
echo ""

echo "========================================================"
echo "Step 2: checking aligned file names"
echo "========================================================"
echo ""

shopt -s nullglob
for file in "${BASE_DIR}"/counts_for_R_pipeline_T*_with_GQ_geneBody.aligned.csv; do
    if [ -f "$file" ]; then
        newname=$(echo "$file" | sed 's/_with_GQ_geneBody\.aligned\.csv/_with_GQ_geneBody_aligned.csv/')
        echo "Renaming: $(basename "$file") -> $(basename "$newname")"
        mv "$file" "$newname"
    fi
done
shopt -u nullglob

echo ""
echo "Current aligned files:"
ls -lh "${BASE_DIR}"/counts_for_R_pipeline_T*_aligned.csv 2>/dev/null | awk '{print $9, $5}' || echo "no aligned files found"
echo ""

ALIGNED_COUNT=0
shopt -s nullglob
aligned_files=("${BASE_DIR}"/counts_for_R_pipeline_T*_aligned.csv)
ALIGNED_COUNT=${#aligned_files[@]}
shopt -u nullglob

echo "Found $ALIGNED_COUNT aligned files"
echo ""

if [ $ALIGNED_COUNT -eq 0 ]; then
    echo "Error: no aligned files found!"
    exit 1
fi

echo "========================================================"
echo "Step 3: integrating tumor-sample ASE data"
echo "========================================================"
echo ""

cd "$WORK_DIR"

if [ -f "/path/to/data_root/wasp/code/integrate_SNPs_with_ASE.R" ]; then
    echo "Running the R script..."
    Rscript /path/to/data_root/wasp/code/integrate_SNPs_with_ASE.R
    echo ""
else
    echo "Warning: integrate_T_samples_ASE_aligned.R not found"
    echo "Make sure the script is in the current directory: $WORK_DIR"
fi

echo ""
echo "========================================================"
echo "Processing complete!"
echo "========================================================"
echo "End time: $(date)"
echo ""

if [ -f "${WORK_DIR}/all_T_samples_SNPs_with_ASE_annotation.csv" ]; then
    echo "Output files:"
    ls -lh "${WORK_DIR}"/all_T_samples_*.csv "${WORK_DIR}"/significant_ASE_genes_list.csv "${WORK_DIR}"/T_samples_SNP_ASE_summary.csv 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
    echo ""
    echo "✓ Data processing pipeline finished!"
else
    echo "Warning: expected output file not found"
    echo "Check the R script log"
fi
