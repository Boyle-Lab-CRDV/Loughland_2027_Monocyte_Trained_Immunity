#!/bin/bash

#SBATCH --job-name="merge_peaks"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G                     
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=zuleima.pava@burnet.edu.au
#SBATCH --time=00-24:00:00

set -euo pipefail

module load miniconda
conda activate atac_seq_env

# Path where input folders are located
INPUT_BASE="/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis"
OUTPUT_BASE="/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/featurecounts"
mkdir -p "$OUTPUT_BASE"

# List of subfolders to process (space-separated)
FOLDERS=("chmi_1_m" "chmi_1_t" "sabah_1_m" "sabah_1_t" "sabah_2_m" "sabah_2_t" "chmi_2_m" "chmi_2_t")

# Loop through each folder
for folder in "${FOLDERS[@]}"; do
    echo "Processing folder: $folder"
    
    INPUT_DIR="$INPUT_BASE/$folder"
    OUTPUT_DIR="$OUTPUT_BASE/${folder}_featurecounts"
    mkdir -p "$OUTPUT_DIR"

    cd "$INPUT_DIR"

    echo "Merging peak files..."
    zcat *.narrowPeak.gz | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > "$OUTPUT_DIR/merged_peaks.bed"

    echo "Generating all_peaks.bed and counting overlaps..."
    zcat *.narrowPeak.gz > "$OUTPUT_DIR/all_peaks.bed"
    bedtools subtract -a "$OUTPUT_DIR/merged_uni_peaks.bed" -b "$OUTPUT_DIR/all_peaks.bed" -c > "$OUTPUT_DIR/merged_uni_peaks_with_counts.bed"

#    echo "Filtering peaks present in >= 4 samples..."
#    awk '$4 >= 4' "$OUTPUT_DIR/merged_uni_peaks_with_counts.bed" | cut -f1-3 > "$OUTPUT_DIR/filtered_peaks.bed"

    echo "Converting to SAF format..."
    awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} \
         {print "peak_"NR, $1, $2, $3, "."}' "$OUTPUT_DIR/merged_uni_peaks_with_counts.bed" > "$OUTPUT_DIR/merged_uni_peaks_with_counts.saf"

    echo "Running featureCounts..."
    featureCounts -p -a "$OUTPUT_DIR/merged_uni_peaks_with_counts.saf" -F SAF \
                  -o "$OUTPUT_DIR/uni_peak_counts.txt" -T 16 "$INPUT_DIR"/*.bam

    echo "Done with $folder"
done
