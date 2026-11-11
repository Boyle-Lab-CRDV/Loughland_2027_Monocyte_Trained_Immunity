#!/bin/bash

#SBATCH --job-name="fastq_QC"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G                     
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=zuleima.pava@burnet.edu.au
#SBATCH --time=00-24:00:00

module load miniconda
conda activate atac_seq_env

BASE_DIR="/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/QC_160226"
FASTQ_CSV="$BASE_DIR/all_fastq_paths2.csv"
REFERENCE="/working_groups/boylelab/shared/ZulyPava/formating4hpc/reference/GRCh38_noalt_as"
OUTDIR="$BASE_DIR/bam"
THREADS=16

mkdir -p "$OUTDIR"

# Set to 1 only if a downstream tool truly requires name-sorted BAM
MAKE_NAME_SORTED=0

echo "Reading samples from: $FASTQ_CSV"
echo "Output dir: $OUTDIR"
echo "Threads: $THREADS"

# Skip header, read: sample_id,r1,r2
tail -n +2 "$FASTQ_CSV" | while IFS=, read -r SAMPLE R1 R2; do
  # Trim possible Windows CR characters
  SAMPLE="${SAMPLE//$'\r'/}"
  R1="${R1//$'\r'/}"
  R2="${R2//$'\r'/}"

  # Basic checks
  if [[ -z "$SAMPLE" || -z "$R1" || -z "$R2" ]]; then
    echo "Skipping malformed line (need sample,r1,r2)"
    continue
  fi
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: FASTQ not found for $SAMPLE"
    echo "  R1: $R1"
    echo "  R2: $R2"
    exit 1
  fi

  echo "=== $SAMPLE ==="

  SORTED_BAM="$OUTDIR/${SAMPLE}.sorted.bam"
  NAME_BAM="$OUTDIR/${SAMPLE}.namesorted.bam"
  LOG="$OUTDIR/${SAMPLE}.bowtie2.log"

  # Align -> BAM -> coordinate sort (no SAM on disk)
  bowtie2 -x "$REFERENCE" -1 "$R1" -2 "$R2" \
    --no-unal --very-sensitive -p "$THREADS" 2> "$LOG" \
  | samtools view -b - \
  | samtools sort --threads "$THREADS" -o "$SORTED_BAM" -

  samtools index "$SORTED_BAM"

  if [[ "$MAKE_NAME_SORTED" -eq 1 ]]; then
    samtools sort -n --threads "$THREADS" -o "$NAME_BAM" "$SORTED_BAM"
  fi

  echo "Done: $SAMPLE"
done

echo "All samples complete."



