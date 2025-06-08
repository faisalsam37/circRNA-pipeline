
#!/bin/bash
set -euo pipefail


# Load config
source config.sh


# Create  output directories
mkdir -p "$FASTQ_DIR" "$TRIMMED_DIR" "$FASTQC_DIR" "$ALIGN_DIR" "$CIRCRNA_DIR" "$LOG_DIR"

# User Input
echo "Available circRNA detection tools:"
echo "  1. circRNA_finder"
echo "  2. find_circ"
echo "  3. clear_quant"
echo "  4. ciri2"
echo ""

# Keep asking until a valid tool is entered
while true; do
  read -p "Enter the name of the circRNA detection tool to use: " CIRC_TOOL
  if [[ "$CIRC_TOOL" == "circRNA_finder" || "$CIRC_TOOL" == "find_circ" || "$CIRC_TOOL" == "clear_quant" || "$CIRC_TOOL" == "ciri2" ]]; then
    break
  else
    echo "Wrong circRNA detection tool, choose one from the list"
  fi
done

read -p "Enter the SRA accession number (e.g., SRR1234567): " SRA_ID


# Setup logging 
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$LOG_DIR/${SRA_ID}_${CIRC_TOOL}_$TIMESTAMP.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "Logging to: $LOG_FILE"


# Check if FASTQ files already exist before downloading and converting
if [[ -f "$FASTQ_DIR/${SRA_ID}_1.fastq.gz" && -f "$FASTQ_DIR/${SRA_ID}_2.fastq.gz" ]] || [[ -f "$FASTQ_DIR/${SRA_ID}.fastq.gz" ]]; then
  echo "FASTQ files for $SRA_ID already exist in $FASTQ_DIR. Skipping download."
else
  echo "Downloading and converting SRA data for $SRA_ID directly with fasterq-dump..."

  if command -v fasterq-dump >/dev/null 2>&1; then
    echo "Using fasterq-dump..."
    fasterq-dump "$SRA_ID" -O "$FASTQ_DIR" --split-files --skip-technical || {
      echo "fasterq-dump failed, falling back to fastq-dump"
      fastq-dump --gzip --split-files -O "$FASTQ_DIR" "$SRA_ID"
    }

    # gzip manually after fasterq-dump
    gzip "$FASTQ_DIR/${SRA_ID}_1.fastq" "$FASTQ_DIR/${SRA_ID}_2.fastq"
  else
    echo "fasterq-dump not found, using fastq-dump..."
    fastq-dump --gzip --split-files -O "$FASTQ_DIR" "$SRA_ID"
  fi

fi

# Force rename downloaded FASTQ files to standardized names
echo "Renaming downloaded FASTQ files to standardized names..."

# Count how many fastq or fastq.gz files are in FASTQ_DIR related to this SRA_ID
FASTQ_FILES=($(ls "$FASTQ_DIR"/*"${SRA_ID}"*.fastq* 2>/dev/null || true))
NUM_FILES=${#FASTQ_FILES[@]}

if [[ $NUM_FILES -eq 0 ]]; then
  echo "Error: No FASTQ files found to rename!"
  exit 1
elif [[ $NUM_FILES -eq 1 ]]; then
  # Single-end read
  echo "Single-end FASTQ detected."
  mv -f "${FASTQ_FILES[0]}" "$FASTQ_DIR/${SRA_ID}.fastq.gz"
elif [[ $NUM_FILES -eq 2 ]]; then
  # Paired-end read
  echo "Paired-end FASTQ detected."
  # Rename first file to _1 and second to _2
  mv -f "${FASTQ_FILES[0]}" "$FASTQ_DIR/${SRA_ID}_1.fastq.gz"
  mv -f "${FASTQ_FILES[1]}" "$FASTQ_DIR/${SRA_ID}_2.fastq.gz"
else
  echo "Warning: More than two FASTQ files found, using first two for paired-end rename."
  mv -f "${FASTQ_FILES[0]}" "$FASTQ_DIR/${SRA_ID}_1.fastq.gz"
  mv -f "${FASTQ_FILES[1]}" "$FASTQ_DIR/${SRA_ID}_2.fastq.gz"
fi


# Step 2: Determine Read Type and Run FastQC 
if [[ -f "$FASTQ_DIR/${SRA_ID}_1.fastq.gz" && -f "$FASTQ_DIR/${SRA_ID}_2.fastq.gz" ]]; then
  MODE="PE"
  echo "Detected paired-end reads."
  fastqc "$FASTQ_DIR/${SRA_ID}_1.fastq.gz" "$FASTQ_DIR/${SRA_ID}_2.fastq.gz" -o "$FASTQC_DIR"

elif [[ -f "$FASTQ_DIR/${SRA_ID}.fastq.gz" ]]; then
  MODE="SE"
  SE_FASTQ="$FASTQ_DIR/${SRA_ID}.fastq.gz"
  echo "Detected single-end reads."
  fastqc "$SE_FASTQ" -o "$FASTQC_DIR"

else
  echo "Error: Could not determine read type for $SRA_ID."
  exit 1
fi

# Unzip
echo "Unzipping FastQC output..."
if compgen -G "$FASTQC_DIR"/*_fastqc.zip > /dev/null; then
  for zipfile in "$FASTQC_DIR"/*_fastqc.zip; do
    unzip -o "$zipfile" -d "$FASTQC_DIR" >/dev/null
  done
else
  echo "Warning: No FastQC ZIP files found to unzip!"
fi

if [[ "$MODE" == "PE" ]]; then
  SUMMARY1="$FASTQC_DIR/${SRA_ID}_1_fastqc/summary.txt"
  SUMMARY2="$FASTQC_DIR/${SRA_ID}_2_fastqc/summary.txt"
else
  SUMMARY1="$FASTQC_DIR/${SRA_ID}_fastqc/summary.txt"
fi


# Step 3: Determine Trimming Parameters Based on FastQC Report
echo "Analyzing FastQC summary to choose trimming strategy..."
if grep -q 'FAIL\|WARN' "$SUMMARY1" || { [[ "$MODE" == "PE" ]] && grep -q 'FAIL\|WARN' "$SUMMARY2"; }; then
  echo "Low quality detected. Using flexible trimming for small dataset."
  TRIM_PARAMS="SLIDINGWINDOW:4:10 LEADING:3 TRAILING:3 MINLEN:20"
else
  echo "Good quality detected. Using light trimming."
  TRIM_PARAMS="SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:36"
fi



# Unzip FastQC reports if zip files exist
for zipfile in "$FASTQC_DIR"/*.zip; do
  unzip -o "$zipfile" -d "$FASTQC_DIR" >/dev/null
done

# Step 4: Trimming with Trimmomatic
if [[ "$MODE" == "PE" ]]; then
  # Check if trimmed files already exist
  if [[ -f "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz" && -f "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz" ]]; then
    echo "Trimmed paired-end files for $SRA_ID already exist. Skipping trimming..."
  else
    echo "Trimming paired-end files for $SRA_ID..."
    java -jar "$TRIMMOMATIC_JAR" PE -threads 8 -phred33 \
      "$FASTQ_DIR/${SRA_ID}_1.fastq.gz" "$FASTQ_DIR/${SRA_ID}_2.fastq.gz" \
      "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz" "$TRIMMED_DIR/${SRA_ID}_1_unpaired.fastq.gz" \
      "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz" "$TRIMMED_DIR/${SRA_ID}_2_unpaired.fastq.gz" \
      ILLUMINACLIP:"$TRIM_ADAPTER_PE":2:30:10 $TRIM_PARAMS
  fi
else
  # Check if trimmed file already exists for single-end
  if [[ -f "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz" ]]; then
    echo "Trimmed single-end file for $SRA_ID already exists. Skipping trimming..."
  else
    echo "Trimming single-end file for $SRA_ID..."
    java -jar "$TRIMMOMATIC_JAR" SE -threads 8 -phred33 \
      "$FASTQ_DIR/${SRA_ID}.fastq.gz" "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz" \
      ILLUMINACLIP:"$TRIM_ADAPTER_SE":2:30:10 $TRIM_PARAMS
  fi
fi



# Step 5: Alignment and circRNA detection
ALIGN_OUT="$ALIGN_DIR/${CIRC_TOOL}_${SRA_ID}"
CIRC_OUT="$CIRCRNA_DIR/${CIRC_TOOL}/${SRA_ID}"
mkdir -p "$ALIGN_OUT" "$CIRC_OUT"

case "$CIRC_TOOL" in
  clear_quant)
    echo "Aligning with STAR for clear_quant..."
    if [[ "$MODE" == "PE" ]]; then
      STAR --genomeDir "$STAR_INDEX" \
        --readFilesIn "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz" "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz" \
        --runThreadN 8 --readFilesCommand gunzip -c \
        --outFileNamePrefix "$ALIGN_OUT/${SRA_ID}_" --outSAMtype BAM SortedByCoordinate
    else
      STAR --genomeDir "$STAR_INDEX" \
        --readFilesIn "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz" \
        --runThreadN 8 --readFilesCommand gunzip -c \
        --outFileNamePrefix "$ALIGN_OUT/${SRA_ID}_" --outSAMtype BAM SortedByCoordinate
    fi
    echo "Running clear_quant..."
    "$CLEAR_PATH/clear_quant" -1 "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz" ${MODE:+-2 "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz"} \
      -g "$REFERENCE_GENOME" -G "$GTF_FILE" -i "$HISAT2_INDEX" -j "$BOWTIE2_INDEX" -o "$CIRC_OUT" -p 8
    ;;

  circRNA_finder)
    echo "Aligning with STAR for circRNA_finder..."
    STAR_OUT="$ALIGN_OUT/star"
    mkdir -p "$STAR_OUT"
    if [[ "$MODE" == "PE" ]]; then
      "$CIRCRNA_FINDER_PATH/runStar.pl" --inFile1 "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz" --inFile2 "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz" --genomeDir "$STAR_INDEX" --outPrefix "$STAR_OUT/${SRA_ID}_"
    else
      "$CIRCRNA_FINDER_PATH/runStar.pl" --inFile1 "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz" --genomeDir "$STAR_INDEX" --outPrefix "$STAR_OUT/${SRA_ID}_"
    fi
    echo "Running postProcessStarAlignment.pl..."
    "$CIRCRNA_FINDER_PATH/postProcessStarAlignment.pl" --starDir "$STAR_OUT" --outDir "$CIRC_OUT" --minLen 100
    ;;

  find_circ)
    echo "Aligning with Bowtie2 for find_circ..."
    if [[ "$MODE" == "PE" ]]; then
      bowtie2 -p 8 -x "$BOWTIE2_INDEX" \
        -1 <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz") \
        -2 <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz") \
        | samtools view -bS - | samtools sort -n -o "$ALIGN_OUT/${SRA_ID}_sorted.bam"
    else
      bowtie2 -p 8 -x "$BOWTIE2_INDEX" \
        -U <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz") \
        | samtools view -bS - | samtools sort -n -o "$ALIGN_OUT/${SRA_ID}_sorted.bam"
    fi
    echo "Extracting unmapped reads and generating anchors..."
    samtools view -hf 4 "$ALIGN_OUT/${SRA_ID}_sorted.bam" | python2 "$FIND_CIRC_PATH/unmapped2anchors.py" - | gzip > "$ALIGN_OUT/${SRA_ID}_anchors.fastq.gz"
    echo "Running find_circ.py..."
    bowtie2 -p 8 -x "$BOWTIE2_INDEX" -U "$ALIGN_OUT/${SRA_ID}_anchors.fastq.gz" | "$FIND_CIRC_PATH/find_circ.py" --genome="$REFERENCE_GENOME" --name "$SRA_ID" > "$CIRC_OUT/splice_sites.bed"
    ;;

 ciri2)
  echo "Aligning with BWA for CIRI2..."
  if [[ "$MODE" == "PE" ]]; then
    bwa mem "$BWA_INDEX" <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_1_trimmed.fastq.gz") <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_2_trimmed.fastq.gz") > "$ALIGN_OUT/${SRA_ID}_ciri2.sam"
  else
    bwa mem "$BWA_INDEX" <(gunzip -c "$TRIMMED_DIR/${SRA_ID}_trimmed.fastq.gz") > "$ALIGN_OUT/${SRA_ID}_ciri2.sam"
  fi
  echo "Running CIRI2..."
  "$CIRI2_PATH/CIRI2.pl" -I "$ALIGN_OUT/${SRA_ID}_ciri2.sam" -O "$CIRC_OUT/ciri2_output.txt" -F "$REFERENCE_GENOME" -A "$GTF_FILE"
  ;;
*)
  echo "Error: Unsupported circRNA tool: $CIRC_TOOL"
  exit 1
  ;;
esac

echo "Pipeline complete for $SRA_ID using $CIRC_TOOL."
