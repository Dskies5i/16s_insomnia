#!/bin/bash

# ============================================================
# NAME: 01b_download_study_sras.sh
# DESCRIPTION: Dynamically downloads FASTQ files using SRA Toolkit
#              and detailed metadata using pysradb
# USAGE: bash 0_preprocessing/scripts/01_download_study_sras_w_metadata.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

set -e
export PATH="/usr/local/sra-toolkit/bin:$PATH"

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."

CONFIG_FILE="$SCRIPT_DIR/../config/sra_projects.tsv"
THREADS=8

echo "=========================================="
echo "Step 1: Checking dependencies..."
echo "=========================================="
if ! command -v prefetch &> /dev/null; then
    echo "❌ ERROR: SRA Toolkit (prefetch) not found!"
    exit 1
fi
if ! command -v fasterq-dump &> /dev/null; then
    echo "❌ ERROR: SRA Toolkit (fasterq-dump) not found!"
    exit 1
fi
if ! command -v pysradb &> /dev/null; then
    echo "❌ ERROR: pysradb not found! Install it via 'pip install pysradb'."
    exit 1
fi
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Config file not found at: $CONFIG_FILE"
    exit 1
fi
echo "✅ All dependencies and config file detected."

# ------------------------------------------------------------
# Sanitization (removes \r from sra_projects.tsv)
# ------------------------------------------------------------
echo "=========================================="
echo "Step 2: Sanitizing config file ($CONFIG_FILE)..."
echo "=========================================="
tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"
echo "✅ Config file cleaned of carriage returns (\r)."

process_project() {
    local DATASET_NAME=$1
    local PROJECT=$2

    local RAW_DATA_DIR="$ROOT_DIR/data/raw/$DATASET_NAME"
    local FASTQ_DIR="$RAW_DATA_DIR/fastq"
    local SRR_LIST="$RAW_DATA_DIR/SRR_Acc_List.txt"
    local METADATA_CSV="$RAW_DATA_DIR/${DATASET_NAME}_metadata_detailed.csv"

    echo "=========================================="
    echo "🔄 Processing dataset: $DATASET_NAME ($PROJECT)"
    echo "=========================================="

    mkdir -p "$FASTQ_DIR"

    # ------------------------------------------------------------
    # 1. Fetch detailed metadata via pysradb
    # ------------------------------------------------------------
    echo "Fetching detailed metadata via pysradb for $PROJECT..."
    pysradb metadata "$PROJECT" --detailed --saveto "$METADATA_CSV"

    if [ ! -s "$METADATA_CSV" ]; then
        echo "❌ ERROR: pysradb failed to generate metadata for $PROJECT. Skipping."
        return 0
    fi
    echo "✅ Detailed metadata saved at: $METADATA_CSV"

    # ------------------------------------------------------------
    # 2. Extract SRR run accessions from pysradb CSV
    # ------------------------------------------------------------
    python3 -c "
import pandas as pd
df = pd.read_csv('$METADATA_CSV')
runs = df['run_accession'].dropna().unique()
with open('$SRR_LIST', 'w') as f:
    for r in sorted(runs):
        f.write(f'{r}\n')
"

    if [ ! -s "$SRR_LIST" ]; then
        echo "❌ ERROR: Could not parse run accessions from $METADATA_CSV. Skipping."
        return 0
    fi

    TOTAL_RUNS=$(wc -l < "$SRR_LIST" | tr -d ' ')
    echo "✅ Found $TOTAL_RUNS runs to process for $DATASET_NAME."

    # ------------------------------------------------------------
    # 3. Download FASTQ Files
    # ------------------------------------------------------------
    while read -r run_id <&3; do
        run_id=$(echo "$run_id" | tr -d '\r' | xargs)
        [ -z "$run_id" ] && continue

        # Skip if compressed FASTQ files already exist
        if ls "$FASTQ_DIR/${run_id}"*.fastq.gz &>/dev/null; then
            echo "⏭️  $run_id already downloaded and compressed. Skipping."
            continue
        fi

        echo "--- Downloading $run_id ($(date +%H:%M:%S)) ---"


        SRA_FILE="$RAW_DATA_DIR/$run_id/$run_id.sra"

        if [ -f "$SRA_FILE" ]; then
            fasterq-dump "$SRA_FILE" --split-files -e "$THREADS" -p -O "$FASTQ_DIR"
            rm -rf "$RAW_DATA_DIR/$run_id"
        else
            fasterq-dump "$run_id" --split-files -e "$THREADS" -p -O "$FASTQ_DIR"
        fi

        gzip -f "$FASTQ_DIR/${run_id}"*.fastq

    done 3< "$SRR_LIST"

    echo "✅ Dataset $DATASET_NAME processed successfully."
}

# Read TSV line by line safely using File Descriptor 4
while IFS=$'\t' read -r DATASET PROJECT <&4 || [ -n "$DATASET" ]; do

    DATASET=$(echo "$DATASET" | tr -d '\r' | xargs)
    PROJECT=$(echo "$PROJECT" | tr -d '\r' | xargs)

    [[ -z "$DATASET" || "$DATASET" == "dataset"* || "$DATASET" =~ ^# ]] && continue

    process_project "$DATASET" "$PROJECT"
done 4< "$CONFIG_FILE"

# #############################################################################
# #############################################################################
# #THIS WORKS - BASELINE
# 
# #############################################################################
# #############################################################################
# #############################################################################
# 
# 
# #!/bin/bash
# 
# # ============================================================
# # NAME: 01b_download_study_sras.sh
# # DESCRIPTION: Dynamically downloads FASTQ files and metadata from NCBI SRA
# #              using NCBI E-utilities API + SRA Toolkit
# # USAGE: bash 0_preprocessing/scripts/01_download_study_sras.sh
# # ============================================================
# 
# set -e
# export PATH="/usr/local/sra-toolkit/bin:$PATH"
# 
# SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# ROOT_DIR="$SCRIPT_DIR/../.."
# 
# CONFIG_FILE="$SCRIPT_DIR/../config/sra_projects.tsv"
# THREADS=8
# 
# echo "=========================================="
# echo "Step 1: Checking dependencies..."
# echo "=========================================="
# if ! command -v prefetch &> /dev/null; then
#     echo "❌ ERROR: SRA Toolkit (prefetch) not found!"
#     exit 1
# fi
# if ! command -v fasterq-dump &> /dev/null; then
#     echo "❌ ERROR: SRA Toolkit (fasterq-dump) not found!"
#     exit 1
# fi
# if ! command -v python3 &> /dev/null; then
#     echo "❌ ERROR: python3 not found!"
#     exit 1
# fi
# if [ ! -f "$CONFIG_FILE" ]; then
#     echo "❌ ERROR: Config file not found at: $CONFIG_FILE"
#     exit 1
# fi
# echo "✅ All dependencies and config file detected."
# 
# # ------------------------------------------------------------
# # Sanitización automática en sitio (elimina \r de sra_projects.tsv)
# # ------------------------------------------------------------
# echo "=========================================="
# echo "Step 2: Sanitizing config file ($CONFIG_FILE)..."
# echo "=========================================="
# tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"
# echo "✅ Config file cleaned of carriage returns (\r)."
# 
# process_project() {
#     local DATASET_NAME=$1
#     local PROJECT=$2
# 
#     local RAW_DATA_DIR="$ROOT_DIR/data/raw/$DATASET_NAME"
#     local FASTQ_DIR="$RAW_DATA_DIR/fastq"
#     local SRR_LIST="$RAW_DATA_DIR/SRR_Acc_List.txt"
#     local METADATA_CSV="$RAW_DATA_DIR/metadata.csv"
# 
#     echo "=========================================="
#     echo "🔄 Processing dataset: $DATASET_NAME ($PROJECT)"
#     echo "=========================================="
# 
#     mkdir -p "$FASTQ_DIR"
# 
#     echo "Fetching run list and metadata from NCBI API for $PROJECT..."
# 
#     # 1. Fetch UIDs from ESearch via stdin pipe to python
#     SEARCH_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${PROJECT}%5BBioProject%5D&retmode=json&retmax=1000"
# 
#     UID_LIST=$(curl -s "$SEARCH_URL" | python3 -c "
# import sys, json
# try:
#     data = json.load(sys.stdin)
#     id_list = data.get('esearchresult', {}).get('idlist', [])
#     print(','.join(id_list))
# except Exception:
#     pass
# ")
# 
#     # Fallback if [BioProject] tag returns empty
#     if [ -z "$UID_LIST" ]; then
#         SEARCH_URL_FALLBACK="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${PROJECT}&retmode=json&retmax=1000"
#         UID_LIST=$(curl -s "$SEARCH_URL_FALLBACK" | python3 -c "
# import sys, json
# try:
#     data = json.load(sys.stdin)
#     id_list = data.get('esearchresult', {}).get('idlist', [])
#     print(','.join(id_list))
# except Exception:
#     pass
# ")
#     fi
# 
#     if [ -z "$UID_LIST" ]; then
#         echo "❌ ERROR: No SRA entries found for $PROJECT. Skipping."
#         return 0
#     fi
# 
#     # 2. Fetch Run Accessions & Generate CSV Metadata directly from ESummary JSON
#     SUMMARY_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=${UID_LIST}&retmode=json"
# 
#     curl -s "$SUMMARY_URL" | python3 -c "
# import sys, json, csv, re
# 
# try:
#     data = json.load(sys.stdin)
#     result = data.get('result', {})
#     uids = result.get('uids', [])
# 
#     accs = set()
#     metadata_rows = []
# 
#     for u in uids:
#         item = result.get(u, {})
#         runs_str = item.get('runs', '')
#         found_runs = re.findall(r'([S|E|D]RR[0-9]+)', runs_str)
# 
#         expself = item.get('expsxml', '')
#         biosample = re.search(r'BioSample=\"([^\"]+)\"', expself)
#         biosample_id = biosample.group(1) if biosample else ''
# 
#         sample_name = re.search(r'SAMN[0-9]+', expself)
#         sample_acc = sample_name.group(0) if sample_name else ''
# 
#         exp_acc = re.search(r'([S|E|D]RX[0-9]+)', expself)
#         exp_id = exp_acc.group(1) if exp_acc else ''
# 
#         organism = re.search(r'ScientificName=\"([^\"]+)\"', expself)
#         org_name = organism.group(1) if organism else ''
# 
#         for run in found_runs:
#             accs.add(run)
#             metadata_rows.append({
#                 'Run': run,
#                 'BioProject': '$PROJECT',
#                 'BioSample': biosample_id,
#                 'SampleName': sample_acc,
#                 'Experiment': exp_id,
#                 'Organism': org_name,
#                 'SRA_UID': u
#             })
# 
#     # Output SRR Acc List
#     with open('$SRR_LIST', 'w') as f_acc:
#         for a in sorted(accs):
#             f_acc.write(f'{a}\n')
# 
#     # Output CSV Metadata
#     if metadata_rows:
#         with open('$METADATA_CSV', 'w', newline='') as f_csv:
#             writer = csv.DictWriter(f_csv, fieldnames=['Run', 'BioProject', 'BioSample', 'SampleName', 'Experiment', 'Organism', 'SRA_UID'])
#             writer.writeheader()
#             # Deduplicate by Run accession
#             seen_runs = set()
#             for row in metadata_rows:
#                 if row['Run'] not in seen_runs:
#                     writer.writerow(row)
#                     seen_runs.add(row['Run'])
# 
# except Exception as e:
#     sys.exit(1)
# "
# 
#     if [ ! -s "$SRR_LIST" ]; then
#         echo "❌ ERROR: Could not parse run accessions for $PROJECT. Skipping."
#         return 0
#     fi
# 
#     TOTAL_RUNS=$(wc -l < "$SRR_LIST" | tr -d ' ')
#     echo "✅ Found $TOTAL_RUNS runs to process for $DATASET_NAME."
# 
#     if [ -s "$METADATA_CSV" ]; then
#         echo "✅ Metadata successfully generated at: $METADATA_CSV"
#     else
#         echo "⚠️ Warning: Metadata file was not generated."
#     fi
# 
#     # ------------------------------------------------------------
#     # 3. Download FASTQ Files
#     # ------------------------------------------------------------
#     while read -r run_id <&3; do
#         run_id=$(echo "$run_id" | tr -d '\r' | xargs)
#         [ -z "$run_id" ] && continue
# 
#         # Skip if compressed paired/single files already exist
#         if ls "$FASTQ_DIR/${run_id}"*.fastq.gz &>/dev/null; then
#             echo "⏭️  $run_id already downloaded and compressed. Skipping."
#             continue
#         fi
# 
#         echo "--- Downloading $run_id ($(date +%H:%M:%S)) ---"
# 
#         prefetch "$run_id" --output-directory "$RAW_DATA_DIR"
# 
#         SRA_FILE="$RAW_DATA_DIR/$run_id/$run_id.sra"
# 
#         if [ -f "$SRA_FILE" ]; then
#             fasterq-dump "$SRA_FILE" --split-files -e "$THREADS" -p -O "$FASTQ_DIR"
#             rm -rf "$RAW_DATA_DIR/$run_id"
#         else
#             fasterq-dump "$run_id" --split-files -e "$THREADS" -p -O "$FASTQ_DIR"
#         fi
# 
#         gzip -f "$FASTQ_DIR/${run_id}"*.fastq
# 
#     done 3< "$SRR_LIST"
# 
#     echo "✅ Dataset $DATASET_NAME processed successfully."
# }
# 
# # Read TSV line by line safely using File Descriptor 4
# while IFS=$'\t' read -r DATASET PROJECT <&4 || [ -n "$DATASET" ]; do
# 
#     DATASET=$(echo "$DATASET" | tr -d '\r' | xargs)
#     PROJECT=$(echo "$PROJECT" | tr -d '\r' | xargs)
# 
#     [[ -z "$DATASET" || "$DATASET" == "dataset"* || "$DATASET" =~ ^# ]] && continue
# 
#     process_project "$DATASET" "$PROJECT"
# done 4< "$CONFIG_FILE"
# 
