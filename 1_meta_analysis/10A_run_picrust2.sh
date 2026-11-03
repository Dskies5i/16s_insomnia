#!/usr/bin/env bash
# NAME: 1_meta_analysis/10A_run_picrust2.sh
#
# COMMENT: PICRUSt2 must run BEFORE 10_function_taxa_ko_cor.R.
#          It produces outputs/tables/merged_KO_metagenome.tsv,
#          the KO prediction table script 10 needs.
#
# IMPORTANT: PICRUSt2 finishes → merged_KO_metagenome.tsv exists → script 10 reads it
#
# USE:
#    conda activate qiime2-amplicon-2024.5
#    bash 1_meta_analysis/10A_run_picrust2.sh

set -euo pipefail

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  Starting PICRUSt2 Execution Pipeline"
echo "════════════════════════════════════════════════════════════"

# ---- 1. Safety: must be in the right conda env -------------------
if [[ "${CONDA_DEFAULT_ENV:-}" != "qiime2-amplicon-2024.5" ]]; then
  echo "❌ Error: this script must be run from the qiime2-amplicon-2024.5 env."
  echo "    Run:  conda activate qiime2-amplicon-2024.5"
  exit 1
fi
echo "✅ Conda env: ${CONDA_DEFAULT_ENV}"

# ---- 2. Install PICRUSt2 if not present --------------------------
if ! command -v picrust2_pipeline.py &> /dev/null; then
  echo "PICRUSt2 not found in current env. Installing via Conda..."
  conda install -c bioconda -c conda-forge picrust2 -y
else
  echo "✅ PICRUSt2 is already installed."
fi

# ---- 3. Check and prepare required input files -------------------
FASTA="data/pre_processed/unified/unified_rep_seqs.fasta"
BIOM="data/pre_processed/unified/feature-table.biom"
TSV="data/pre_processed/unified/unified_otu_table.tsv"

if [[ ! -f "$FASTA" ]]; then
  echo "❌ Error: missing required input FASTA: $FASTA"
  exit 1
fi

# Auto-convert TSV to BIOM if BIOM is missing but TSV exists
if [[ ! -f "$BIOM" ]]; then
  if [[ -f "$TSV" ]]; then
    echo "⚠️  feature-table.biom not found, but unified_otu_table.tsv exists."
    echo "    Converting TSV to BIOM format automatically..."
    biom convert -i "$TSV" -o "$BIOM" --table-type="OTU table" --to-binary
    echo "✅ Successfully converted TSV to BIOM."
  else
    echo "❌ Error: neither feature-table.biom nor unified_otu_table.tsv found."
    exit 1
  fi
fi

echo "✅ Input files verified:"
echo "    FASTA: $FASTA"
echo "    BIOM : $BIOM"

# ---- 4. Sanity check: BIOM should have samples -------------------
NUM_SAMPLES=$(biom summarize-table -i "$BIOM" \
  | awk '/Num samples:/ {print $3}' | tr -d ',')
echo "    BIOM samples: ${NUM_SAMPLES}"
if [[ "${NUM_SAMPLES}" -lt 1 ]]; then
  echo "❌ Error: BIOM reports 0 samples. Refusing to run."
  exit 1
fi

# ---- Checkpoint: Skip if output already exists -------------------
FINAL_OUTPUT="outputs/tables/merged_KO_metagenome.tsv"

if [[ -f "$FINAL_OUTPUT" ]]; then
  echo ""
  echo "⏩ [SALTADO] PICRUSt2 output already exists at: $FINAL_OUTPUT"
  echo "   (Se omite la computación pesada de PICRUSt2 y se continúa)."
  echo ""
  exit 0
fi

# ---- 5. Run PICRUSt2 ---------------------------------------------
echo ""
echo "Running PICRUSt2 pipeline (this may take 15–60 minutes)..."
echo ""

picrust2_pipeline.py \
  -s "$FASTA" \
  -i "$BIOM" \
  -o picrust2_out \
  -p 4

# ---- 6. Automatically find and move the KO prediction table ------
echo ""
echo "Formatting and moving output..."
mkdir -p outputs/tables

# Search dynamically for the unstratified KO prediction table inside picrust2_out
FOUND_KO=$(find picrust2_out -name "pred_metagenome_unstrat.tsv.gz" | head -n 1)

if [[ -z "$FOUND_KO" || ! -f "$FOUND_KO" ]]; then
  echo "❌ Error: Could not find pred_metagenome_unstrat.tsv.gz anywhere inside picrust2_out/"
  exit 1
fi

echo "Found KO file at: $FOUND_KO"
cp "$FOUND_KO" outputs/tables/pred_metagenome_unstrat.tsv.gz
gunzip -f outputs/tables/pred_metagenome_unstrat.tsv.gz
mv outputs/tables/pred_metagenome_unstrat.tsv outputs/tables/merged_KO_metagenome.tsv

echo ""
echo "✅ Successfully created outputs/tables/merged_KO_metagenome.tsv!"
echo ""
head -2 outputs/tables/merged_KO_metagenome.tsv

# ---- 7. Verify Sample ID Alignment -------------------------------
echo ""
echo "════════════════════════════════════════════════════════════"
echo "  Verifying Sample ID Alignment between KO & Phyloseq"
echo "════════════════════════════════════════════════════════════"

echo "First 3 Sample IDs in merged_KO_metagenome.tsv:"
head -1 outputs/tables/merged_KO_metagenome.tsv | tr '\t' '\n' | sed -n '2,4p'

if [[ -f "outputs/models/physeq_filtered.rds" ]]; then
  echo ""
  echo "First 3 Sample IDs in physeq_filtered.rds:"
  Rscript -e 'ps <- readRDS("outputs/models/physeq_filtered.rds"); cat(head(sample_names(ps), 3), sep = "\n")'
else
  echo "⚠️ Warning: outputs/models/physeq_filtered.rds not found to cross-check sample IDs."
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  PICRUSt2 complete. Next:"
echo "    Rscript 1_meta_analysis/10_function_taxa_ko_cor.R"
echo "════════════════════════════════════════════════════════════"





# # OLD WORKING VERSION WITHOUT SKIP CHECK
# #!/usr/bin/env bash
# # NAME: 1_meta_analysis/10A_run_picrust2.sh
# #
# # COMMENT: PICRUSt2 must run BEFORE 10_function_taxa_ko_cor.R.
# #          It produces outputs/tables/merged_KO_metagenome.tsv,
# #          the KO prediction table script 10 needs.
# #
# # IMPORTANT: PICRUSt2 finishes → merged_KO_metagenome.tsv exists → script 10 reads it
# #
# # USE:
# #   conda activate qiime2-amplicon-2024.5
# #   bash 1_meta_analysis/10A_run_picrust2.sh
# 
# set -euo pipefail
# 
# echo ""
# echo "════════════════════════════════════════════════════════════"
# echo "  Starting PICRUSt2 Execution Pipeline"
# echo "════════════════════════════════════════════════════════════"
# 
# # ---- 1. Safety: must be in the right conda env -------------------
# if [[ "${CONDA_DEFAULT_ENV:-}" != "qiime2-amplicon-2024.5" ]]; then
#   echo "❌ Error: this script must be run from the qiime2-amplicon-2024.5 env."
#   echo "    Run:  conda activate qiime2-amplicon-2024.5"
#   exit 1
# fi
# echo "✅ Conda env: ${CONDA_DEFAULT_ENV}"
# 
# # ---- 2. Install PICRUSt2 if not present --------------------------
# if ! command -v picrust2_pipeline.py &> /dev/null; then
#   echo "PICRUSt2 not found in current env. Installing via Conda..."
#   conda install -c bioconda -c conda-forge picrust2 -y
# else
#   echo "✅ PICRUSt2 is already installed."
# fi
# 
# # ---- 3. Check and prepare required input files -------------------
# FASTA="data/pre_processed/unified/unified_rep_seqs.fasta"
# BIOM="data/pre_processed/unified/feature-table.biom"
# TSV="data/pre_processed/unified/unified_otu_table.tsv"
# 
# if [[ ! -f "$FASTA" ]]; then
#   echo "❌ Error: missing required input FASTA: $FASTA"
#   exit 1
# fi
# 
# # Auto-convert TSV to BIOM if BIOM is missing but TSV exists
# if [[ ! -f "$BIOM" ]]; then
#   if [[ -f "$TSV" ]]; then
#     echo "⚠️  feature-table.biom not found, but unified_otu_table.tsv exists."
#     echo "    Converting TSV to BIOM format automatically..."
#     biom convert -i "$TSV" -o "$BIOM" --table-type="OTU table" --to-binary
#     echo "✅ Successfully converted TSV to BIOM."
#   else
#     echo "❌ Error: neither feature-table.biom nor unified_otu_table.tsv found."
#     exit 1
#   fi
# fi
# 
# echo "✅ Input files verified:"
# echo "    FASTA: $FASTA"
# echo "    BIOM : $BIOM"
# 
# # ---- 4. Sanity check: BIOM should have samples -------------------
# NUM_SAMPLES=$(biom summarize-table -i "$BIOM" \
#   | awk '/Num samples:/ {print $3}' | tr -d ',')
# echo "    BIOM samples: ${NUM_SAMPLES}"
# if [[ "${NUM_SAMPLES}" -lt 1 ]]; then
#   echo "❌ Error: BIOM reports 0 samples. Refusing to run."
#   exit 1
# fi
# 
# # ---- 5. Run PICRUSt2 ---------------------------------------------
# echo ""
# echo "Running PICRUSt2 pipeline (this may take 15–60 minutes)..."
# echo ""
# 
# picrust2_pipeline.py \
#   -s "$FASTA" \
#   -i "$BIOM" \
#   -o picrust2_out \
#   -p 4
# 
# # ---- 6. Automatically find and move the KO prediction table ------
# echo ""
# echo "Formatting and moving output..."
# mkdir -p outputs/tables
# 
# # Search dynamically for the unstratified KO prediction table inside picrust2_out
# FOUND_KO=$(find picrust2_out -name "pred_metagenome_unstrat.tsv.gz" | head -n 1)
# 
# if [[ -z "$FOUND_KO" || ! -f "$FOUND_KO" ]]; then
#   echo "❌ Error: Could not find pred_metagenome_unstrat.tsv.gz anywhere inside picrust2_out/"
#   exit 1
# fi
# 
# echo "Found KO file at: $FOUND_KO"
# cp "$FOUND_KO" outputs/tables/pred_metagenome_unstrat.tsv.gz
# gunzip -f outputs/tables/pred_metagenome_unstrat.tsv.gz
# mv outputs/tables/pred_metagenome_unstrat.tsv outputs/tables/merged_KO_metagenome.tsv
# 
# echo ""
# echo "✅ Successfully created outputs/tables/merged_KO_metagenome.tsv!"
# echo ""
# head -2 outputs/tables/merged_KO_metagenome.tsv
# 
# # ---- 7. Verify Sample ID Alignment -------------------------------
# echo ""
# echo "════════════════════════════════════════════════════════════"
# echo "  Verifying Sample ID Alignment between KO & Phyloseq"
# echo "════════════════════════════════════════════════════════════"
# 
# echo "First 3 Sample IDs in merged_KO_metagenome.tsv:"
# head -1 outputs/tables/merged_KO_metagenome.tsv | tr '\t' '\n' | sed -n '2,4p'
# 
# if [[ -f "outputs/models/physeq_filtered.rds" ]]; then
#   echo ""
#   echo "First 3 Sample IDs in physeq_filtered.rds:"
#   Rscript -e 'ps <- readRDS("outputs/models/physeq_filtered.rds"); cat(head(sample_names(ps), 3), sep = "\n")'
# else
#   echo "⚠️ Warning: outputs/models/physeq_filtered.rds not found to cross-check sample IDs."
# fi
# 
# echo ""
# echo "════════════════════════════════════════════════════════════"
# echo "  PICRUSt2 complete. Next:"
# echo "    Rscript 1_meta_analysis/10_function_taxa_ko_cor.R"
# echo "════════════════════════════════════════════════════════════"
# 



















