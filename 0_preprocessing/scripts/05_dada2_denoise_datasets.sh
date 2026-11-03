# # ============================================================
# # ============================================================
# # ============================================================
  # # uses corrected qiime2 processed table from metadata
# # ============================================================
# # ============================================================
# # ============================================================


#!/bin/bash

# ============================================================
# NAME: 05_dada2_denoise_datasets.sh
# DESCRIPTION: Denoising and error correction with DADA2 per dataset
#              reading parameters dynamically from config/datasets_config.tsv
# USAGE: bash 0_preprocessing/scripts/05_dada2_denoise_datasets.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."

CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
QIIME2_DIR="$ROOT_DIR/data/pre_processed/qiime2_workspace"
THREADS=4  # Optimized for macOS ARM stability

echo "=========================================="
echo "Checking environment and configuration..."
echo "=========================================="

if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME 2 not detected. Activate environment: conda activate qiime2-amplicon-2024.5"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Configuration file not found at: $CONFIG_FILE"
    exit 1
fi

# Sanitize config file in place
tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"

process_study() {
    local STUDY_NAME=$1
    local TF=$2
    local TR=$3
    local LF=$4
    local LR=$5
    local REGION=$6

    local STUDY_DIR="$QIIME2_DIR/$STUDY_NAME"
    local RAW_FASTQ="$ROOT_DIR/data/raw/$STUDY_NAME/fastq"

    echo "------------------------------------------"
    echo "🔄 PROCESSING DADA2 FOR: $STUDY_NAME ($REGION)"
    echo "------------------------------------------"

    if [ -f "$STUDY_DIR/table.qza" ] && [ -f "$STUDY_DIR/rep_seqs.qza" ]; then
        echo "✅ $STUDY_NAME already has DADA2 completed. Skipping..."
        return 0
    fi

    if [ ! -f "$STUDY_DIR/demux_trimmed.qza" ]; then
        echo "❌ ERROR: Missing $STUDY_DIR/demux_trimmed.qza. Run import/trimming step first."
        exit 1
    fi

    # Determine if dataset is Paired-End or Single-End
    IS_PAIRED=false
    if compgen -G "$RAW_FASTQ/*_2.fastq.gz" > /dev/null && [ "$TR" -gt 0 ]; then
        IS_PAIRED=true
    fi

    if [ "$IS_PAIRED" = true ]; then
        echo "🧬 Denoising Paired-End (Trunc F: $TF | Trunc R: $TR | Trim F: $LF | Trim R: $LR)..."
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$STUDY_DIR/demux_trimmed.qza" \
            --p-trunc-len-f "$TF" \
            --p-trunc-len-r "$TR" \
            --p-trim-left-f "$LF" \
            --p-trim-left-r "$LR" \
            --o-table "$STUDY_DIR/table.qza" \
            --o-representative-sequences "$STUDY_DIR/rep_seqs.qza" \
            --o-denoising-stats "$STUDY_DIR/denoising_stats.qza" \
            --p-n-threads "$THREADS"
    else
        echo "🧬 Denoising Single-End (Trunc Len: $TF | Trim Left: $LF)..."
        qiime dada2 denoise-single \
            --i-demultiplexed-seqs "$STUDY_DIR/demux_trimmed.qza" \
            --p-trunc-len "$TF" \
            --p-trim-left "$LF" \
            --o-table "$STUDY_DIR/table.qza" \
            --o-representative-sequences "$STUDY_DIR/rep_seqs.qza" \
            --o-denoising-stats "$STUDY_DIR/denoising_stats.qza" \
            --p-n-threads "$THREADS"
    fi

    echo "✅ Dataset $STUDY_NAME processed successfully."
}

# Read TSV line by line safely using process substitution
while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do
    
    STUDY=$(echo "$STUDY" | xargs)
    TF=$(echo "$TF" | xargs)
    TR=$(echo "$TR" | xargs)
    LF=$(echo "$LF" | xargs)
    LR=$(echo "$LR" | xargs)
    REGION=$(echo "$REGION" | xargs)

    [[ -z "$STUDY" || "$STUDY" =~ ^# ]] && continue

    process_study "$STUDY" "$TF" "$TR" "$LF" "$LR" "$REGION"
done < <(tail -n +2 "$CONFIG_FILE")

echo "=========================================="
echo "🎉 All active datasets processed with DADA2."
echo "=========================================="


# # ============================================================
# # ============================================================
# # ============================================================
# # =this one works but does not use a qiime2 corrected table from metadata, uses raw metadata
# # ============================================================
# # ============================================================
# # ============================================================

#!/bin/bash

# ============================================================
# NOMBRE: 05_dada2_denoise_datasets.sh
# DESCRIPCIÓN: Denoising y corrección de errores con DADA2 por dataset 
#              leyendo parámetros dinámicamente desde config/datasets_config.tsv
# USO: bash 0_preprocessing/scripts/05_dada2_denoise_datasets.sh
# ============================================================

###### Rcpp fails often, correct inside environment...
### conda activate qiime2-amplicon-2024.5
### corrext: Rscript -e "install.packages('Rcpp', repos='https://cloud.r-project.org')"
### Rscript -e "library(dada2); cat('✅ DADA2 loaded successfully!\n')"

set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."

CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
QIIME2_DIR="$ROOT_DIR/data/pre_processed/qiime2_workspace"
THREADS=4  # Reducido a 4 para estabilidad de memoria en macOS ARM

echo "=========================================="
echo "Verificando entorno informático y configuración..."
echo "=========================================="

if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno con: conda activate qiime2-amplicon-2024.5"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Archivo de configuración no encontrado en: $CONFIG_FILE"
    exit 1
fi

# Sanitización automática en sitio (elimina \r de datasets_config.tsv)
tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"

process_study() {
    local STUDY_NAME=$1
    local TF=$2
    local TR=$3
    local LF=$4
    local LR=$5
    local REGION=$6

    local STUDY_DIR="$QIIME2_DIR/$STUDY_NAME"
    local RAW_FASTQ="$ROOT_DIR/data/raw/$STUDY_NAME/fastq"

    echo "------------------------------------------"
    echo "🔄 PROCESANDO DADA2 PARA: $STUDY_NAME ($REGION)"
    echo "------------------------------------------"

    if [ -f "$STUDY_DIR/table.qza" ] && [ -f "$STUDY_DIR/rep_seqs.qza" ]; then
        echo "✅ $STUDY_NAME ya tiene DADA2 completado. Saltando..."
        return 0
    fi

    if [ ! -f "$STUDY_DIR/demux_trimmed.qza" ]; then
        echo "❌ ERROR: Falta $STUDY_DIR/demux_trimmed.qza. Ejecuta el paso previo de importación/trimming."
        exit 1
    fi

    # Determinar si el dataset es Paired-End o Single-End
    IS_PAIRED=false
    if compgen -G "$RAW_FASTQ/*_2.fastq.gz" > /dev/null && [ "$TR" -gt 0 ]; then
        IS_PAIRED=true
    fi

    if [ "$IS_PAIRED" = true ]; then
        echo "🧬 Denoising Paired-End con DADA2 (Trunc F: $TF | Trunc R: $TR | Trim F: $LF | Trim R: $LR)..."
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs "$STUDY_DIR/demux_trimmed.qza" \
            --p-trunc-len-f "$TF" \
            --p-trunc-len-r "$TR" \
            --p-trim-left-f "$LF" \
            --p-trim-left-r "$LR" \
            --o-table "$STUDY_DIR/table.qza" \
            --o-representative-sequences "$STUDY_DIR/rep_seqs.qza" \
            --o-denoising-stats "$STUDY_DIR/denoising_stats.qza" \
            --p-n-threads "$THREADS"
    else
        echo "🧬 Denoising Single-End con DADA2 (Trunc Len: $TF | Trim Left: $LF)..."
        qiime dada2 denoise-single \
            --i-demultiplexed-seqs "$STUDY_DIR/demux_trimmed.qza" \
            --p-trunc-len "$TF" \
            --p-trim-left "$LF" \
            --o-table "$STUDY_DIR/table.qza" \
            --o-representative-sequences "$STUDY_DIR/rep_seqs.qza" \
            --o-denoising-stats "$STUDY_DIR/denoising_stats.qza" \
            --p-n-threads "$THREADS"
    fi

    echo "✅ Dataset $STUDY_NAME procesado correctamente."
}

# Leer TSV línea por línea omitiendo la cabecera
tail -n +2 "$CONFIG_FILE" | while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do
    
    # Sanitizar variables en memoria
    STUDY=$(echo "$STUDY" | xargs)
    TF=$(echo "$TF" | xargs)
    TR=$(echo "$TR" | xargs)
    LF=$(echo "$LF" | xargs)
    LR=$(echo "$LR" | xargs)
    REGION=$(echo "$REGION" | xargs)

    [[ -z "$STUDY" || "$STUDY" =~ ^# ]] && continue

    process_study "$STUDY" "$TF" "$TR" "$LF" "$LR" "$REGION"
done

echo "=========================================="
echo "🎉 Todos los datasets activos han sido procesados con DADA2."
echo "=========================================="

# # ============================================================
# # ============================================================
# # ============================================================
# # = OLDER SCRIPT
# # ============================================================
# # ============================================================
# # ============================================================

# 
# #!/bin/bash
# 
# # ============================================================
# # NOMBRE: 05_dada2_denoise_datasets.sh
# # DESCRIPCIÓN: Denoising y corrección de errores con DADA2 por dataset.
# # USO: bash 0_preprocessing/scripts/05_denoise_datasets.sh
# # ============================================================
# 
# echo "=========================================="
# echo "Running 05_dada2_denoise_datasets.sh"
# echo "=========================================="
# 
# set -e
# 
# SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# ROOT_DIR="$SCRIPT_DIR/../.."
# 
# QIIME2_DIR="$ROOT_DIR/data/processed/qiime2_workspace"
# THREADS=8
# 
# echo "=========================================="
# echo "Verificando entorno..."
# echo "=========================================="
# 
# if ! command -v qiime &> /dev/null; then
#     echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno con: conda activate qiime2-amplicon-2024.5"
#     exit 1
# fi
# 
# ### ------------------------------------------------------------
# ## DATASET INFO
# ### ------------------------------------------------------------
# 
# DATASETS=(
#   "chen_2022"
#   "haimov_2022"
#   "barone_2019"
# )
# 
# # Parámetros DADA2 por estudio [Chen, Haimov, Barone]
# # Ajustados tras inspeccionar 04_viz_quality_plot.sh (.qzv)
# TRUNC_LEN_F=(235 200 240)
# TRUNC_LEN_R=(230 150 220)
# TRIM_LEFT_F=(0 0 0)
# TRIM_LEFT_R=(0 0 0)
# 
# echo "=========================================="
# echo "Verificando entorno informático..."
# echo "=========================================="
# 
# if ! command -v qiime &> /dev/null; then
#     echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno: conda activate qiime2-amplicon-2024.5"
#     exit 1
# fi
# 
# #if [ ! -f "$SILVA_CLASSIFIER" ]; then
# #    echo "❌ ERROR: Clasificador SILVA no encontrado en: $SILVA_CLASSIFIER"
# #    echo "   Ejecuta primero: bash 0_preprocessing/scripts/03_download_silva_classifier.sh"
# #    exit 1
# #fi
# 
# #mkdir -p "$QIIME2_DIR" "$MANIFESTS_DIR" "$EXPORT_DIR"
# process_study() {
#     local STUDY_NAME=$1
#     local TF=$2
#     local TR=$3
#     local LF=$4
#     local LR=$5
# 
#     local STUDY_DIR="$QIIME2_DIR/$STUDY_NAME"
# 
#     echo "------------------------------------------"
#     echo "🔄 PROCESANDO DADA2 PARA: $STUDY_NAME"
#     echo "   Truncado F: $TF | Truncado R: $TR"
#     echo "------------------------------------------"
# 
#     if [ -f "$STUDY_DIR/table.qza" ] && [ -f "$STUDY_DIR/rep_seqs.qza" ]; then
#         echo "✅ $STUDY_NAME ya tiene DADA2 completado. Saltando..."
#         return 0
#     fi
# 
#     if [ ! -f "$STUDY_DIR/demux_trimmed.qza" ]; then
#         echo "❌ ERROR: Falta $STUDY_DIR/demux_trimmed.qza. Ejecuta 03b_import_and_trim.sh primero."
#         exit 1
#     fi
# 
#     echo "🧬 Denoising con DADA2..."
#     qiime dada2 denoise-paired \
#         --i-demultiplexed-seqs "$STUDY_DIR/demux_trimmed.qza" \
#         --p-trunc-len-f "$TF" \
#         --p-trunc-len-r "$TR" \
#         --p-trim-left-f "$LF" \
#         --p-trim-left-r "$LR" \
#         --o-table "$STUDY_DIR/table.qza" \
#         --o-representative-sequences "$STUDY_DIR/rep_seqs.qza" \
#         --o-denoising-stats "$STUDY_DIR/denoising_stats.qza" \
#         --p-n-threads "$THREADS" \
#         --verbose
# 
#     echo "✅ Dataset $STUDY_NAME procesado correctamente."
# }
# 
# for i in "${!DATASETS[@]}"; do
#     process_study "${DATASETS[$i]}" "${TRUNC_LEN_F[$i]}" "${TRUNC_LEN_R[$i]}" "${TRIM_LEFT_F[$i]}" "${TRIM_LEFT_R[$i]}"
# done
# 
# echo "=========================================="
# echo "🎉 Todos los datasets han sido denoised por separado."
# echo "Siguiente paso: Ejecuta 'bash 0_preprocessing/scripts/06_merge_cluster_export.sh'"
# echo "=========================================="