#!/bin/bash

# ============================================================
# NOMBRE: 04_qiime_import_and_trim.sh
# DESCRIPCIÓN: Importa FASTQs y remueve primers con Cutadapt (Paired y Single-End)
# USO: bash 0_preprocessing/scripts/04_qiime_import_and_trim.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

echo "RUNNING:  bash 0_preprocessing/scripts/03b_qiime_import_and_trim.sh"
set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."

CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
QIIME2_DIR="$ROOT_DIR/data/pre_processed/qiime2_workspace"
MANIFESTS_DIR="$QIIME2_DIR/manifests"

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

mkdir -p "$QIIME2_DIR" "$MANIFESTS_DIR"

# Leer configuración e iterar sanitizando variables
tail -n +2 "$CONFIG_FILE" | while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do
    
    # Limpieza estricta de espacios y \r
    STUDY=$(echo "$STUDY" | tr -d '\r' | xargs)
    FWD=$(echo "$FWD" | tr -d '\r' | xargs)
    REV=$(echo "$REV" | tr -d '\r' | xargs)
    REGION=$(echo "$REGION" | tr -d '\r' | xargs)

    [[ -z "$STUDY" || "$STUDY" =~ ^# ]] && continue

    RAW_FASTQ="$ROOT_DIR/data/raw/$STUDY/fastq"
    STUDY_DIR="$QIIME2_DIR/$STUDY"
    MANIFEST="$MANIFESTS_DIR/${STUDY}_manifest.tsv"

    mkdir -p "$STUDY_DIR"

    echo "=========================================="
    echo "🔄 Procesando importación y Cutadapt para: $STUDY ($REGION)"
    echo "=========================================="

    if [ ! -d "$RAW_FASTQ" ]; then
        echo "❌ ERROR: No existe el directorio FASTQ: $RAW_FASTQ"
        continue
    fi

    # Detectar si es Paired-End buscando _2.fastq.gz
    IS_PAIRED=false
    if compgen -G "$RAW_FASTQ/*_2.fastq.gz" > /dev/null; then
        IS_PAIRED=true
    fi

    # 1. Generación de Manifest
    if [ ! -f "$MANIFEST" ]; then
        echo "📄 Generando manifest (Paired-End: $IS_PAIRED)..."
        
        if [ "$IS_PAIRED" = true ]; then
            printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANIFEST"
            for f in "$RAW_FASTQ"/*_1.fastq.gz; do
                if [ -f "$f" ]; then
                    SAMPLE=$(basename "$f" _1.fastq.gz)
                    FWD_PATH=$(realpath "$f")
                    REV_PATH=$(realpath "$RAW_FASTQ/${SAMPLE}_2.fastq.gz")
                    printf "%s\t%s\t%s\n" "$SAMPLE" "$FWD_PATH" "$REV_PATH" >> "$MANIFEST"
                fi
            done
        else
            printf "sample-id\tabsolute-filepath\n" > "$MANIFEST"
            for f in "$RAW_FASTQ"/*.fastq.gz; do
                if [ -f "$f" ]; then
                    SAMPLE=$(basename "$f")
                    SAMPLE="${SAMPLE%.fastq.gz}"
                    SAMPLE="${SAMPLE%_1}"
                    FWD_PATH=$(realpath "$f")
                    printf "%s\t%s\n" "$SAMPLE" "$FWD_PATH" >> "$MANIFEST"
                fi
            done
        fi

        if [ ! -s "$MANIFEST" ] || [ $(wc -l < "$MANIFEST") -lt 2 ]; then
            echo "❌ ERROR: Manifest para $STUDY está vacío. Revisa: $RAW_FASTQ"
            exit 1
        fi
        echo "✅ Manifest creado en: $MANIFEST"
    fi

    # 2. Importación a QIIME 2
    if [ ! -f "$STUDY_DIR/demux.qza" ]; then
        echo "📥 Importando a QIIME 2..."
        if [ "$IS_PAIRED" = true ]; then
            qiime tools import \
                --type 'SampleData[PairedEndSequencesWithQuality]' \
                --input-path "$MANIFEST" \
                --output-path "$STUDY_DIR/demux.qza" \
                --input-format PairedEndFastqManifestPhred33V2
        else
            qiime tools import \
                --type 'SampleData[SequencesWithQuality]' \
                --input-path "$MANIFEST" \
                --output-path "$STUDY_DIR/demux.qza" \
                --input-format SingleEndFastqManifestPhred33V2
        fi
    fi

    # 3. Cutadapt
    if [ ! -f "$STUDY_DIR/demux_trimmed.qza" ]; then
        echo "✂️ Recortando primers con Cutadapt..."
        if [ "$IS_PAIRED" = true ]; then
            qiime cutadapt trim-paired \
                --i-demultiplexed-sequences "$STUDY_DIR/demux.qza" \
                --p-front-f "$FWD" \
                --p-front-r "$REV" \
                --p-cores 4 \
                --o-trimmed-sequences "$STUDY_DIR/demux_trimmed.qza" \
                --quiet
        else
            qiime cutadapt trim-single \
                --i-demultiplexed-sequences "$STUDY_DIR/demux.qza" \
                --p-front "$FWD" \
                --p-cores 4 \
                --o-trimmed-sequences "$STUDY_DIR/demux_trimmed.qza" \
                --quiet
        fi
        echo "✅ Cutadapt completado para $STUDY"
    else
        echo "✅ $STUDY ya cuenta con demux_trimmed.qza"
    fi
done

# #!/bin/bash
# 
# # ============================================================
# # NOMBRE: 03b_qiime_import_and_trim.sh
# # DESCRIPCIÓN: Importa FASTQs y remueve primers con Cutadapt leyendo desde config/datasets_config.tsv
# # USO: bash 0_preprocessing/scripts/03b_qiime_import_and_trim.sh
# # ============================================================
# 
# echo "RUNNING:  bash 0_preprocessing/scripts/03b_qiime_import_and_trim.sh"
# set -e
# 
# SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# ROOT_DIR="$SCRIPT_DIR/../.."
# 
# CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
# QIIME2_DIR="$ROOT_DIR/data/processed/qiime2_workspace"
# MANIFESTS_DIR="$QIIME2_DIR/manifests"
# 
# echo "=========================================="
# echo "Verificando entorno informático y configuración..."
# echo "=========================================="
# 
# if ! command -v qiime &> /dev/null; then
#     echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno con: conda activate qiime2-amplicon-2024.5"
#     exit 1
# fi
# 
# if [ ! -f "$CONFIG_FILE" ]; then
#     echo "❌ ERROR: Archivo de configuración no encontrado en: $CONFIG_FILE"
#     exit 1
# fi
# 
# mkdir -p "$QIIME2_DIR" "$MANIFESTS_DIR"
# 
# # Leer configuración saltando la cabecera
# tail -n +2 "$CONFIG_FILE" | while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION; do
#     
#     # Omitir líneas vacías o comentadas (#)
#     [[ -z "$STUDY" || "$STUDY" =~ ^# ]] && continue
# 
#     RAW_FASTQ="$ROOT_DIR/data/raw/$STUDY/fastq"
#     STUDY_DIR="$QIIME2_DIR/$STUDY"
#     MANIFEST="$MANIFESTS_DIR/${STUDY}_manifest.tsv"
# 
#     mkdir -p "$STUDY_DIR"
# 
#     echo "=========================================="
#     echo "🔄 Procesando importación y Cutadapt para: $STUDY ($REGION)"
#     echo "=========================================="
# 
#     # 1. Manifest
#     if [ ! -f "$MANIFEST" ]; then
#         echo "📄 Generando manifest..."
#         printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANIFEST"
#         for f in "$RAW_FASTQ"/*_1.fastq.gz; do
#             if [ -f "$f" ]; then
#                 SAMPLE=$(basename "$f" _1.fastq.gz)
#                 FWD_PATH=$(realpath "$f")
#                 REV_PATH=$(realpath "$RAW_FASTQ/${SAMPLE}_2.fastq.gz")
#                 printf "%s\t%s\t%s\n" "$SAMPLE" "$FWD_PATH" "$REV_PATH" >> "$MANIFEST"
#             fi
#         done
# 
#         # Validar Manifest
#         if [ ! -s "$MANIFEST" ] || [ $(wc -l < "$MANIFEST") -lt 2 ]; then
#             echo "❌ ERROR: Manifest para $STUDY está vacío o incompleto. Revisa la ruta: $RAW_FASTQ"
#             exit 1
#         fi
#     fi
# 
#     # 2. Import
#     if [ ! -f "$STUDY_DIR/demux.qza" ]; then
#         echo "📥 Importando a QIIME 2..."
#         qiime tools import \
#             --type 'SampleData[PairedEndSequencesWithQuality]' \
#             --input-path "$MANIFEST" \
#             --output-path "$STUDY_DIR/demux.qza" \
#             --input-format PairedEndFastqManifestPhred33V2
#     fi
# 
#     # 3. Cutadapt
#     if [ ! -f "$STUDY_DIR/demux_trimmed.qza" ]; then
#         echo "✂️ Recortando primers con Cutadapt (FWD: $FWD | REV: $REV)..."
#         qiime cutadapt trim-paired \
#             --i-demultiplexed-sequences "$STUDY_DIR/demux.qza" \
#             --p-front-f "$FWD" \
#             --p-front-r "$REV" \
#             --p-cores 4 \
#             --o-trimmed-sequences "$STUDY_DIR/demux_trimmed.qza" \
#             --quiet
#         echo "✅ Cutadapt completado para $STUDY"
#     else
#         echo "✅ $STUDY ya cuenta con demux_trimmed.qza"
#     fi
# done

# #======================================================
# #======================================================
# #======================================================
################ OlDER VERSION
# #======================================================
# #======================================================
# #======================================================
# #======================================================
# 
# #======================================================
# #!/bin/bash
# 
# # ============================================================
# # NOMBRE: 03b_qiime_import_and_trim.sh
# # DESCRIPCIÓN: Importa FASTQs y remueve primers con Cutadapt por dataset
# # USO: bash 0_preprocessing/scripts/03b_qiime_import_and_trim.sh
# # ============================================================
# 
# echo "=========================================="
# echo "Running 03b_import_and_trim.sh"
# echo "=========================================="
# 
# 
# set -e
# 
# SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# ROOT_DIR="$SCRIPT_DIR/../.."
# 
# QIIME2_DIR="$ROOT_DIR/data/processed/qiime2_workspace"
# MANIFESTS_DIR="$QIIME2_DIR/manifests"
# 
# echo "=========================================="
# echo "Verificando entorno informático..."
# echo "=========================================="
# 
# if ! command -v qiime &> /dev/null; then
#     echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno con: conda activate qiime2-amplicon-2024.5"
#     exit 1
# fi
# 
# mkdir -p "$QIIME2_DIR" "$MANIFESTS_DIR"
# 
# 
# ### ------------------------------------------------------------
# ## DATASET INFO
# ### ------------------------------------------------------------
# 
# DATASETS=("chen_2022" "haimov_2022")
# 
# FWD_PRIMERS=(
#   "CCTAYGGGRBGCASCAG"   # Chen (V3-V4)
#   "GTGCCAGCMGCCGCGGTAA" # Haimov (V4)
# )
# 
# REV_PRIMERS=(
#   "GGACTACNNGGGTATCTAAT"  # Chen
#   "GGACTACHVGGGTWTCTAAT" # Haimov
# 
# )
# 
# for i in "${!DATASETS[@]}"; do
#     STUDY="${DATASETS[$i]}"
#     FWD="${FWD_PRIMERS[$i]}"
#     REV="${REV_PRIMERS[$i]}"
# 
#     RAW_FASTQ="$ROOT_DIR/data/raw/$STUDY/fastq"
#     STUDY_DIR="$QIIME2_DIR/$STUDY"
#     MANIFEST="$MANIFESTS_DIR/${STUDY}_manifest.tsv"
# 
#     mkdir -p "$STUDY_DIR"
# 
#     echo "=========================================="
#     echo "🔄 Procesando importación y Cutadapt para: $STUDY"
#     echo "=========================================="
# 
#     # 1. Manifest
#     if [ ! -f "$MANIFEST" ]; then
#         echo "📄 Generando manifest..."
#         printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANIFEST"
#         for f in "$RAW_FASTQ"/*_1.fastq.gz; do
#             if [ -f "$f" ]; then
#                 SAMPLE=$(basename "$f" _1.fastq.gz)
#                 FWD_PATH=$(realpath "$f")
#                 REV_PATH=$(realpath "$RAW_FASTQ/${SAMPLE}_2.fastq.gz")
#                 printf "%s\t%s\t%s\n" "$SAMPLE" "$FWD_PATH" "$REV_PATH" >> "$MANIFEST"
#             fi
#         done
# 
#         # Validar Manifest
#         if [ ! -s "$MANIFEST" ] || [ $(wc -l < "$MANIFEST") -lt 2 ]; then
#             echo "❌ ERROR: Manifest para $STUDY está vacío. Revisa la ruta: $RAW_FASTQ"
#             exit 1
#         fi
#     fi
# 
#     # 2. Import
#     if [ ! -f "$STUDY_DIR/demux.qza" ]; then
#         echo "📥 Importando a QIIME 2..."
#         qiime tools import \
#             --type 'SampleData[PairedEndSequencesWithQuality]' \
#             --input-path "$MANIFEST" \
#             --output-path "$STUDY_DIR/demux.qza" \
#             --input-format PairedEndFastqManifestPhred33V2
#     fi
# 
#     # 3. Cutadapt
#     if [ ! -f "$STUDY_DIR/demux_trimmed.qza" ]; then
#         echo "✂️ Recortando primers con Cutadapt..."
#         qiime cutadapt trim-paired \
#             --i-demultiplexed-sequences "$STUDY_DIR/demux.qza" \
#             --p-front-f "$FWD" \
#             --p-front-r "$REV" \
#             --p-cores 4 \
#             --o-trimmed-sequences "$STUDY_DIR/demux_trimmed.qza" \
#             --quiet
#         echo "✅ Cutadapt completado para $STUDY"
#     else
#         echo "✅ $STUDY ya cuenta con demux_trimmed.qza"
#     fi
# done
# 
# 
