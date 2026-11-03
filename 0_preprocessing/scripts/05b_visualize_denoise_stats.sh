#!/bin/bash

# ============================================================
# NOMBRE: 05b_visualize_denoise_stats.sh
# DESCRIPCIÓN: Generación automática de visualizaciones QIIME 2 (.qzv) 
#              para evaluar denoising_stats y feature-tables por dataset
# USO: bash 0_preprocessing/scripts/05b_visualize_denoise_stats.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."

CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
QIIME2_DIR="$ROOT_DIR/data/pre_processed/qiime2_workspace"
METADATA_DIR="$ROOT_DIR/data/pre_processed/metadata"

echo "=========================================="
echo "Verificando entorno y configuración..."
echo "=========================================="

if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno con: conda activate qiime2-amplicon-2024.5"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Archivo de configuración no encontrado en: $CONFIG_FILE"
    exit 1
fi

# Sanitización de caracteres de fin de línea (\r) de Windows
tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"

echo "=========================================="
echo "📊 Generando visualizaciones (.qzv) para datasets activos..."
echo "=========================================="

while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do
    
    # Limpiar espacios en blanco
    STUDY=$(echo "$STUDY" | xargs)
    
    # Ignorar líneas vacías, comentarios (#) o encabezados
    [[ -z "$STUDY" || "$STUDY" =~ ^# || "$STUDY" == "dataset_id" || "$STUDY" == "study_id" ]] && continue

    STUDY_DIR="$QIIME2_DIR/$STUDY"
    META_FILE="$METADATA_DIR/${STUDY}_metadata_q2.tsv"
    
    echo "------------------------------------------"
    echo "🔍 Procesando visualizaciones para: $STUDY"
    echo "------------------------------------------"

    if [ ! -d "$STUDY_DIR" ]; then
        echo "⚠️ WARNING: La carpeta $STUDY_DIR no existe. Saltando..."
        continue
    fi

    # 1. Visualización de estadísticas de Denoising (DADA2)
    if [ -f "$STUDY_DIR/denoising_stats.qza" ]; then
        if [ -f "$STUDY_DIR/denoising_stats.qzv" ]; then
            echo "  ⏩ denoising_stats.qzv ya existe. Saltando..."
        else
            echo "  📈 Generando denoising_stats.qzv..."
            qiime metadata tabulate \
                --m-input-file "$STUDY_DIR/denoising_stats.qza" \
                --o-visualization "$STUDY_DIR/denoising_stats.qzv"
        fi
    else
        echo "  ⚠️ Falta denoising_stats.qza en $STUDY_DIR"
    fi

    # 2. Resumen de la Tabla de Abundancias (Feature Table with Metadata)
    if [ -f "$STUDY_DIR/table.qza" ]; then
        if [ -f "$STUDY_DIR/table.qzv" ]; then
            echo "  ⏩ table.qzv ya existe. Saltando..."
        else
            echo "  📋 Generando table.qzv..."
            
            # Check if metadata TSV exists to attach it
            if [ -f "$META_FILE" ]; then
                qiime feature-table summarize \
                    --i-table "$STUDY_DIR/table.qza" \
                    --m-sample-metadata-file "$META_FILE" \
                    --o-visualization "$STUDY_DIR/table.qzv"
            else
                echo "  ⚠️ Metadata no encontrada en $META_FILE. Generando table.qzv sin metadata..."
                qiime feature-table summarize \
                    --i-table "$STUDY_DIR/table.qza" \
                    --o-visualization "$STUDY_DIR/table.qzv"
            fi
        fi
    else
        echo "  ⚠️ Falta table.qza en $STUDY_DIR"
    fi

    # 3. Resumen de Secuencias Representativas
    if [ -f "$STUDY_DIR/rep_seqs.qza" ]; then
        if [ -f "$STUDY_DIR/rep_seqs.qzv" ]; then
            echo "  ⏩ rep_seqs.qzv ya existe. Saltando..."
        else
            echo "  🧬 Generando rep_seqs.qzv..."
            qiime feature-table tabulate-seqs \
                --i-data "$STUDY_DIR/rep_seqs.qza" \
                --o-visualization "$STUDY_DIR/rep_seqs.qzv"
        fi
    else
        echo "  ⚠️ Falta rep_seqs.qza en $STUDY_DIR"
    fi

    echo "  ✅ Estado actualizado para: $STUDY"

done < <(tail -n +2 "$CONFIG_FILE")

echo "=========================================="
echo "🎉 Proceso de visualización completado."
echo "Puedes arrastrar los archivos .qzv generados a https://view.qiime2.org"
echo "=========================================="

















# #!/bin/bash
# 
# # ============================================================
# # NOMBRE: 05b_visualize_denoise_stats.sh
# # DESCRIPCIÓN: Generación automática de visualizaciones QIIME 2 (.qzv) 
# #              para evaluar denoising_stats y feature-tables por dataset
# # USO: bash 0_preprocessing/scripts/05b_visualize_denoise_stats.sh
# # ============================================================
# 
# set -e
# 
# SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# ROOT_DIR="$SCRIPT_DIR/../.."
# 
# CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
# QIIME2_DIR="$ROOT_DIR/data/processed/qiime2_workspace"
# 
# echo "=========================================="
# echo "Verificando entorno y configuración..."
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
# # Sanitización de caracteres de fin de línea (\r) de Windows
# tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"
# 
# echo "=========================================="
# echo "📊 Generando visualizaciones (.qzv) para datasets activos..."
# echo "=========================================="
# 
# while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do
#     
#     # Limpiar espacios en blanco
#     STUDY=$(echo "$STUDY" | xargs)
#     
#     # Ignorar líneas vacías, comentarios (#) o encabezados
#     [[ -z "$STUDY" || "$STUDY" =~ ^# || "$STUDY" == "dataset_id" || "$STUDY" == "study_id" ]] && continue
# 
#     STUDY_DIR="$QIIME2_DIR/$STUDY"
#     
#     echo "------------------------------------------"
#     echo "🔍 Procesando visualizaciones para: $STUDY"
#     echo "------------------------------------------"
# 
#     if [ ! -d "$STUDY_DIR" ]; then
#         echo "⚠️ WARNING: La carpeta $STUDY_DIR no existe. Saltando..."
#         continue
#     fi
# 
#     # 1. Visualización de estadísticas de Denoising (DADA2)
#     if [ -f "$STUDY_DIR/denoising_stats.qza" ]; then
#         if [ -f "$STUDY_DIR/denoising_stats.qzv" ]; then
#             echo "  ⏩ denoising_stats.qzv ya existe. Saltando..."
#         else
#             echo "  📈 Generando denoising_stats.qzv..."
#             qiime metadata tabulate \
#                 --m-input-file "$STUDY_DIR/denoising_stats.qza" \
#                 --o-visualization "$STUDY_DIR/denoising_stats.qzv"
#         fi
#     else
#         echo "  ⚠️ Falta denoising_stats.qza en $STUDY_DIR"
#     fi
# 
#     # 2. Resumen de la Tabla de Abundancias (Feature Table)
#     if [ -f "$STUDY_DIR/table.qza" ]; then
#         if [ -f "$STUDY_DIR/table.qzv" ]; then
#             echo "  ⏩ table.qzv ya existe. Saltando..."
#         else
#             echo "  📋 Generando table.qzv..."
#             qiime feature-table summarize \
#                 --i-table "$STUDY_DIR/table.qza" \
#                 --o-visualization "$STUDY_DIR/table.qzv"
#         fi
#     else
#         echo "  ⚠️ Falta table.qza en $STUDY_DIR"
#     fi
# 
#     # 3. Resumen de Secuencias Representativas
#     if [ -f "$STUDY_DIR/rep_seqs.qza" ]; then
#         if [ -f "$STUDY_DIR/rep_seqs.qzv" ]; then
#             echo "  ⏩ rep_seqs.qzv ya existe. Saltando..."
#         else
#             echo "  🧬 Generando rep_seqs.qzv..."
#             qiime feature-table tabulate-seqs \
#                 --i-data "$STUDY_DIR/rep_seqs.qza" \
#                 --o-visualization "$STUDY_DIR/rep_seqs.qzv"
#         fi
#     else
#         echo "  ⚠️ Falta rep_seqs.qza en $STUDY_DIR"
#     fi
# 
#     echo "  ✅ Estado actualizado para: $STUDY"
# 
# done < "$CONFIG_FILE"
# 
# echo "=========================================="
# echo "🎉 Proceso de visualización completado."
# echo "Puedes arrastrar los archivos .qzv generados a https://view.qiime2.org"
# echo "=========================================="