#!/bin/bash

# ============================================================
# NOMBRE: 02_metadata_qiime2_format.sh
# DESCRIPCIÓN: Localiza los metadatos CSV detallados generados en el paso 01b 
#              y ejecuta el script de Python encargado de formatearlos e 
#              inicializarlos para su compatibilidad con QIIME 2.
#
# PROPÓSITO DEL FLUJO:
#   - Lee dinámicamente los estudios activos desde `datasets_config.tsv`.
#   - Procesa cada archivo CSV crudo a través del motor de formateo.
#   - Inserta la cabecera especial `#q2:types` para que QIIME 2 identifique 
#     correctamente variables numéricas (como la edad) y categóricas 
#     (como los grupos de estudio) durante las pruebas estadísticas.
#
# USO: bash 0_preprocessing/scripts/02_metadata_qiime2_format.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================
set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."

CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
METADATA_OUT_DIR="$ROOT_DIR/data/pre_processed/metadata"
PYTHON_SCRIPT="$SCRIPT_DIR/utils/02_metadata_qiime2_format_engine.py"

mkdir -p "$METADATA_OUT_DIR"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: No se encontró el archivo de configuración: $CONFIG_FILE"
    exit 1
fi

tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"

echo "=========================================="
echo "📄 Generando archivos de metadata para QIIME 2..."
echo "=========================================="

while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do

    STUDY=$(echo "$STUDY" | xargs)
    [[ -z "$STUDY" || "$STUDY" =~ ^# ]] && continue

    # Updated to match 01b output path
    RAW_META_CSV="$ROOT_DIR/data/raw/$STUDY/${STUDY}_metadata_detailed.csv"

    # Clean output name without extra suffixes
    QIIME2_META_TSV="$METADATA_OUT_DIR/${STUDY}_metadata_q2.tsv"

    echo "------------------------------------------"
    echo "🔍 Procesando metadata para: $STUDY"
    echo "------------------------------------------"

    if [ ! -f "$RAW_META_CSV" ]; then
        echo "⚠️ WARNING: No se encontró $RAW_META_CSV. Verifique el paso 01b."
        continue
    fi

    python3 "$PYTHON_SCRIPT" "$RAW_META_CSV" "$QIIME2_META_TSV"

done < <(tail -n +2 "$CONFIG_FILE")

echo "=========================================="
echo "🎉 Todos los archivos de metadata formateados en: $METADATA_OUT_DIR"
echo "=========================================="