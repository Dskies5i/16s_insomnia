#!/bin/bash

# ============================================================
# NOMBRE: 04_viz_quality_plot.sh
# Resumen: Genera resúmenes .qzv de calidad para inspección visual 
# USO: 
#     1. conda activate qiime2-amplicon-2024.5
#     2. bash 0_preprocessing/scripts/04_viz_quality_plot.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

echo 'Running 0_preprocessing/scripts/04_viz_quality_plot.sh'
set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."

CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
QIIME2_DIR="$ROOT_DIR/data/processed/qiime2_workspace"

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

echo "=========================================="
echo "📊 GENERANDO ARCHIVOS .QZV DE CALIDAD"
echo "=========================================="

# Process each dataset line with sanitization and fallback for missing trailing newline
tail -n +2 "$CONFIG_FILE" | while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do
    
    # Strip carriage returns and outer whitespace
    STUDY=$(echo "$STUDY" | tr -d '\r' | xargs)

    [[ -z "$STUDY" || "$STUDY" =~ ^# ]] && continue

    STUDY_DIR="$QIIME2_DIR/$STUDY"
    
    if [ -f "$STUDY_DIR/demux_trimmed.qzv" ]; then
        echo "✅ $STUDY ya cuenta con demux_trimmed.qzv"
    elif [ -f "$STUDY_DIR/demux_trimmed.qza" ]; then
        echo "📊 Resumiendo calidad para $STUDY..."
        qiime demux summarize \
          --i-data "$STUDY_DIR/demux_trimmed.qza" \
          --o-visualization "$STUDY_DIR/demux_trimmed.qzv"
        echo "✅ Generado: $STUDY_DIR/demux_trimmed.qzv"
    else
        echo "⚠️  Atención: $STUDY_DIR/demux_trimmed.qza no existe. Ejecuta Cutadapt primero."
    fi
done



#### ========================================
## OLDER VERSION
#### ========================================

#!/bin/bash

# ============================================================
# NOMBRE: 04_viz_quality_plot.sh
# Resumen: Genera resúmenes .qzv de calidad para inspección visual 
#              - para ajustar correctamente los parametros de qiime2 para cada estudio, es necesario vizualizar la calidad por nucleótido (Phred Score / Q-score) a lo largo de cada ciclo de secuenciación.
#              - https://view.qiime2.org/
#              - Esa gráfica de calidad es la que te permite calcular con precisión matemática cuánto debes truncar sin romper el ensamble de las lecturas (merge).
#              - Calidad de las bases (Eje Y): Muestra el valor de Phred Score (Q). Buscas que la mediana de calidad se mantenga por encima de Q25 - Q30 (99.7% – 99.9% de precisión).
#              - Posición en el amplicón (Eje X): Los ciclos de secuenciación de inicio a fin (ej. 1 a 250 pb).
#              - Caída al final de la lectura: Verás que hacia los últimos nucleótidos (extremo 3'), la calidad decae. Esto ocurre siempre de forma más marcada en la lectura reversa (R2).
#              - El objetivo es cortar las bases con baja calidad en el extremo 3' mediante --p-trunc-len-f y --p-trunc-len-r, pero asegurando suficiente solapamiento para que DADA2 pueda unir R1 y R2.
#             
# USO: 
#    1. conda activate qiime2-amplicon-2024.5
#    2. bash 0_preprocessing/scripts/04_viz_quality_plot.sh
# ============================================================
# set -e
# 
# SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# ROOT_DIR="$SCRIPT_DIR/../.."
# QIIME2_DIR="$ROOT_DIR/data/processed/qiime2_workspace"
# DATASETS=("chen_2022" "haimov_2022")
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
# 
# echo "=========================================="
# echo "📊 GENERANDO ARCHIVOS .QZV DE CALIDAD"
# echo "=========================================="
# 
# for STUDY in "${DATASETS[@]}"; do
#     STUDY_DIR="$QIIME2_DIR/$STUDY"
#     if [ -f "$STUDY_DIR/demux_trimmed.qza" ]; then
#         echo "📊 Resumiendo calidad para $STUDY..."
#         qiime demux summarize \
#           --i-data "$STUDY_DIR/demux_trimmed.qza" \
#           --o-visualization "$STUDY_DIR/demux_trimmed.qzv"
#         echo "✅ Generado: $STUDY_DIR/demux_trimmed.qzv"
#     else
#         echo "⚠️  Atención: $STUDY_DIR/demux_trimmed.qza no existe. Ejecuta Cutadapt primero."
#     fi
# done