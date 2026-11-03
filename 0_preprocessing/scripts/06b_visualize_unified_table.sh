#!/bin/bash

# ============================================================
# NOMBRE: 07_visualize_unified_table.sh
# DESCRIPCIÓN: Genera archivos de visualización (.qzv) para evaluar la 
#              integridad de la tabla unificada (99% VSEARCH), secuencias y taxonomía.
# USO: bash 0_preprocessing/scripts/06b_visualize_unified_table.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."
UNIFIED_DIR="$ROOT_DIR/data/pre_processed/qiime2_workspace/unified"

echo "=========================================="
echo "Verificando archivos unificados (.qza)..."
echo "=========================================="

if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno con: conda activate qiime2-amplicon-2024.5"
    exit 1
fi

# Verificar artefactos de entrada obligatorios
REQUIRED_FILES=(
    "$UNIFIED_DIR/final_unified_table.qza"
    "$UNIFIED_DIR/unified_rep_seqs.qza"
    "$UNIFIED_DIR/taxonomy.qza"
)

for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "❌ ERROR: No se encontró el archivo requerido: $file"
        echo "   Ejecuta primero: bash 0_preprocessing/scripts/06_merge_cluster_export.sh"
        exit 1
    fi
done

echo "=========================================="
echo "📊 Generando visualizaciones unificadas (.qzv)..."
echo "=========================================="

# 1. Resumen de la Tabla de Abundancia Unificada
echo "📈 [1/3] Generando final_unified_table.qzv..."
qiime feature-table summarize \
    --i-table "$UNIFIED_DIR/final_unified_table.qza" \
    --o-visualization "$UNIFIED_DIR/final_unified_table.qzv"

# 2. Resumen de Secuencias Representativas (99% Identity)
echo "🧬 [2/3] Generando unified_rep_seqs.qzv..."
qiime metadata tabulate \
    --m-input-file "$UNIFIED_DIR/unified_rep_seqs.qza" \
    --o-visualization "$UNIFIED_DIR/unified_rep_seqs.qzv"
    
# 3. Resumen de Clasificación Taxonómica SILVA
echo "🏷️ [3/3] Generando taxonomy.qzv..."
qiime metadata tabulate \
    --m-input-file "$UNIFIED_DIR/taxonomy.qza" \
    --o-visualization "$UNIFIED_DIR/taxonomy.qzv"

echo "=========================================="
echo "🎉 ¡Visualizaciones creadas exitosamente!"
echo "📁 Ubicación: $UNIFIED_DIR"
echo "🌐 Arrastra los archivos .qzv a https://view.qiime2.org para su inspección."
echo "=========================================="