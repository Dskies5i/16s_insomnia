#!/bin/bash

# ============================================================
# NOMBRE: 06_merge_cluster_export.sh
# DESCRIPCIÓN: Unificación VSEARCH 99%, clasificación SILVA y exportación a R
#              filtrando dinámicamente datasets activos desde datasets_config.tsv
#
# ¿QUÉ HACE ESTE SCRIPT (META-ANÁLISIS)?
# En lugar de procesar cada estudio por separado y fusionar tablas incompatibles 
# al final, este script realiza un flujo unificado y robusto:
#   1. Fusión de Tablas y Secuencias: Combina las tablas de conteo (table.qza) 
#      y secuencias representativas (rep_seqs.qza) de todos los estudios activos 
#      en datasets_config.tsv.
#   2. Clustering Global VSEARCH (99%): Agrupa todas las secuencias mediante 
#      de novo clustering al 99% de identidad. Esto genera un diccionario de 
#      características (OTUs/ASVs) único y compartido para todos los datasets.
#   3. Clasificación Taxonómica SILVA Unificada: Ejecuta el clasificador SILVA 
#      una sola vez sobre el archivo maestro de secuencias unificadas, produciendo 
#      un único taxonomy.qza para todas las muestras.
#   4. Exportación para R: Exporta las matrices y archivos FASTA limpios a 
#      data/pre_processed/unified/ para cargarlos directamente en R (phyloseq).
#
# USO: bash 0_preprocessing/scripts/06_merge_cluster_export.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================
# ------------------------------------------------------------
# ------------------------------------------------------------
# NOTA DE ENTORNO (macOS ARM / Conda):
# Si DADA2/Rcpp falló previamente en 05, asegurar compilación correcta:
# Rscript -e "install.packages('Rcpp', repos='https://cloud.r-project.org')"
# ------------------------------------------------------------
# ------------------------------------------------------------

set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
ROOT_DIR="$SCRIPT_DIR/../.."

CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
QIIME2_DIR="$ROOT_DIR/data/pre_processed/qiime2_workspace"
    # ******* qiime2_workspace/unified IS THE INTERNAL WORKING directory where QIIME 2 processes and modifies the data.
    # merged_table_raw.qza & merged_rep_seqs.qza (merged raw tables)
    # final_unified_table.qza & unified_rep_seqs.qza (after VSEARCH clustering)
    # taxonomy.qza (SILVA taxonomic classifications)
UNIFIED_DIR="$QIIME2_DIR/unified"
    # ******** data/processed/unified -- THIS IS THE HUMAN READABLE DIR containing clean text tables and FASTA files for downstream statistical analysis in R.
    # unified_otu_table.tsv and unified_otu_table.csv (abundance tables)
    # unified_rep_seqs.fasta (representative sequences)
    # Taxonomy mapping files (taxonomy.tsv)
EXPORT_DIR="$ROOT_DIR/data/pre_processed/unified"
SILVA_CLASSIFIER="$ROOT_DIR/resources/classifiers/silva-138-99-nb-classifier.qza"
THREADS=4 # Reducido a 4 para estabilidad de memoria RAM en macOS ARM

echo "=========================================="
echo "Verificando entorno y archivos de entrada..."
echo "=========================================="

if ! command -v qiime &> /dev/null; then
    echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno con: conda activate qiime2-amplicon-2024.5"
    exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ ERROR: Archivo de configuración no encontrado en: $CONFIG_FILE"
    exit 1
fi

if [ ! -f "$SILVA_CLASSIFIER" ]; then
    echo "❌ ERROR: Clasificador SILVA no encontrado en: $SILVA_CLASSIFIER"
    echo "   Ejecuta primero: bash 0_preprocessing/scripts/03_download_silva_classifier.sh"
    exit 1
fi

# Sanitización en sitio del archivo de configuración (remueve \r de Windows)
tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"

mkdir -p "$UNIFIED_DIR" "$EXPORT_DIR"

TABLE_ARGS=()
SEQ_ARGS=()

echo "🔍 Buscando salidas procesadas para estudios activos en datasets_config.tsv..."

# 1. Validación de existencia de archivos (Uso de process substitution para evitar subshell)
while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do
    
    STUDY=$(echo "$STUDY" | xargs)
    [[ -z "$STUDY" || "$STUDY" =~ ^# ]] && continue

    STUDY_DIR="$QIIME2_DIR/$STUDY"
    
    if [ -f "$STUDY_DIR/table.qza" ] && [ -f "$STUDY_DIR/rep_seqs.qza" ]; then
        echo "   ✔ Registrado: $STUDY"
        TABLE_ARGS+=(--i-tables "$STUDY_DIR/table.qza")
        SEQ_ARGS+=(--i-data "$STUDY_DIR/rep_seqs.qza")
    else
        echo "   ❌ ERROR: Falta table.qza o rep_seqs.qza para $STUDY en $STUDY_DIR."
        echo "      Ejecuta primero: bash 0_preprocessing/scripts/05_dada2_denoise_datasets.sh"
        exit 1
    fi
done < <(tail -n +2 "$CONFIG_FILE")

if [ ${#TABLE_ARGS[@]} -eq 0 ]; then
    echo "❌ ERROR: No se encontraron tablas válidas para fusionar."
    exit 1
fi

# 2. Merge
if [ -f "$UNIFIED_DIR/merged_table_raw.qza" ] && [ -f "$UNIFIED_DIR/merged_rep_seqs.qza" ]; then
    echo "=========================================="
    echo "⏭️ El merge ya existe. Saltando paso..."
    echo "=========================================="
else
    echo "=========================================="
    echo "🔀 Fusionando tablas de abundancia y secuencias..."
    echo "=========================================="
    qiime feature-table merge "${TABLE_ARGS[@]}" --o-merged-table "$UNIFIED_DIR/merged_table_raw.qza"
    qiime feature-table merge-seqs "${SEQ_ARGS[@]}" --o-merged-data "$UNIFIED_DIR/merged_rep_seqs.qza"
fi

# 3. VSEARCH 99%
if [ -f "$UNIFIED_DIR/final_unified_table.qza" ] && [ -f "$UNIFIED_DIR/unified_rep_seqs.qza" ]; then
    echo "=========================================="
    echo "⏭️ El clustering por VSEARCH ya existe. Saltando paso..."
    echo "=========================================="
else
    echo "=========================================="
    echo "🔬 Clustering de novo con VSEARCH al 99% de identidad..."
    echo "=========================================="
    qiime vsearch cluster-features-de-novo \
        --i-sequences "$UNIFIED_DIR/merged_rep_seqs.qza" \
        --i-table "$UNIFIED_DIR/merged_table_raw.qza" \
        --p-perc-identity 0.99 \
        --p-threads "$THREADS" \
        --o-clustered-table "$UNIFIED_DIR/final_unified_table.qza" \
        --o-clustered-sequences "$UNIFIED_DIR/unified_rep_seqs.qza"
fi

# 4. SILVA Taxonomy Classification
if [ -f "$UNIFIED_DIR/taxonomy.qza" ]; then
    echo "=========================================="
    echo "⏭️ La taxonomía SILVA ya existe. Saltando paso..."
    echo "=========================================="
else
    echo "=========================================="
    echo "🏷️ Asignando taxonomía con el clasificador SILVA..."
    echo "=========================================="
    qiime feature-classifier classify-sklearn \
        --i-classifier "$SILVA_CLASSIFIER" \
        --i-reads "$UNIFIED_DIR/unified_rep_seqs.qza" \
        --p-n-jobs "$THREADS" \
        --o-classification "$UNIFIED_DIR/taxonomy.qza"
fi

# 5. Exports for R
if [ -f "$EXPORT_DIR/unified_otu_table.tsv" ] && [ -f "$EXPORT_DIR/unified_otu_table.csv" ] && [ -f "$EXPORT_DIR/unified_rep_seqs.fasta" ] && [ -f "$EXPORT_DIR/taxonomy.tsv" ]; then
    echo "=========================================="
    echo "⏭️ Los archivos exportados para R ya existen. Saltando paso..."
    echo "=========================================="
else
    echo "=========================================="
    echo "📦 Exportando matrices y archivos FASTA para R..."
    echo "=========================================="
    qiime tools export --input-path "$UNIFIED_DIR/final_unified_table.qza" --output-path "$EXPORT_DIR"

    # Convert BIOM to clean TSV/CSV without BIOM comment headers
    biom convert -i "$EXPORT_DIR/feature-table.biom" -o "$EXPORT_DIR/unified_otu_table_raw.tsv" --to-tsv
    grep -v '^# Construct' "$EXPORT_DIR/unified_otu_table_raw.tsv" | sed 's/^#OTU ID/OTU_ID/' > "$EXPORT_DIR/unified_otu_table.tsv"
    sed 's/\t/,/g' "$EXPORT_DIR/unified_otu_table.tsv" > "$EXPORT_DIR/unified_otu_table.csv"
    rm "$EXPORT_DIR/unified_otu_table_raw.tsv"

    qiime tools export --input-path "$UNIFIED_DIR/unified_rep_seqs.qza" --output-path "$EXPORT_DIR"
    [ -f "$EXPORT_DIR/dna-sequences.fasta" ] && mv "$EXPORT_DIR/dna-sequences.fasta" "$EXPORT_DIR/unified_rep_seqs.fasta"

    qiime tools export --input-path "$UNIFIED_DIR/taxonomy.qza" --output-path "$EXPORT_DIR"
fi

echo "=========================================="
echo "🎉 Proceso finalizado con éxito en: $EXPORT_DIR"
echo "=========================================="



# #===============================================================================
# #===============================================================================
# #===============================================================================
# #!/bin/bash
# 
# # ============================================================
# # NOMBRE: 06_merge_cluster_export.sh
# # DESCRIPCIÓN: Unificación VSEARCH 99%, clasificación SILVA y exportación a R
# #              filtrando dinámicamente datasets activos desde datasets_config.tsv
# # USO: bash 0_preprocessing/scripts/06_merge_cluster_export.sh
# #
# # NOTA DE ENTORNO (macOS ARM / Conda):
# # Si DADA2/Rcpp falló previamente en 05, asegurar compilación correcta:
# # Rscript -e "install.packages('Rcpp', repos='https://cloud.r-project.org')"
# # ============================================================
# 
# set -e
# 
# SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# ROOT_DIR="$SCRIPT_DIR/../.."
# 
# CONFIG_FILE="$SCRIPT_DIR/../config/datasets_config.tsv"
# QIIME2_DIR="$ROOT_DIR/data/pre_processed/qiime2_workspace"
# UNIFIED_DIR="$QIIME2_DIR/unified"
# EXPORT_DIR="$ROOT_DIR/data/pre_processed/unified"
# SILVA_CLASSIFIER="$ROOT_DIR/resources/classifiers/silva-138-99-nb-classifier.qza"
# THREADS=4 # Reducido a 4 para estabilidad de memoria RAM en macOS ARM
# 
# echo "=========================================="
# echo "Verificando entorno y archivos de entrada..."
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
# if [ ! -f "$SILVA_CLASSIFIER" ]; then
#     echo "❌ ERROR: Clasificador SILVA no encontrado en: $SILVA_CLASSIFIER"
#     echo "   Ejecuta primero: bash 0_preprocessing/scripts/03_download_silva_classifier.sh"
#     exit 1
# fi
# 
# # Sanitización en sitio del archivo de configuración (remueve \r de Windows)
# tr -d '\r' < "$CONFIG_FILE" > "${CONFIG_FILE}.tmp" && mv "${CONFIG_FILE}.tmp" "$CONFIG_FILE"
# 
# mkdir -p "$UNIFIED_DIR" "$EXPORT_DIR"
# 
# TABLE_ARGS=()
# SEQ_ARGS=()
# 
# echo "🔍 Buscando salidas procesadas para estudios activos en datasets_config.tsv..."
# 
# # 1. Validación de existencia de archivos (Uso de process substitution para evitar subshell)
# while IFS=$'\t' read -r STUDY FWD REV TF TR LF LR REGION || [ -n "$STUDY" ]; do
#     
#     STUDY=$(echo "$STUDY" | xargs)
#     [[ -z "$STUDY" || "$STUDY" =~ ^# ]] && continue
# 
#     STUDY_DIR="$QIIME2_DIR/$STUDY"
#     
#     if [ -f "$STUDY_DIR/table.qza" ] && [ -f "$STUDY_DIR/rep_seqs.qza" ]; then
#         echo "   ✔ Registrado: $STUDY"
#         TABLE_ARGS+=(--i-tables "$STUDY_DIR/table.qza")
#         SEQ_ARGS+=(--i-data "$STUDY_DIR/rep_seqs.qza")
#     else
#         echo "   ❌ ERROR: Falta table.qza o rep_seqs.qza para $STUDY en $STUDY_DIR."
#         echo "      Ejecuta primero: bash 0_preprocessing/scripts/05_dada2_denoise_datasets.sh"
#         exit 1
#     fi
# done < <(tail -n +2 "$CONFIG_FILE")
# 
# if [ ${#TABLE_ARGS[@]} -eq 0 ]; then
#     echo "❌ ERROR: No se encontraron tablas válidas para fusionar."
#     exit 1
# fi
# 
# # 2. Merge
# echo "=========================================="
# echo "🔀 Fusionando tablas de abundancia y secuencias..."
# echo "=========================================="
# 
# qiime feature-table merge "${TABLE_ARGS[@]}" --o-merged-table "$UNIFIED_DIR/merged_table_raw.qza"
# qiime feature-table merge-seqs "${SEQ_ARGS[@]}" --o-merged-data "$UNIFIED_DIR/merged_rep_seqs.qza"
# 
# # 3. VSEARCH 99%
# echo "=========================================="
# echo "🔬 Clustering de novo con VSEARCH al 99% de identidad..."
# echo "=========================================="
# qiime vsearch cluster-features-de-novo \
#     --i-sequences "$UNIFIED_DIR/merged_rep_seqs.qza" \
#     --i-table "$UNIFIED_DIR/merged_table_raw.qza" \
#     --p-perc-identity 0.99 \
#     --p-threads "$THREADS" \
#     --o-clustered-table "$UNIFIED_DIR/final_unified_table.qza" \
#     --o-clustered-sequences "$UNIFIED_DIR/unified_rep_seqs.qza"
# 
# # 4. SILVA Taxonomy Classification
# echo "=========================================="
# echo "🏷️ Asignando taxonomía con el clasificador SILVA..."
# echo "=========================================="
# qiime feature-classifier classify-sklearn \
#     --i-classifier "$SILVA_CLASSIFIER" \
#     --i-reads "$UNIFIED_DIR/unified_rep_seqs.qza" \
#     --p-n-jobs "$THREADS" \
#     --o-classification "$UNIFIED_DIR/taxonomy.qza"
# 
# # 5. Exports for R
# echo "=========================================="
# echo "📦 Exportando matrices y archivos FASTA para R..."
# echo "=========================================="
# qiime tools export --input-path "$UNIFIED_DIR/final_unified_table.qza" --output-path "$EXPORT_DIR"
# 
# # Convert BIOM to clean TSV/CSV without BIOM comment headers
# biom convert -i "$EXPORT_DIR/feature-table.biom" -o "$EXPORT_DIR/unified_otu_table_raw.tsv" --to-tsv
# grep -v '^# Construct' "$EXPORT_DIR/unified_otu_table_raw.tsv" | sed 's/^#OTU ID/OTU_ID/' > "$EXPORT_DIR/unified_otu_table.tsv"
# sed 's/\t/,/g' "$EXPORT_DIR/unified_otu_table.tsv" > "$EXPORT_DIR/unified_otu_table.csv"
# rm "$EXPORT_DIR/unified_otu_table_raw.tsv"
# 
# qiime tools export --input-path "$UNIFIED_DIR/unified_rep_seqs.qza" --output-path "$EXPORT_DIR"
# [ -f "$EXPORT_DIR/dna-sequences.fasta" ] && mv "$EXPORT_DIR/dna-sequences.fasta" "$EXPORT_DIR/unified_rep_seqs.fasta"
# 
# qiime tools export --input-path "$UNIFIED_DIR/taxonomy.qza" --output-path "$EXPORT_DIR"
# 
# echo "=========================================="
# echo "🎉 Exportación completada con éxito en: $EXPORT_DIR"
# echo "=========================================="



# # ============================================================
# # ============================================================
# # = OLDER VERSION
# # ============================================================
# # ============================================================
# #!/bin/bash
# 
# # ============================================================
# # NOMBRE: 06_merge_cluster_export.sh
# # DESCRIPCIÓN: Unificación VSEARCH 99%, clasificación SILVA y exportación a R
# # USO: bash 0_preprocessing/scripts/06_merge_cluster_export.sh
# # ============================================================
# 
# set -e
# 
# SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
# ROOT_DIR="$SCRIPT_DIR/../.."
# 
# QIIME2_DIR="$ROOT_DIR/data/processed/qiime2_workspace"
# UNIFIED_DIR="$QIIME2_DIR/unified"
# EXPORT_DIR="$ROOT_DIR/data/processed/unified"
# SILVA_CLASSIFIER="$ROOT_DIR/resources/classifiers/silva-138-99-nb-classifier.qza"
# THREADS=8
# 
# echo "=========================================="
# echo "Verificando entorno y archivos de entrada..."
# echo "=========================================="
# 
# if ! command -v qiime &> /dev/null; then
#     echo "❌ ERROR: QIIME 2 no detectado. Activa el entorno con: conda activate qiime2-amplicon-2024.5"
#     exit 1
# fi
# 
# if [ ! -f "$SILVA_CLASSIFIER" ]; then
#     echo "❌ ERROR: Clasificador SILVA no encontrado en: $SILVA_CLASSIFIER"
#     echo "   Ejecuta primero: bash 0_preprocessing/scripts/03_download_silva_classifier.sh"
#     exit 1
# fi
# 
# mkdir -p "$UNIFIED_DIR" "$EXPORT_DIR"
# 
# ## ---------------------------------------------------
# ## Dynamically gather all processed dataset tables and rep_seqs
# ## ---------------------------------------------------
# 
# TABLE_ARGS=()
# SEQ_ARGS=()
# 
# for study_dir in "$QIIME2_DIR"/*/; do
#     study_name=$(basename "$study_dir")
#     
#     # Exclude the unified target directory itself
#     if [ "$study_name" != "unified" ] && [ "$study_name" != "manifests" ]; then
#         if [ -f "$study_dir/table.qza" ] && [ -f "$study_dir/rep_seqs.qza" ]; then
#             TABLE_ARGS+=(--i-tables "$study_dir/table.qza")
#             SEQ_ARGS+=(--i-sequences "$study_dir/rep_seqs.qza")
#             echo " Found processed data for: $study_name"
#         fi
#     fi
# done
# ## ---------------------------------------------------
# # 1. Merge
# ## ---------------------------------------------------
# qiime feature-table merge "${TABLE_ARGS[@]}" --o-merged-table "$UNIFIED_DIR/merged_table_raw.qza"
# qiime feature-table merge-seqs "${SEQ_ARGS[@]}" --o-merged-sequences "$UNIFIED_DIR/merged_rep_seqs.qza"
# 
# ## ---------------------------------------------------
# ## 2. VSEARCH 99%
# ## ---------------------------------------------------
# qiime vsearch cluster-features-de-novo \
#     --i-sequences "$UNIFIED_DIR/merged_rep_seqs.qza" \
#     --i-table "$UNIFIED_DIR/merged_table_raw.qza" \
#     --p-perc-identity 0.99 \
#     --p-threads "$THREADS" \
#     --o-clustered-table "$UNIFIED_DIR/final_unified_table.qza" \
#     --o-clustered-sequences "$UNIFIED_DIR/unified_rep_seqs.qza"
# 
# ## ---------------------------------------------------
# # 3. SILVA Taxonomy Classification
# ## ---------------------------------------------------
# qiime feature-classifier classify-sklearn \
#     --i-classifier "$SILVA_CLASSIFIER" \
#     --i-reads "$UNIFIED_DIR/unified_rep_seqs.qza" \
#     --p-n-jobs "$THREADS" \
#     --o-classification "$UNIFIED_DIR/taxonomy.qza"
#     
# ## ---------------------------------------------------
# # 4. Exports for R
# ## ---------------------------------------------------
# qiime tools export --input-path "$UNIFIED_DIR/final_unified_table.qza" --output-path "$EXPORT_DIR"
# biom convert -i "$EXPORT_DIR/feature-table.biom" -o "$EXPORT_DIR/unified_otu_table.tsv" --to-tsv
# sed 's/\t/,/g' "$EXPORT_DIR/unified_otu_table.tsv" > "$EXPORT_DIR/unified_otu_table.csv"
# 
# qiime tools export --input-path "$UNIFIED_DIR/unified_rep_seqs.qza" --output-path "$EXPORT_DIR"
# [ -f "$EXPORT_DIR/dna-sequences.fasta" ] && mv "$EXPORT_DIR/dna-sequences.fasta" "$EXPORT_DIR/unified_rep_seqs.fasta"
# 
# qiime tools export --input-path "$UNIFIED_DIR/taxonomy.qza" --output-path "$EXPORT_DIR"