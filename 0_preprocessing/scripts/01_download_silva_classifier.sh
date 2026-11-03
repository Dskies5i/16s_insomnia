#!/bin/bash

# ============================================================
# NAME: 01_download_silva_classifier.sh
# DESCRIPTION: Downloads pre-trained Full-Length SILVA taxonomic classifier for QIIME 2 2024.5
# USAGE: bash 0_preprocessing/scripts/01_download_silva_classifier.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

echo "Running 03_download_silva_classifier.sh"

set -e

PROJECT_ROOT=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
TARGET_DIR="$PROJECT_ROOT/resources/classifiers"

mkdir -p "$TARGET_DIR"

# Correct official URL and filename for QIIME 2 2024.5 (scikit-learn 1.4.2)
CLASSIFIER_NAME="silva-138-99-nb-classifier.qza"
URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
OUTPUT_FILE="$TARGET_DIR/$CLASSIFIER_NAME"

echo "=========================================="
echo "Target Directory: $TARGET_DIR"
echo "=========================================="

if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    echo "✅ Full-length classifier already exists: $OUTPUT_FILE"
else
    echo "📥 Downloading SILVA 138 Full-Length classifier (~500MB)..."
    
    if command -v wget &> /dev/null; then
        wget -O "$OUTPUT_FILE" "$URL"
    elif command -v curl &> /dev/null; then
        curl -L -o "$OUTPUT_FILE" "$URL"
    else
        echo "❌ ERROR: Neither wget nor curl is installed."
        exit 1
    fi

    if [ ! -s "$OUTPUT_FILE" ]; then
        echo "❌ ERROR: Download failed or resulted in an empty file."
        rm -f "$OUTPUT_FILE"
        exit 1
    fi
    echo "✅ Download complete!"
fi

echo "=========================================="
echo "🎉 Full-length classifier ready at:"
echo "   $OUTPUT_FILE"
echo "=========================================="

# In case it needs permission
# chmod +x 0_preprocessing/scripts/03_download_silva_classifier.sh

# RUN THE SCRIPT:
# bash 0_preprocessing/scripts/03_download_silva_classifier.sh