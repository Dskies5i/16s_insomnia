#!/bin/bash
# ============================================================
# NAME: 01_install_qiime2.sh
# DESCRIPTION: Installs QIIME 2 (Amplicon) on Apple Silicon macOS
# USAGE: bash 0_preprocessing/scripts/01_install_qiime2.sh
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

# ============================================================
echo "Running  02_install_qiime2.sh"
# ============================================================

set -e

# --- 1. Initialize Conda for Bash Scripting ---
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")

if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "❌ ERROR: Conda environment initialization script not found."
    exit 1
fi

ENV_NAME="qiime2-amplicon-2024.5"

# --- 2. Create QIIME 2 Environment if it doesn't exist ---
if ! conda env list | grep -q "$ENV_NAME"; then
    echo "=========================================="
    echo "Installing QIIME 2 ($ENV_NAME) for macOS..."
    echo "This may take 10-15 minutes."
    echo "=========================================="
    
    ENV_FILE="qiime2-amplicon-2024.5-py39-osx-conda.yml"
    echo "Downloading QIIME 2 environment configuration..."
    curl -L -sS "https://data.qiime2.org/distro/amplicon/$ENV_FILE" -o /tmp/qiime2.yml
    
    if [ ! -s /tmp/qiime2.yml ]; then
        echo "❌ ERROR: Failed to download $ENV_FILE"
        exit 1
    fi
    
    # Force classic solver temporarily during environment creation to avoid libmamba conflicts
    CONDA_SOLVER=classic CONDA_SUBDIR=osx-64 conda env create -n "$ENV_NAME" --file /tmp/qiime2.yml
    rm -f /tmp/qiime2.yml
else
    echo "✅ QIIME 2 environment '$ENV_NAME' already exists."
fi

# --- 3. Activate Environment and Verify ---
echo "=========================================="
echo "Activating $ENV_NAME..."
echo "=========================================="
conda activate "$ENV_NAME"

echo "Checking QIIME 2 version:"
qiime --version

echo "=========================================="
echo "🎉 QIIME 2 is installed and ready to use!"
echo "Note: Before running future steps in your terminal, run:"
echo "      conda activate $ENV_NAME"
echo "=========================================="

Pro Tip: Once your script completes and prints 🎉 QIIME 2 is installed and ready to use!, you can free up several gigabytes of cached installer files by running:

# Only deletes the leftover compressed installer archives and temporary index files
conda clean --all -y



## ---------------------------------------------------
#### TROUBLESHOOT
## ---------------------------------------------------
####1.- WGET NOT WORKING:
## Restore base conda subdir configuration to native Apple Silicon
# conda config --env --remove-key subdir 2>/dev/null || true
# CONDA_SUBDIR=osx-arm64 conda install -n base --force-reinstall openssl ca-certificates certifi conda -y
# CONDA_SOLVER=classic CONDA_SUBDIR=osx-arm64 conda install -n base --force-reinstall openssl ca-certificates certifi conda -y
# Verify: conda info
## --> did not work

## ---------------------------------------------------
#### 2.- Push comes to shovel --- reinstall conda:
# brew reinstall --cask miniconda
## -- Innitialize shell profile 
# conda init zsh
# source ~/.zshrc

## -- Verify repair:
# conda info
# conda list openssl


## ---------------------------------------------------

