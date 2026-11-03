#!/bin/bash

# ============================================================
# NAME: 00_install_sra_tools.sh
# DESCRIPTION: Downloads and configures SRA Toolkit on macOS permanently
#              if it is not already installed.
# USAGE: bash 0_preprocessing/scripts/00_install_sra_tools.sh
# ENV: conda activate qiime2-amplicon-2024.5
# 
# ============================================================

set -e

# Force PATH inclusion so Conda subshells always see the SRA binaries
export PATH="/usr/local/sra-toolkit/bin:$PATH"

echo "=========================================="
echo "Step 1: Checking for existing installation..."
echo "=========================================="

if command -v prefetch &> /dev/null && command -v fasterq-dump &> /dev/null; then
    echo "✅ SRA Toolkit is already installed!"
    echo "📍 prefetch location: $(which prefetch)"
    echo "📍 version: $(prefetch --version 2>&1 | head -n 1)"
    echo "⏭️ Skipping installation."
    exit 0
fi

echo "⚠️ SRA Toolkit not detected. Proceeding with installation..."

echo "=========================================="
echo "Step 2: Downloading SRA Toolkit for macOS..."
echo "=========================================="
cd ~/Downloads
curl -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz

echo "=========================================="
echo "Step 3: Extracting archive..."
echo "=========================================="
tar -xzf sratoolkit.current-mac64.tar.gz
EXTRACTED_DIR=$(tar -tf sratoolkit.current-mac64.tar.gz | head -n 1 | cut -f1 -d"/")

echo "=========================================="
echo "Step 4: Moving binaries to /usr/local/sra-toolkit..."
echo "=========================================="
sudo rm -rf /usr/local/sra-toolkit
sudo mv "$EXTRACTED_DIR" /usr/local/sra-toolkit
rm -f sratoolkit.current-mac64.tar.gz

echo "=========================================="
echo "Step 5: Updating shell profile PATH..."
echo "=========================================="
# Clean up any legacy PATH lines
sed -i '' '/sra-toolkit/d' ~/.zshrc 2>/dev/null || true
sed -i '' '/sra-toolkit/d' ~/.zprofile 2>/dev/null || true

# Append permanent path export
echo 'export PATH="/usr/local/sra-toolkit/bin:$PATH"' >> ~/.zshrc
echo 'export PATH="/usr/local/sra-toolkit/bin:$PATH"' >> ~/.zprofile

echo "=========================================="
echo "Step 6: Verifying installation..."
echo "=========================================="
if command -v prefetch &> /dev/null && command -v fasterq-dump &> /dev/null; then
    echo "🎉 SRA Toolkit installed successfully!"
    echo "📍 prefetch location: $(which prefetch)"
    echo "📍 fasterq-dump location: $(which fasterq-dump)"
else
    echo "❌ Installation failed or binaries are not executable."
    exit 1
fi

#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#### $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




# ## ------- Installing SRA toolkit MACBOOK
# 
# ## 00_install_sra_tools.sh
# echo "1.- Installing SRA toolkit MACBOOK"
# cd ~/Downloads
# curl -O https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
# tar -vxzf sratoolkit.current-mac64.tar.gz
# 
# # Move it to a permanent location
# echo "2.- Move it to a permanent location"
# #cd ~/Downloads
# sudo mv ~/Downloads/sratoolkit.3.4.1-mac-x86_64 /usr/local/sra-toolkit
# 
# # Remove the old Downloads path from .zshrc
# echo "3.- Remove the old Downloads path from .zshrc"
# sed -i '' '/Downloads\/sratoolkit/d' ~/.zshrc
# 
# # Add the new permanent path
# echo "4.- Add the new permanent path"
# echo 'export PATH="/usr/local/sra-toolkit/bin:$PATH"' >> ~/.zprofile
# 
# # Reload the configuration
# echo "5.- Reload the configuration"
# source ~/.zprofile
# 
# #Verify
# echo "6.- verify"
# which fasterq-dump
# 
# #expected result: "/usr/local/sra-toolkit/bin/prefetch"
# 
