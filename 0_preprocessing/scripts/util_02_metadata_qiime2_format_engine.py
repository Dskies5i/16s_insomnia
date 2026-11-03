

#!/usr/bin/env python3

#util_02_metadata_qiime2_format_engine.py
import sys
import pandas as pd

def generate_qiime2_metadata(metadata_csv: str, output_tsv: str):
    df = pd.read_csv(metadata_csv)
    
    # 1. CRITICAL: Primary sample-id MUST be run_accession (SRR...)
    if "run_accession" not in df.columns:
        raise KeyError("❌ ERROR: 'run_accession' column missing from input CSV!")
    
    # Extract study_id from filename (e.g., haimov_2022)
    study_name = output_tsv.split("/")[-1].replace("_metadata_q2.tsv", "")
    df["study_id"] = study_name

    # Select standard columns present in SRA metadata
    desired_cols = [
        "run_accession", 
        "study_id", 
        "library_name", 
        "sample_accession", 
        "sex", 
        "age", 
        "tissue", 
        "biomaterial_provider"
    ]
    
    selected_cols = [col for col in desired_cols if col in df.columns]
    metadata = df[selected_cols].copy()
    
    # Rename run_accession to sample-id for QIIME 2 compatibility
    metadata.rename(columns={"run_accession": "sample-id"}, inplace=True)
    metadata.drop_duplicates(subset=["sample-id"], inplace=True)
    
    # Add QIIME 2 type header row (#q2:types)
    type_row = {}
    for col in metadata.columns:
        if col == "sample-id":
            type_row[col] = "#q2:types"
        elif col in ["age"]:
            type_row[col] = "numeric"
        else:
            type_row[col] = "categorical"
            
    types_df = pd.DataFrame([type_row])
    final_metadata = pd.concat([types_df, metadata], ignore_index=True)
    
    # Format and save as TSV
    final_metadata.to_csv(output_tsv, sep="\t", index=False)
    print(f"✅ Saved QIIME 2 metadata: {output_tsv}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 02_prepare_qiime2_metadata.py <input_csv> <output_tsv>")
        sys.exit(1)
    generate_qiime2_metadata(sys.argv[1], sys.argv[2])
