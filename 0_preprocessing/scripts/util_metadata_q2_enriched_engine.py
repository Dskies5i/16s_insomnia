#!/usr/bin/env python3

#util_metadata_q2_enriched_engine.py

"""
utils/metadata_q2_enriched_engine.py

Core processing engine to standardise, enrich, and format SRA metadata into 
QIIME 2 compatible metadata files (#q2:types) for microbiome meta-analyses.
"""

import os
import pandas as pd


def parse_prefix_group(series: pd.Series, prefix_map: dict, default: str = "Unknown") -> pd.Series:
    """Map string prefixes in SRA identifiers to group categories.
    
    Example prefix_map: {"A": "Insomnia", "C": "Control"}
    """
    s = series.astype(str).str.strip().str.upper()
    
    def match_prefix(val):
        for prefix, group_label in prefix_map.items():
            if val.startswith(prefix.upper()):
                return group_label
        return default
        
    return s.apply(match_prefix)


def standardize_series(series: pd.Series, custom_regex: str = None) -> pd.Series:
    """Generic ID cleaner with optional regex override for joining datasets."""
    s = series.astype(str).str.replace(r"^X", "", regex=True).str.replace(r"\.", "-", regex=True)
    if custom_regex:
        return s.str.extract(custom_regex)[0]
    return s.str.extract(r"^([A-Za-z0-9]+(?:-[0-9]+)?)")[0]


def enrich_study_metadata(
    sra_path: str,
    supp_path: str = None,
    output_path: str = None,
    id_col_sra: str = "library_name",
    id_col_supp: str = None,
    custom_regex: str = None,
    clinical_cols: list = None,
    numeric_cols: list = None,
    group_prefix_map: dict = None,
    value_mappings: dict = None
):
    """
    Standardises SRA metadata, optionally merges external supplementary clinical data,
    harmonises phenotypic variables (including multi-class insomnia subtypes), and formats
    the final table with QIIME 2 directive headers.
    """
    if not os.path.exists(sra_path):
        print(f"⚠️ Skipping: missing SRA metadata at {sra_path}")
        return

    sra_df = pd.read_csv(sra_path, sep="\t", comment="#")

    # -------------------------------------------------------------------------
    # 1. LOAD & MERGE DATA
    # -------------------------------------------------------------------------
    # DIRECT MODE: No supplementary table required (e.g., Chen 2022)
    if supp_path is None or not os.path.exists(supp_path):
        print(f"ℹ️ Processing {os.path.basename(sra_path)} directly from SRA attributes...")
        merged_df = sra_df.copy()
        
        if group_prefix_map:
            merged_df["group"] = parse_prefix_group(merged_df[id_col_sra], group_prefix_map)

    # JOIN MODE: External supplementary table available (e.g., Haimov 2022)
    else:
        print(f"ℹ️ Merging {os.path.basename(sra_path)} with supplementary file {os.path.basename(supp_path)}...")
        supp_df = pd.read_csv(supp_path)
        supp_df.columns = supp_df.columns.str.strip()
        supp_id = id_col_supp or ("SampleID" if "SampleID" in supp_df.columns else supp_df.columns[0])

        sra_df["join_key"] = standardize_series(sra_df[id_col_sra], custom_regex)
        supp_df["join_key"] = standardize_series(supp_df[supp_id], custom_regex)

        supp_unique = supp_df.drop_duplicates(subset=["join_key"]).copy()

        if clinical_cols:
            target_cols = ["join_key"] + [c for c in clinical_cols if c in supp_unique.columns]
            supp_unique = supp_unique[target_cols]

        merged_df = pd.merge(sra_df, supp_unique, on="join_key", how="left")
        merged_df.drop(columns=["join_key"], inplace=True, errors="ignore")

    # -------------------------------------------------------------------------
    # 2. STANDARDISE VALUES, COLUMNS & MACRO-GROUPS
    # -------------------------------------------------------------------------
    if value_mappings:
        # Step A: Replace specific categorical values (e.g., "feces" -> "stool", "Ins" -> "Insomnia")
        for col, col_map in value_mappings.get("replace_values", {}).items():
            if col in merged_df.columns:
                merged_df[col] = merged_df[col].replace(col_map)
        
        # Step B: Rename columns to standardised naming (e.g., "diagnosis" -> "group_detailed")
        if "rename_cols" in value_mappings:
            merged_df.rename(columns=value_mappings["rename_cols"], inplace=True)

        # Step C: Derive primary macro 'group' column (e.g., collapse Paradoxical Insomnia & Insomnia -> Insomnia)
        if "derive_macro_group" in value_mappings:
            cfg = value_mappings["derive_macro_group"]
            src_col = cfg.get("source_col")
            mapping = cfg.get("mapping", {})
            if src_col and src_col in merged_df.columns:
                merged_df["group"] = merged_df[src_col].map(mapping).fillna("Unknown")

    # -------------------------------------------------------------------------
    # 3. GUARANTEE DUAL-COLUMN HARMONISATION ACROSS ALL STUDIES
    # -------------------------------------------------------------------------
    # If a study only has a binary 'group', populate 'group_detailed' with the same values.
    if "group" in merged_df.columns and "group_detailed" not in merged_df.columns:
        merged_df["group_detailed"] = merged_df["group"]
        
    # If a study defined 'group_detailed' but not 'group', set default 'group' as 'group_detailed'
    elif "group_detailed" in merged_df.columns and "group" not in merged_df.columns:
        merged_df["group"] = merged_df["group_detailed"]

    # -------------------------------------------------------------------------
    # 4. FORMAT QIIME 2 DIRECTIVE HEADER (#q2:types) AND SAVE
    # -------------------------------------------------------------------------
    numeric_set = set(numeric_cols or ["age", "BMI", "ISI", "SE", "GDS_SUM", "SAST_SUM"])
    type_row = {}
    
    for col in merged_df.columns:
        if col in ["sample-id", "run_accession"]:
            type_row[col] = "#q2:types"
        elif col in numeric_set or pd.api.types.is_numeric_dtype(merged_df[col]):
            type_row[col] = "numeric"
        else:
            type_row[col] = "categorical"

    final_df = pd.concat([pd.DataFrame([type_row]), merged_df], ignore_index=True)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    final_df.to_csv(output_path, sep="\t", index=False)
    print(f"✅ Saved QIIME 2 metadata: {os.path.basename(output_path)}")










# ===========================================================
# OLDER VERSION (works wonders but could not detect if the original metadata file already has all the identifiers needed)
#
# #!/usr/bin/env python3
# import os
# import sys
# import pandas as pd
# 
# def standardize_series(series, custom_regex=None):
#     """Generic ID cleaner with optional regex override."""
#     s = series.astype(str).str.replace(r"^X", "", regex=True).str.replace(r"\.", "-", regex=True)
#     if custom_regex:
#         return s.str.extract(custom_regex)[0]
#     return s.str.extract(r"^([A-Za-z0-9]+(?:-[0-9]+)?)")[0]
# 
# def enrich_study_metadata(
#     sra_path: str,
#     supp_path: str,
#     output_path: str,
#     id_col_sra: str = "library_name",
#     id_col_supp: str = None,
#     custom_regex: str = None,
#     clinical_cols: list = None,
#     numeric_cols: list = None
# ):
#     if not os.path.exists(sra_path) or not os.path.exists(supp_path):
#         print(f"⚠️ Skipping enrichment: missing input files for {output_path}")
#         return
# 
#     # Load data
#     sra_df = pd.read_csv(sra_path, sep="\t", comment="#")
#     supp_df = pd.read_csv(supp_path)
#     supp_df.columns = supp_df.columns.str.strip()
# 
#     supp_id = id_col_supp or ("SampleID" if "SampleID" in supp_df.columns else supp_df.columns[0])
# 
#     # Clean IDs to build join keys
#     sra_df["join_key"] = standardize_series(sra_df[id_col_sra], custom_regex)
#     supp_df["join_key"] = standardize_series(supp_df[supp_id], custom_regex)
# 
#     # Deduplicate
#     supp_unique = supp_df.drop_duplicates(subset=["join_key"]).copy()
# 
#     # Subset columns
#     if clinical_cols:
#         target_cols = ["join_key"] + [c for c in clinical_cols if c in supp_unique.columns]
#         supp_unique = supp_unique[target_cols]
# 
#     # Merge
#     merged_df = pd.merge(sra_df, supp_unique, on="join_key", how="left")
#     merged_df.drop(columns=["join_key"], inplace=True)
# 
#     # Format QIIME 2 directive header
#     numeric_set = set(numeric_cols or ["age", "BMI", "ISI", "SE", "GDS_SUM", "SAST_SUM"])
#     type_row = {}
#     for col in merged_df.columns:
#         if col == "sample-id":
#             type_row[col] = "#q2:types"
#         elif col in numeric_set or pd.api.types.is_numeric_dtype(merged_df[col]):
#             type_row[col] = "numeric"
#         else:
#             type_row[col] = "categorical"
# 
#     final_df = pd.concat([pd.DataFrame([type_row]), merged_df], ignore_index=True)
#     os.makedirs(os.path.dirname(output_path), exist_ok=True)
#     final_df.to_csv(output_path, sep="\t", index=False)
#     print(f"✅ Enriched: {os.path.basename(output_path)}")
