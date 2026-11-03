#!/usr/bin/env python3

# ============================================================
# NOMBRE: 03_merge_study_w_downloaded_metadata.py
# DESCRIPCIÓN: Combina los metadatos crudos de cada estudio (CSV o clínicos 
#              suplementarios) con los metadatos descargados de SRA.
#              Utiliza diccionarios de configuración (`STUDY_CONFIGS`) para 
#              mapear kwargs hacia el motor central de enriquecimiento.
#
# PROPÓSITO DEL FLUJO:
#   - Toma atributos SRA y archivos clínicos complementarios.
#   - Limpia y estandariza los IDs de muestra para que coincidan con las 
#     corridas de secuenciación.
#   - Armoniza fenotipos y grupos de estudio (ej. Insomnio vs. Control).
#   - Genera archivos finales enriquecidos (`_format_enriched.tsv`) listos con la 
#     fila de directivas `#q2:types` requerida por QIIME 2 para reconocer 
#     columnas numéricas y categóricas en los análisis estadísticos posteriores.
#
# USO: python 0_preprocessing/scripts/03_merge_study_w_downloaded_metadata.py
# ENV: conda activate qiime2-amplicon-2024.5
# ============================================================

import os
from utils.metadata_q2_enriched_engine import enrich_study_metadata

# ==============================================================================
# 0. CONFIGURACIÓN DE RUTAS BASE
# ==============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))

# Funciones auxiliares para simplificar la lectura de rutas en los diccionarios
path_meta = lambda filename: os.path.join(PROJECT_ROOT, "data/pre_processed/metadata", filename)
path_raw  = lambda study, filename: os.path.join(PROJECT_ROOT, f"data/raw/{study}", filename)

# ==============================================================================
# 1. CONFIGURACIÓN DE ESTUDIOS (STUDY_CONFIGS)
# ==============================================================================
STUDY_CONFIGS = [
    
    # --- Study 1: Haimov 2022 (Modo Combinado / Join Mode) ---
    {
        "name": "haimov_2022",
        "sra_path": path_meta("haimov_2022_metadata_q2.tsv"),
        "supp_path": path_raw("haimov_2022", "table_s1.csv"),
        "output_path": path_meta("haimov_2022_metadata_q2_format_enriched.tsv"),
        "id_col_sra": "library_name",
        "id_col_supp": "SampleID",
        "custom_regex": r"^([A-Za-z0-9]+(?:-[0-9]+)?)",
        "clinical_cols": ["BMI", "ISI", "ISI_cat", "SE", "SE_cat", "GDS_SUM", "SAST_SUM"],
        "numeric_cols": ["age", "BMI", "ISI", "SE", "GDS_SUM", "SAST_SUM"],
        "value_mappings": {
            "replace_values": {
                "SE_cat": {"Ins": "Insomnia", "nIns": "Control"},
                "tissue": {"feces": "stool"}
            },
            "rename_cols": {
                "SE_cat": "group"
            }
        }
    },
    
    # --- Study 2: Chen 2022 (Modo Directo SRA / Direct Mode) ---
    {
        "name": "chen_2022",
        "sra_path": path_meta("chen_2022_metadata_q2.tsv"),
        "supp_path": None,  # Sin archivo complementario externo
        "output_path": path_meta("chen_2022_metadata_q2_format_enriched.tsv"),
        "id_col_sra": "library_name",
        "numeric_cols": ["age"],
        "group_prefix_map": {
            "A": "Insomnia",
            "B": "Insomnia",  # Mapeado temporalmente; ajustar si se requiere subgrupo específico
            "C": "Insomnia"
        },
        "value_mappings": {
            "replace_values": {
                "tissue": {"Intestinal tract": "stool"}
            },
            "rename_cols": {
                "group": "group_detailed"
            },
            "derive_macro_group": {
                "source_col": "group_detailed",
                "mapping": {
                    "Insomnia": "Insomnia",
                    "Control": "Control"
                }
            }
        }
    
    },
    
    # --- Study 3: Liub 2019 (Modo Directo SRA / Derivado por Prefijo) ---
    {
        "name": "liub_2019",
        "sra_path": path_meta("liub_2019_metadata_q2.tsv"),
        "supp_path": None,  # No se requiere archivo complementario externo
        "output_path": path_meta("liub_2019_metadata_q2_format_enriched.tsv"),
        "id_col_sra": "library_name",
        # Usamos una función o mapeo de prefijos para asignar el grupo basándonos en library_name
        "group_prefix_map": {
            "insomnia": "Insomnia",
            "normal": "Control"
        },
        "value_mappings": {
            "rename_cols": {
                "group": "group"
            }
        }
    }
    
    
]

# ==============================================================================
# 2. EJECUCIÓN DEL PROCESO POR LOTE
# ==============================================================================
def main():
    print("🔍 Iniciando enriquecimiento de metadatos por lote...")
    
    for config in STUDY_CONFIGS:
        study_name = config.get("name", "Estudio Desconocido")
        print(f"\n📂 Procesando estudio: {study_name}...")
        
        # Extraer argumentos excluyendo la clave interna 'name'
        kwargs = {k: v for k, v in config.items() if k != "name"}
        
        # Ejecutar motor central
        enrich_study_metadata(**kwargs)

    print("\n🎉 ¡Proceso de enriquecimiento de metadatos finalizado!")


if __name__ == "__main__":
    main()

# # ============================================================
# # NOMBRE: 03_merge_study_w_downloaded_metadata.py
# # DESCRIPCIÓN: Combina los metadatos crudos de cada estudio (CSV o clínicos 
# #              suplementarios) con los metadatos descargados de SRA.
# #              Utiliza diccionarios de configuración (`STUDY_CONFIGS`) para 
# #              mapear kwargs hacia el motor central de enriquecimiento.
# #
# # PROPÓSITO DEL FLUJO:
# #   - Toma atributos SRA y archivos clínicos complementarios.
# #   - Limpia y estandariza los IDs de muestra para que coincidan con las 
# #     corridas de secuenciación.
# #   - Armoniza fenotipos y grupos de estudio (ej. Insomnio vs. Control).
# #   - Genera archivos finales enriquecidos (`_enriched.tsv`) listos con la 
# #     fila de directivas `#q2:types` requerida por QIIME 2 para reconocer 
# #     columnas numéricas y categóricas en los análisis estadísticos posteriores.
# #
# # USO: python 0_preprocessing/scripts/03_merge_study_w_downloaded_metadata.py
# # ============================================================
# 
# 
# #!/usr/bin/env python3
# import os
# #conda install -n qiime2-amplicon-2024.5 -c conda-forge pandas    ---- OR: pip install pandas
# # python -c "import pandas; print(pandas.__version__)"
# #import pandas as pd
# from utils.metadata_q2_enriched_engine import enrich_study_metadata
# 
# SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))
# 
# # ==============================================================================
# # 1. THIS SAVES KWARGS DICTIONARIES
# # ==============================================================================
# STUDY_CONFIGS = [
#       # --- Study 1: Haimov 2022 ---   
#     {
#         "name": "haimov_2022",
#         "sra_path": os.path.join(PROJECT_ROOT, "data/pre_processed/metadata/haimov_2022_metadata_q2.tsv"),
#         "supp_path": os.path.join(PROJECT_ROOT, "data/raw/haimov_2022/table_s1.csv"),
#         "output_path": os.path.join(PROJECT_ROOT, "data/pre_processed/metadata/haimov_2022_metadata_q2_format_enriched.tsv"),
#         "id_col_sra": "library_name",
#         "id_col_supp": "SampleID",
#         "custom_regex": r"^([A-Za-z0-9]+(?:-[0-9]+)?)",
#         "clinical_cols": ["BMI", "ISI", "ISI_cat", "SE", "SE_cat", "GDS_SUM", "SAST_SUM"],
#         "numeric_cols": ["age", "BMI", "ISI", "SE", "GDS_SUM", "SAST_SUM"],
#         "value_mappings": {
#             # Step A: Standardise Ins / nIns values
#             "replace_values": {
#                 "SE_cat": {
#                     "Ins": "Insomnia",
#                     "nIns": "Control"
#                 },
#                 "tissue": {"feces": "stool"}
#             },
#             # Step B: Rename SE_cat to group
#             "rename_cols": {
#                 "SE_cat": "group"
#             }
#         }
#     },
#     
#      # --- Study 2: Chen 2022 ---
#     {
#         "name": "chen_2022",
#         "sra_path": os.path.join(PROJECT_ROOT, "data/pre_processed/metadata/chen_2022_metadata_q2.tsv"),
#         "supp_path": None,  # Direct SRA mode
#         "output_path": os.path.join(PROJECT_ROOT, "data/pre_processed/metadata/chen_2022_metadata_q2_format_enriched.tsv"),
#         "id_col_sra": "library_name",
#         "numeric_cols": ["age"],
#         # Define A, B, and C prefixes
#         "group_prefix_map": {
#             "A": "Insomnia",
#             "B": "Insomnia",  # Replace 'Subgroup_B' once your investigation confirms the exact phenotype (e.g., Paradoxical_Insomnia)
#             "C": "Insomnia"
#         },
#         "value_mappings": {
#             "replace_values": {
#                 "tissue": {"Intestinal tract": "stool"}
#             },
#             # Convert direct group assignment to group_detailed so macro mapping can collapse it
#             "rename_cols": {
#                 "group": "group_detailed"
#             },
#             "derive_macro_group": {
#                 "source_col": "group_detailed",
#                 "mapping": {
#                     "Insomnia": "Insomnia",
#                     #"Subgroup_B": "Insomnia",  # Maps Subgroup B into primary Insomnia macro-group (or change as needed)
#                     "Control": "Control"
#                 }
#             }
#         }
#     }
#     
#     #   # --- Example template for 3-group studies (e.g., DSM-5 Paradoxical Insomnia cohorts)
#     # {
#     #     "name": "three_group_study_202X",
#     #     "sra_path": os.path.join(PROJECT_ROOT, "data/pre_processed/metadata/study3_metadata_q2.tsv"),
#     #     "supp_path": os.path.join(PROJECT_ROOT, "data/raw/study3/supplement.csv"),
#     #     "output_path": os.path.join(PROJECT_ROOT, "data/pre_processed/metadata/study3_metadata_q2_format_enriched.tsv"),
#     #     "id_col_sra": "library_name",
#     #     "id_col_supp": "SampleID",
#     #     "value_mappings": {
#     #         "replace_values": {
#     #             "diagnosis_col": {
#     #                 "Psychophysiological Insomnia": "Insomnia",
#     #                 "Paradoxical Insomnia": "Paradoxical_Insomnia",
#     #                 "Healthy Control": "Control"
#     #             }
#     #         },
#     #         "rename_cols": {
#     #             "diagnosis_col": "group_detailed"
#     #         },
#     #         "derive_macro_group": {
#     #             "source_col": "group_detailed",
#     #             "mapping": {
#     #                 "Insomnia": "Insomnia",
#     #                 "Paradoxical_Insomnia": "Insomnia",
#     #                 "Control": "Control"
#     #             }
#     #         }
#     #     }
#     # }
# ]
# 
# # ==============================================================================
# # 2. THIS IS WHERE THE DICTIONARY GETS CONVERTED TO KWARGS
# # ==============================================================================
# 
# def main():
#     print("🔍 Starting batch metadata enrichment...")
#     for config in STUDY_CONFIGS:
#         study_name = config.get("name", "Unnamed Study")
#         print(f"\nProcessing study: {study_name}...")
#         
#         # Step A: Take the study dictionary and drop the internal 'name' key
#         kwargs = {k: v for k, v in config.items() if k != "name"}
#         
#         # Step B: Pass the kwargs dictionary into the engine function
#         enrich_study_metadata(**kwargs)
# 
# 
# if __name__ == "__main__":
#     main()
