# ============================================================
# NOMBRE: 01_build_phyloseq.R
# DESCRIPCIÓN: Construcción, filtrado y normalización del objeto 
#              phyloseq global a partir de las matrices unificadas 
#              de QIIME 2 (VSEARCH 99%) y metadatos maestros.
#
# ¿QUÉ HACE ESTE SCRIPT?
#    1. Lee la tabla OTU/ASV unificada, la taxonomía SILVA y la 
#       tabla maestra de metadatos desde data/pre_processed/unified/.
#    2. Construye el objeto phyloseq integrando conteos, taxonomía y muestras.
#    3. Aplica filtros de calidad (elimina taxones no deseados como 'g__other' 
#       o 'Unassigned' y remueve sumas en cero).
#    4. Fija y estandariza los niveles de los grupos de estudio (Control vs. Insomnia).
#    5. Guarda el objeto phyloseq resultante y exporta las tablas filtradas 
#       de soporte a outputs/models/ y outputs/tables/.
#
# ARCHIVOS DE SALIDA (OUTPUTS):
#    - outputs/models/physeq_filtered.rds       (Objeto phyloseq completo y filtrado en R)
#    - outputs/tables/tax_filtered.csv          (Tabla taxonómica limpia y filtrada)
#    - outputs/tables/otu_filtered.csv          (Matriz de conteos filtrada)
#    - outputs/tables/metadata_filtered.csv     (Tabla maestra de metadatos filtrada)
#
# USO: Rscript 1_meta_analysis/01_build_phyloseq.R
#
# ENVIRONMENT: conda activate qiime2-amplicon-2024.5
# phyloseq is instaleld inside "qiime2-amplicon-2024.5" environment
#
# WORKS? -- YES
# ============================================================
source("1_meta_analysis/00_setup.R")
library(phyloseq)

ps <- build_phyloseq()

if ("Genus" %in% colnames(tax_table(ps))) {
  ps <- subset_taxa(ps, !grepl("g__other|Unassigned", tax_table(ps)[, "Genus"]))
}
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)
ps <- enforce_group_levels(ps, cfg$columns$group, cfg$groups$levels)

save_rds_safe(ps, file.path(OUT_DIR, "models/physeq_filtered.rds"))

tax  <- as.data.frame(tax_table(ps))    %>% tibble::rownames_to_column("OTU")
otu  <- as.data.frame(otu_table(ps))    %>% tibble::rownames_to_column("OTU")
meta <- as.data.frame(sample_data(ps))  %>% tibble::rownames_to_column("Sample")

write_csv_safe(tax,  file.path(OUT_DIR, "tables/tax_filtered.csv"))
write_csv_safe(otu,  file.path(OUT_DIR, "tables/otu_filtered.csv"))
write_csv_safe(meta, file.path(OUT_DIR, "tables/metadata_filtered.csv"))

cat("✅ Objeto phyloseq construido y guardado en outputs/models/\n")
cat("   Taxa:", ntaxa(ps), "| Samples:", nsamples(ps), "\n")


