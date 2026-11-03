#00_setup.R
suppressPackageStartupMessages({
  library(yaml)
  library(fs)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# Como corres desde la raíz (16s_insomnia/), el config está directamente en config/config.yaml
cfg <- yaml::read_yaml("config/config.yaml")

set.seed(cfg$project$seed)

# La carpeta de outputs está directamente en outputs/ en la raíz
ROOT_DIR <- "."
OUT_DIR <- file.path(ROOT_DIR, cfg$outputs$dir)

# Crear directorios de salida de manera segura
dir_create(file.path(OUT_DIR, "tables"))
dir_create(file.path(OUT_DIR, "figures"))
dir_create(file.path(OUT_DIR, "models"))
dir_create(file.path(OUT_DIR, "logs"))

# Cargar funciones auxiliares desde la carpeta 'utils/' en la raíz
source("utils/utils_io.R")
source("utils/utils_phyloseq.R")
source("utils/utils_plot.R")
source("utils/utils_stats.R")
source("utils/utils_helper_functions.R")