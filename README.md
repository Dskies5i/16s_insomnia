No puedo generar ni entregar archivos PDF directamente desde este chat (solo texto). Lo que sí puedo hacer es darte el contenido en formato **Markdown limpio tipo README.md** listo para copiar y pegar, y decirte cómo convertirlo a PDF en 10 segundos.

**Cómo obtener el PDF:**
1. Pega el contenido de abajo en un archivo llamado `README.md` (por ejemplo en RStudio: File → New File → Markdown).
2. En RStudio: botón **Knit** → elige **PDF** (requiere tener instalado `rmarkdown` + LaTeX, o `tinytex::install_tinytex()`).
3. Alternativa rápida sin LaTeX: pega el Markdown en [https://www.markdowntopdf.com](https://www.markdowntopdf.com) o en VS Code con la extensión "Markdown PDF".

Aquí está el contenido completo, con `#`, `##`, `###`, negritas y bloques de código, listo para copiar:

---

# 🧬 16s_insomnia: Microbiota Analysis in Insomnia

Bienvenido al repositorio oficial del pipeline bioinformático y de metaanálisis para el estudio del rol de la microbiota intestinal en el insomnio. Este proyecto replica y expande el diseño metodológico de Park et al. (2026), integrando un flujo robusto que abarca desde el preprocesamiento de lecturas crudas de secuenciación de amplicomes 16S rRNA hasta análisis avanzados de diversidad, aprendizaje automático (Boruta), modelado multivariado (MaAsLin2, db-RDA, LEfSe) y predicción funcional mediante PICRUSt2.

---

## 🏗️ Estructura General del Proyecto

```text
16s_insomnia/
├── 0_preprocessing/                      # Módulo de preprocesamiento y gestión de metadatos (Qiime2/DADA2)
├── 1_meta_analysis/                      # Módulo de análisis estadístico avanzado y machine learning en R
├── config/                               # Archivos de configuración central (config.yaml)
├── data/                                 # Almacenamiento de datos crudos, preprocesados y unificados
├── outputs/                              # Resultados generados (tablas, figuras, modelos y logs)
├── resources/                            # Clasificadores SILVA y recursos externos
├── utils/                                # Funciones auxiliares de soporte en R (I/O, graficación, phyloseq)
├── setup_and_run.sh                      # Script maestro de despliegue y sincronización automática
└── README.md                             # Documentación principal del proyecto
```

---

## ⚙️ Entorno de Trabajo y Prerrequisitos

Todo el pipeline se ejecuta bajo un entorno virtual centralizado de **Conda** optimizado para bioinformática (`qiime2-amplicon-2024.5`).

### Instalación Rápida con el Script Maestro (`setup_and_run.sh`)

Para desplegar y configurar automáticamente el proyecto tanto en tu PC con **Linux Mint** como en **macOS**, utiliza el asistente automatizado ubicado en la raíz:

1. Abre tu terminal en la carpeta principal del proyecto (`16s_insomnia/`).
2. Otorga permisos de ejecución (solo la primera vez):

```bash
chmod +x setup_and_run.sh
```

3. Ejecuta el script maestro:

```bash
./setup_and_run.sh
```

Este script detectará tu sistema, validará el entorno Conda `qiime2-amplicon-2024.5` y sincronizará de forma transparente todas las dependencias de R utilizando **`renv`** (`renv.lock`).

---

## 📥 0. Módulo de Preprocesamiento (FASTQ → Tablas Unificadas)

> **Nota importante sobre el flujo:** Este pipeline está diseñado para integrar inicialmente **3 estudios principales** (con la capacidad prevista de incorporar otros **5 estudios adicionales** de forma modular). Debido a esta integración progresiva y a la heterogeneidad de los metadatos públicos, el preprocesamiento no cuenta con un único script de ejecución global (`run_all`), requiriendo una supervisión y ajuste manual de los metadatos por cohorte.

### 📋 Resumen del Flujo de Preprocesamiento

1. **Adquisición:** Descarga de archivos FASTQ crudos y metadatos desde SRA utilizando SRA Toolkit.

2. **Armonización:** Normalización de variables clínicas clave (escalas ISI, IMC, regiones geográficas y estandarización de los grupos de estudio en `Insomnia` vs `Control`).

3. **Procesamiento QIIME 2 / DADA2:**
   - Remoción de primers específicos de región mediante *Cutadapt*.
   - Corrección de errores y *denoising* con *DADA2* (ajustando dinámicamente parámetros de truncamiento *Paired-End* o *Single-End*).
   - Filtrado de muestras de baja profundidad (< 10,000 lecturas retenidas).
   - Clasificación taxonómica mediante un clasificador Naive Bayes entrenado en la base de referencia **SILVA 138.1** (identidad al 99%).

4. **Integración Global:** Fusión de tablas de características mediante *VSEARCH* (*de novo* clustering al 99% de identidad) para generar un vocabulario unificado de ASVs/OTUs compartido entre todas las cohortes, exportando matrices limpias listas para el análisis en R.

Para una revisión detallada de cada script, parámetros de configuración (`datasets_config.tsv`) y diagramas de ejecución de esta fase, consulta la documentación específica en [0_preprocessing/README.md](0_preprocessing/README.md) o revisa la guía interna de la cohorte.

---

## 📊 1. Módulo de Metaanálisis y Análisis Funcional en R (`1_meta_analysis/`)

Una vez consolidada la matriz unificada (`unified_otu_table.tsv`) y el maestro de metadatos (`master_metadata_unified.tsv`), el flujo estadístico avanzado se ejecuta de manera secuencial mediante el script maestro:

```bash
conda activate qiime2-amplicon-2024.5
Rscript 1_meta_analysis/run_all.R
```

El pipeline de metaanálisis ejecuta de forma automatizada los siguientes componentes (respaldados por el diseño optimizado de Park et al.):

- **01. Construcción de Objeto Phyloseq (`01_build_phyloseq.R`):** Carga las tablas unificadas, filtra secuencias no deseadas (mitocondrias, cloroplastos o asignaciones nulas) y fija los niveles de los factores de estudio.
- **02. Diversidad Alfa (`02_alpha_diversity.R`):** Calcula índices de riqueza y diversidad (Shannon, Simpson) utilizando gráficos de violín/caja compuestos mediante *patchwork* y pruebas estadísticas t-test.
- **03. Diversidad Beta y PERMANOVA (`03_beta_diversity_permanova.R`):** Evalúa la disimilitud de Bray-Curtis ajustando por covariables clínicas multivariadas (edad, sexo, cohorte/estudio) mediante *vegan::adonis2*.
- **04. Razón Firmicutes / Bacteroidota (`04_fb_ratio.R`):** Cuantifica de forma robusta la proporción de los principales filos microbianos asociados al sueño.
- **05 & 06. Selección de Características por Machine Learning - Boruta (`05_boruta.R` & `06_boruta_pca_heatmap.R`):** Identifica ASVs discriminantes clave e integra análisis de componentes principales (PCA) junto a mapas de calor normalizados por Z-score.
- **07. Ordenación Restringida db-RDA y Envfit (`07_dbrda_envfit_fix.R`):** Ejecuta selección de modelos escalonados (*ordiR2step*) y correlaciona la matriz comunitaria con variables ambientales/clínicas corrigiendo desalineaciones de p-valores.
- **08. Biomezcladores diferenciales LEfSe (`08_lefse.R`):** Detecta biomarcadores taxonómicos robustos utilizando la API moderna de *microbiomeMarker* y genera gráficos de barras con puntuaciones LDA firmadas.
- **09. Abundancia Diferencial MaAsLin2 (`09_maaslin2_taxa.R`):** Aplica modelos lineales multivariados evitando colapsos gráficos en versiones modernas de ggplot2 mediante un bypass visual explícito.
- **10 al 15. Predicción Funcional (PICRUSt2 y Redes Metabólicas):**
  - Convierte y procesa perfiles metagenómicos predichos hacia ortólogos KEGG (KOs).
  - Calcula correlaciones vectorizadas de Spearman entre taxones clave y funciones metabólicas.
  - Genera diagramas de acorde (*Chord Diagrams* optimizados con límites de memoria RAM) y diagramas de puntos (*Dotplots*) para la visualización de rutas enriquecidas.

---

## 📌 Puntos Clave para Reportes y Tesis

1. **Robustez en Preprocesamiento:** El uso de una base unificada por VSEARCH (99%) previene sesgos derivados de llamadas de ASVs independientes por estudio.

2. **Control de Versiones y Compatibilidad:** La desactivación de motores gráficos internos obsoletos en MaAsLin2 y el manejo defensivo de memoria en diagramas circulares garantizan una ejecución libre de fallos OOM (Out Of Memory) en equipos locales.

3. **Reproducibilidad Total:** Gracias al versionamiento de entornos virtuales y la gestión automatizada de librerías mediante *renv*, cualquier análisis puede replicarse de forma idéntica en distintas arquitecturas (Linux Mint / macOS).

---

**Nota final:** Si al pegarlo en RStudio ves que los `#` y `**` no se renderizan, guarda el archivo con extensión `.md` (no `.txt`) y usa la vista "Visual Markdown" o Knit a PDF/HTML.