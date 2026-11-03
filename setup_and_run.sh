#!/usr/bin/env bash
# ==============================================================================
#    NAME: setup_and_run.sh
# SCRIPT MAESTRO DE DESPLIEGUE Y EJECUCIÓN (Linux Mint / macOS)
# Proyecto: 16s_insomnia
# ==============================================================================
# DESCRIPTION:
#   Este script maestro automatiza el despliegue, la validación de entornos de 
#   Conda y la sincronización estricta de dependencias de R en proyectos de 
#   microbioma. Actúa como un puente inteligente adaptativo entre sistemas operativos.
#
# WHEN TO USE (¿Cuándo usarlo?):
#   - Al clonar o copiar por primera vez el proyecto en una nueva computadora o laptop.
#   - Cada vez que cambies de equipo (por ejemplo, al migrar tu repositorio de tu Mac 
#     a tu PC con Linux Mint) para asegurar que el entorno esté listo inmediatamente.
#   - Al retomar el trabajo diario en tu equipo actual para verificar directorios 
#     y dependencias antes de lanzar los análisis.
#
# WHY IS THIS GOOD FOR MIGRATING TO ANOTHER PC/LAPTOP (¿Por qué es ideal para migrar?):
#   - Reproducibilidad cross-platform: Utiliza la tecnología de 'renv' y el archivo 
#   - 'renv.lock' para congelar y restaurar las versiones exactas de todos los paquetes 
#     de R, evitando conflictos de compatibilidad entre macOS y Linux Mint.
#   - Automatización total: Detecta de forma autónoma si el proyecto se configura por 
#     primera vez (creando el lockfile) o si se está abriendo en un equipo nuevo 
#     (restaurando las librerías automáticamente sin intervención manual tediosa).
#   - Comprobación de seguridad: Valida la presencia obligatoria del entorno virtual 
#     de Conda (qiime2-amplicon-2024.5) y genera la estructura completa de directorios 
#     de salida requerida por el pipeline.
#
# USAGE (Uso básico):
#   1. Dar permisos de ejecución (solo la primera vez):
#      chmod +x setup_and_run.sh
#   2. Ejecutar el script:
#      ./setup_and_run.sh
# ==============================================================================
# COMANDOS DE SOPORTE MANUAL (Opcionales si se requiere forzar actualización):
#   - conda activate qiime2-amplicon-2024.5
#   - Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("BiocVersion", update = FALSE, ask = FALSE); renv::snapshot(prompt = FALSE)'
#   - set -eo pipefail
# ==============================================================================

# Colores para mensajes en consola
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # Sin color

echo -e "${BLUE}==================================================${NC}"
echo -e "${BLUE}  🚀 INICIANDO CONFIGURACIÓN DEL ENTORNO 16S_INSOMNIA${NC}"
echo -e "${BLUE}==================================================${NC}"

# 1. Verificar si Conda está instalado
if ! command -v conda &> /dev/null; then
    echo -e "${RED}❌ Error: Conda no está instalado o no se encuentra en el PATH.${NC}"
    exit 1
fi

# Inicializar conda para el script de bash
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# 2. Verificar o activar el entorno qiime2-amplicon-2024.5
ENV_NAME="qiime2-amplicon-2024.5"
if ! conda env list | grep -q "${ENV_NAME}"; then
    echo -e "${RED}❌ El entorno de Conda '${ENV_NAME}' no existe.${NC}"
    echo -e "Por favor créalo primero usando tu archivo .yml o la configuración adecuada."
    exit 1
fi

echo -e "${GREEN}✅ Activando entorno Conda: ${ENV_NAME}${NC}"
conda activate "${ENV_NAME}"

# 3. Verificar directorios clave de la estructura
echo -e "${BLUE}📁 Verificando estructura de directorios...${NC}"
mkdir -p outputs/tables outputs/figures outputs/models outputs/logs outputs/function/figures outputs/function/tables outputs/boruta outputs/lefse outputs/maaslin2 data/pre_processed/unified resources/classifiers

# 4. Inicializar o Restaurar entorno R con renv de forma automática
echo -e "${BLUE}🔄 Sincronizando paquetes de R mediante renv...${NC}"
Rscript -e '
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
if (file.exists("renv.lock")) {
  message("📂 renv.lock encontrado. Restaurando entorno...")
  renv::restore(prompt = FALSE)
} else {
  message("⚠️ No se encontró renv.lock. Inicializando renv por primera vez...")
  renv::init(restart = FALSE)
}
'

echo -e "${GREEN}✅ Entorno R configurado correctamente.${NC}"

# 5. Menú de ejecución limpio (sin duplicados)
echo -e "\n${BLUE}==================================================${NC}"
echo -e "${BLUE}  ⚙️  ¿Qué deseas ejecutar?${NC}"
echo -e "  1) Ver guía / Recordatorio de Preprocesamiento (Manual)"
echo -e "  2) Ejecutar Metaanálisis y Análisis Funcional completo (R)"
echo -e "${BLUE}==================================================${NC}"
read -p "Selecciona una opción [1-2]: " option

case $option in
    1)
        echo -e "${BLUE}ℹ️  Módulo de Preprocesamiento (0_preprocessing):${NC}"
        echo -e "Recuerda que el preprocesamiento se ejecuta de forma manual cohorte por cohorte"
        echo -e "debido a la integración progresiva de estudios y ajuste de metadatos."
        echo -e "Revisa los scripts en ${GREEN}0_preprocessing/scripts/${NC} según lo requieras."
        ;;
    2)
        echo -e "${BLUE}▶️ Ejecutando pipeline completo de metaanálisis (R)...${NC}"
        Rscript 1_meta_analysis/run_all.R
        echo -e "${GREEN}🎉 ¡Metaanálisis y análisis funcional completados con éxito!${NC}"
        ;;
    *)
        echo -e "${RED}Opción inválida. Saliendo.${NC}"
        exit 1
        ;;
esac