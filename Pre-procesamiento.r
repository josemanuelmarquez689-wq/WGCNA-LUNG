# --- FASE 1: PRE-PROCESAMIENTO (TCGA TUMOR VS ADYACENTE) ---


library(TCGAbiolinks)
library(SummarizedExperiment)

# --- PASO 1: Consultar los datos en el GDC ---
# Filtramos por Proyecto (LUAD) y por el pipeline oficial 'STAR - Counts'
query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts" 
)

# 2. Descargar los datos a tu carpeta de trabajo
# Se creará una carpeta llamada 'GDCdata' automáticamente
GDCdownload(query_luad)

# --- PASO 3: Preparar el objeto en R ---
data_luad <- GDCprepare(query_luad)

# 4. Extraer la Matriz de Expresión en escala TPM
# El TPM es ideal para WGCNA porque ya está normalizado
matriz_tpm <- assay(data_luad, "tpm_unstrand")<- assay(data_luad, "tpm_unstrand")

# 5. Identificar y Separar Grupos (Tumor vs Adyacente)
# Usamos los metadatos internos de TCGAbiolinks
metadatos <- colData(data_luad)

# TP = Tumor Primario (Código 01) | NT = Normal Tissue / Adyacente (Código 11)
datExpr_Tumor_raw <- matriz_tpm[, metadatos$shortLetterCode == "TP"]
datExpr_Sano_raw  <- matriz_tpm[, metadatos$shortLetterCode == "NT"]

# 6. Verificación de seguridad
print(paste("Muestras de Tumor (TP) detectadas:", ncol(datExpr_Tumor_raw)))
print(paste("Muestras Adyacentes (NT) detectadas:", ncol(datExpr_Sano_raw)))

# Guardamos el avance en un archivo .RData para no repetir la descarga
save(datExpr_Tumor_raw, datExpr_Sano_raw, file = "Datos_TCGA_LUAD_Fase1.RData")


#Separamos por estadio

# 1. Obtener los nombres de las muestras que REALMENTE están en tu matriz de tumor
muestras_en_matriz <- colnames(datExpr_Tumor_final)

# 2. Filtrar los metadatos para que SOLO tengan esas muestras
metadatos_filtrados <- colData(data_luad)[muestras_en_matriz, ]

# 3. Ahora sí, crear los índices sobre los metadatos sincronizados
# Buscamos 'Stage I' (incluyendo IA e IB)
indices_reales_I <- grepl("Stage I$|Stage IA|Stage IB", metadatos_filtrados$ajcc_pathologic_stage)

# 4. Extraer la matriz final para el Estadio I
datExpr_StageI <- datExpr_Tumor_final[, indices_reales_I]

# 5. Verificación de éxito
print(paste("Muestras de Estadio I sincronizadas:", ncol(datExpr_StageI)))




