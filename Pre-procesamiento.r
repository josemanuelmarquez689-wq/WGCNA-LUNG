# --- FASE 1: PRE-PROCESAMIENTO (TCGA TUMOR VS ADYACENTE) ---

# 1. Cargar librerías
library(WGCNA)

datos_expresion <- read.delim("LUNG_subtypes_summary.tsv", header = TRUE, row.names = 1)

# Filtramos usando los códigos de vial de TCGA
indices_tumor <- grepl("-01A|-01B", nombres_muestras) # Tumor primario
indices_sano  <- grepl("-11A|-11B", nombres_muestras) # Tejido adyacente
d
# 3. Crear las sub-matrices
datExpr_Tumor_raw <- datos_expresion[, indices_tumor]
datExpr_Sano_raw  <- datos_expresion[, indices_sano]

# 4. Transponer para cumplir con el formato WGCNA (Muestras como filas)
datExpr_Tumor <- as.data.frame(t(datExpr_Tumor_raw))
datExpr_Sano  <- as.data.frame(t(datExpr_Sano_raw))

# 5. Verificación de seguridad
print(paste("Muestras de Tumor (tcgaPrimary):", nrow(datExpr_Tumor)))
print(paste("Muestras Adyacentes (tcgaAdjacent):", nrow(datExpr_Sano)))

# 6. Limpieza de genes con WGCNA
library(WGCNA)
gsg_t <- goodSamplesGenes(datExpr_Tumor, verbose = 3)
gsg_s <- goodSamplesGenes(datExpr_Sano, verbose = 3)

if (!gsg_t$allOK) datExpr_Tumor = datExpr_Tumor[gsg_t$goodSamples, gsg_t$goodGenes]
if (!gsg_s$allOK) datExpr_Sano = datExpr_Sano[gsg_s$goodSamples, gsg_s$goodGenes]

save(datExpr_Tumor, datExpr_Sano, file = "Datos_Filtrados_TCGA.RData")