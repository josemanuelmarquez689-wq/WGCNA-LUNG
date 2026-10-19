# ==============================================================================
# ROBUSTEZ DEL MODULO: CONSERVACION POR JACCARD MAXIMO (red de 24 modulos)
# Extraido de la parte validada de ../WGCNA/Ataque_masivo_2.R (barrido fino +
# ataques por % + Matriz_Conservacion_Real + criticidad).
# Para cada ataque (modulo x 1-15 hubs / 5-10-20-90%):
#   1. Se quitan los hubs (orden por kWithin) y se reconstruye la red desde cero:
#      adyacencia -> TOM -> hclust -> cutreeDynamic (deepSplit 2, minClusterSize 30).
#   2. Matriz de Jaccard: modulos originales x modulos de la red alterada,
#      usando solo los genes sobrevivientes y sin "grey".
#   3. Conservacion del modulo atacado = Jaccard MAXIMO de su fila
#      (el modulo de la red alterada con el que mejor coincide).
# Power = 12: es el que reproduce los resultados guardados de la red de 24
# modulos (BLACK 1H = 0.244318); Ataque_masivo_2.R hoy dice 9, pero eso es de
# otra red.
# ==============================================================================
.libPaths(c("/home/jmarquez/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
suppressMessages({ library(WGCNA); library(pheatmap); library(parallel) })
disableWGCNAThreads()
setwd("/STORAGE/genut/Manuel/directed_attack_network")
load("data/datExpr_colores_24Modulos.RData")

POWER_OPTIMO <- 12; TYPE_NETWORK <- "signed"
DEEP_SPLIT <- 2; MIN_CLUSTER_SIZE <- 30
NUCLEOS <- as.integer(Sys.getenv("NUCLEOS", "20"))
PRUEBA  <- Sys.getenv("PRUEBA", "0") == "1"
DIRECTORIO_SALIDA <- if (PRUEBA) "results/prueba_conservacion" else "results/conservacion"
dir.create(DIRECTORIO_SALIDA, recursive = TRUE, showWarnings = FALSE)

colores_base <- colores_base[colnames(datExpr_5000)]   # mismo orden que datExpr
adj_base <- adjacency(datExpr_5000, power = POWER_OPTIMO, type = TYPE_NETWORK)
Alldegrees <- intramodularConnectivity(adj_base, colores_base)
genes <- colnames(datExpr_5000)

modulos <- sort(setdiff(unique(colores_base), "grey"))
columnas <- c(paste0(1:15, " H"), "5%", "10%", "20%", "90%")
if (PRUEBA) { modulos <- c("black", "purple"); columnas <- c("1 H", "3 H", "5%") }
tamanos_modulos <- table(colores_base)

hubs_de <- function(mod, col) {
  g <- names(colores_base[colores_base == mod])
  df <- Alldegrees[rownames(Alldegrees) %in% g, ]
  df <- df[order(df$kWithin, decreasing = TRUE), ]
  n <- if (grepl("%", col)) max(1, round(length(g) * as.numeric(sub("%", "", col)) / 100))
       else as.integer(sub(" H", "", col))
  rownames(df)[1:min(n, nrow(df))]
}

reconstruir <- function(rem) {
  keep <- which(!(genes %in% rem))
  dissTOM <- 1 - TOMsimilarity(adj_base[keep, keep], TOMType = TYPE_NETWORK, verbose = 0)
  arbol <- hclust(as.dist(dissTOM), method = "average")
  mods <- cutreeDynamic(dendro = arbol, distM = dissTOM, deepSplit = DEEP_SPLIT,
                        pamRespectsDendro = FALSE, minClusterSize = MIN_CLUSTER_SIZE, verbose = 0)
  setNames(labels2colors(mods), genes[keep])
}

matriz_jaccard <- function(colores_pert) {
  comunes <- intersect(names(colores_base), names(colores_pert))
  cb <- colores_base[comunes]; cp <- colores_pert[comunes]
  mb <- sort(unique(cb[cb != "grey"])); mp <- sort(unique(cp[cp != "grey"]))
  J <- matrix(0, length(mb), length(mp), dimnames = list(mb, mp))
  for (a in mb) { ga <- names(cb[cb == a])
    for (b in mp) { gb <- names(cp[cp == b])
      u <- length(union(ga, gb)); J[a, b] <- if (u > 0) length(intersect(ga, gb)) / u else 0 } }
  J
}

trabajos <- expand.grid(col = columnas, mod = modulos, stringsAsFactors = FALSE)
cat(sprintf("Re-clusterizando %d redes alteradas con %d nucleos...\n", nrow(trabajos), NUCLEOS))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(trabajos)), function(t) {
  mod <- trabajos$mod[t]; rem <- hubs_de(mod, trabajos$col[t])
  cp <- reconstruir(rem); J <- matriz_jaccard(cp)
  fila <- J[mod, ]
  list(J = J, fila = data.frame(Modulo = toupper(mod), Ataque = trabajos$col[t], N_Removidos = length(rem),
       Conservacion = max(fila), Mejor_Match = names(fila)[which.max(fila)],
       N_Modulos_Red_Alterada = length(unique(cp[cp != "grey"])),
       Hubs_Removidos = paste(rem, collapse = ", ")))
}, mc.cores = NUCLEOS, mc.preschedule = FALSE)
fallos <- sapply(res, function(r) inherits(r, "try-error") || is.null(r))
if (any(fallos)) stop(sprintf("%d ataques fallaron", sum(fallos)))
cat(sprintf("Tiempo: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

df <- do.call(rbind, lapply(res, `[[`, "fila"))
matrices <- setNames(lapply(res, `[[`, "J"), paste(df$Modulo, df$Ataque, sep = "_"))
saveRDS(matrices, file.path(DIRECTORIO_SALIDA, "Matrices_Jaccard_Por_Ataque.rds"))
write.csv(df, file.path(DIRECTORIO_SALIDA, "Conservacion_Largo.csv"), row.names = FALSE)

M <- matrix(NA, length(modulos), length(columnas), dimnames = list(toupper(modulos), columnas))
for (r in seq_len(nrow(df))) M[df$Modulo[r], df$Ataque[r]] <- df$Conservacion[r]
write.csv(M, file.path(DIRECTORIO_SALIDA, "Matriz_Conservacion_Real.csv"))

# ---- Jaccard maximo por modulo (el ataque donde mejor se conserva) ----
jmax <- data.frame(Modulo = rownames(M), Tamano = as.integer(tamanos_modulos[tolower(rownames(M))]),
                   Jaccard_Maximo = round(apply(M, 1, max), 3),
                   Ataque = colnames(M)[apply(M, 1, which.max)])
jmax <- jmax[order(-jmax$Jaccard_Maximo), ]
write.csv(jmax, file.path(DIRECTORIO_SALIDA, "Jaccard_Maximo_Por_Modulo.csv"), row.names = FALSE)

# ---- Criticidad: mediana de la conservacion en hubs 3-7 ----
ventana <- intersect(c("3 H", "4 H", "5 H", "6 H", "7 H"), columnas)
crit_val <- apply(M[, ventana, drop = FALSE], 1, median, na.rm = TRUE)
criticidad <- cut(crit_val, c(-Inf, 0.3, 0.5, 0.7, Inf), right = FALSE,
                  labels = c("Extremadamente Critico", "Critico", "Moderado", "Robusto"))
write.csv(data.frame(Modulo = rownames(M), Mediana_Hubs_3a7 = round(crit_val, 3), Criticidad = criticidad,
                     Conservacion_1H = round(M[, "1 H"], 3))[order(crit_val), ],
          file.path(DIRECTORIO_SALIDA, "Criticidad_Modulos.csv"), row.names = FALSE)

# ---- Heatmap ----
orden <- order(crit_val)
Mh <- M[orden, , drop = FALSE]
anot <- data.frame(Criticidad = criticidad[orden])
rownames(Mh) <- rownames(anot) <- sprintf("%s (%d)", rownames(Mh), tamanos_modulos[tolower(rownames(Mh))])
paleta <- colorRampPalette(c("#8B0000", "#CC0000", "#FF4500", "#FF8C00", "#FFA500",
  "#FFD700", "#FFFF00", "#ADFF2F", "#7CFC00", "#32CD32", "#00AA00", "#006400"))(100)
pdf(file.path(DIRECTORIO_SALIDA, "Heatmap_Conservacion_Modulo.pdf"),
    width = max(16, (ncol(Mh) * 35 + 400) / 72), height = max(6, (nrow(Mh) * 30 + 240) / 72))
tryCatch(pheatmap(Mh, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
         number_format = "%.2f", number_color = "black", fontsize_number = 8,
         color = paleta, breaks = seq(0, 1, length.out = 101), angle_col = 45,
         cellwidth = 35, cellheight = 30, fontsize = 10, border_color = "grey60",
         annotation_row = anot,
         annotation_colors = list(Criticidad = c("Extremadamente Critico" = "#8B0000", "Critico" = "#CC0000",
                                                 "Moderado" = "#FFA500", "Robusto" = "#006400")),
         main = "Conservacion del modulo atacado: Jaccard MAXIMO contra los modulos de la red re-clusterizada\n(Criticidad = mediana en hubs 3-7 | ROJO = el modulo colapsa | VERDE = resiste)"),
         finally = dev.off())

print(jmax, row.names = FALSE)
cat("LISTO\n")
