# ==============================================================================
# CONTROL ESTRICTO: ataque a hubs vs. quitar el MISMO numero de genes AL AZAR
# DEL MISMO MODULO. Responde: ¿los hubs de X danan la red mas que cualquier gen de X?
# Reutiliza los ataques de 01_robustez_global_integridad.R (Resultados_Integridad_Largo.csv)
# y el mismo grafo (TOM >= umbral top 1% de la red original).
# Nota: en 90% casi todo el modulo se quita en ambos casos -> z ~ 0 esperado.
# ==============================================================================
.libPaths(c("/home/jmarquez/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
suppressMessages({ library(WGCNA); library(igraph); library(pheatmap); library(parallel) })
disableWGCNAThreads()
setwd("/STORAGE/genut/Manuel/directed_attack_network")
load("data/datExpr_colores_24Modulos.RData")

POWER_OPTIMO <- 12; TYPE_NETWORK <- "signed"
N_REP   <- as.integer(Sys.getenv("N_REP", "30"))
NUCLEOS <- as.integer(Sys.getenv("NUCLEOS", "24"))
PRUEBA  <- Sys.getenv("PRUEBA", "0") == "1"
DIR_INT <- if (PRUEBA) "results/prueba_integridad" else "results/integridad"
DIRECTORIO_SALIDA <- if (PRUEBA) "results/prueba_control_modulo" else DIR_INT
dir.create(DIRECTORIO_SALIDA, showWarnings = FALSE)
set.seed(2027)

A <- read.csv(file.path(DIR_INT, "Resultados_Integridad_Largo.csv"), check.names = FALSE)
orig <- read.csv(file.path(DIR_INT, "Red_Original_Metricas.csv"))
UMBRAL <- orig$Umbral_TOM
genes <- colnames(datExpr_5000); N <- length(genes)
adj_base <- adjacency(datExpr_5000, power = POWER_OPTIMO, type = TYPE_NETWORK)

metricas <- function(TOM) {
  n <- ncol(TOM)
  g <- graph_from_adjacency_matrix((TOM >= UMBRAL) * 1, mode = "undirected", diag = FALSE)
  cmp <- components(g)
  c(Eficiencia = global_efficiency(g) * n * (n - 1) / (N * (N - 1)),
    LCC = max(cmp$csize), N_Componentes = sum(cmp$csize > 1))
}
red_sin <- function(rem) {
  keep <- which(!(genes %in% rem))
  metricas(TOMsimilarity(adj_base[keep, keep], TOMType = TYPE_NETWORK, verbose = 0))
}

# combinaciones unicas (modulo, n) de los ataques
comb <- unique(A[, c("Modulo", "n")])
if (PRUEBA) comb <- comb[comb$n %in% c(1, 15), ]
azar <- comb[rep(seq_len(nrow(comb)), each = N_REP), ]
azar$rep <- rep(seq_len(N_REP), nrow(comb))
azar$rem <- lapply(seq_len(nrow(azar)), function(i)
  sample(names(colores_base[colores_base == tolower(azar$Modulo[i])]), azar$n[i]))
cat(sprintf("Controles del mismo modulo: %d redes (%d combinaciones modulo x n, %d rep), %d nucleos\n",
            nrow(azar), nrow(comb), N_REP, NUCLEOS))

t0 <- Sys.time()
res <- mclapply(azar$rem, red_sin, mc.cores = NUCLEOS, mc.preschedule = FALSE)
fallos <- sapply(res, function(r) inherits(r, "try-error") || is.null(r))
if (any(fallos)) stop(sprintf("%d controles fallaron", sum(fallos)))
R <- cbind(azar[, c("Modulo", "n", "rep")], do.call(rbind, res))
cat(sprintf("Tiempo: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
R$Eficiencia_Rel <- R$Eficiencia / orig$Eficiencia
R$LCC_Rel <- R$LCC / orig$LCC
write.csv(R, file.path(DIRECTORIO_SALIDA, "Controles_MismoModulo.csv"), row.names = FALSE)

# ---- z y p empirico por ataque ----
A <- merge(A, comb, by = c("Modulo", "n"))
for (i in seq_len(nrow(A))) {
  r <- R[R$Modulo == A$Modulo[i] & R$n == A$n[i], ]
  A$EfMod_media[i] <- mean(r$Eficiencia_Rel); A$EfMod_sd[i] <- sd(r$Eficiencia_Rel)
  A$P_Mod[i] <- (1 + sum(r$Eficiencia_Rel <= A$Eficiencia_Rel[i])) / (1 + N_REP)
}
# sd = 0 (p. ej. n = tamano del modulo: todas las muestras iguales) -> z no definido
A$Z_Mod <- ifelse(A$EfMod_sd > 0, (A$EfMod_media - A$Eficiencia_Rel) / A$EfMod_sd, NA)
A$Dano_Hubs <- 1 - A$Eficiencia_Rel; A$Dano_Azar_Mod <- 1 - A$EfMod_media
A$Razon_Dano <- A$Dano_Hubs / A$Dano_Azar_Mod
write.csv(A, file.path(DIRECTORIO_SALIDA, "Resultados_Control_MismoModulo_Largo.csv"), row.names = FALSE)

# ---- Heatmaps ----
tam <- table(colores_base)
columnas <- c(paste0(1:15, " H"), "5%", "10%", "20%", "90%")
columnas <- columnas[columnas %in% A$Ataque]
mods <- sort(unique(A$Modulo))
mk <- function(v) { m <- matrix(NA, length(mods), length(columnas), dimnames = list(mods, columnas))
  for (r in seq_len(nrow(A))) if (A$Ataque[r] %in% columnas) m[A$Modulo[r], A$Ataque[r]] <- A[[v]][r]; m }
div <- colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"))(100)
sub_t <- sprintf("Grafo: aristas TOM >= %.4f (top 1%% de la red original)", UMBRAL)

graficar <- function(m, etq, archivo, titulo, lim, centro) {
  mc <- pmin(pmax(m, centro - lim), centro + lim)
  o <- order(rowMeans(mc, na.rm = TRUE), decreasing = TRUE)
  m <- m[o, , drop = FALSE]; mc <- mc[o, , drop = FALSE]; etq <- etq[o, , drop = FALSE]
  rownames(mc) <- rownames(m) <- sprintf("%s (%d)", rownames(m), tam[tolower(rownames(m))])
  write.csv(m, sub("\\.pdf$", ".csv", archivo))
  pdf(archivo, width = max(18, (ncol(m) * 40 + 400) / 72), height = max(7, (nrow(m) * 30 + 260) / 72))
  tryCatch(pheatmap(mc, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = etq,
           number_color = "black", fontsize_number = 7, color = div, na_col = "grey85",
           breaks = seq(centro - lim, centro + lim, length.out = 101), main = titulo, angle_col = 45,
           cellwidth = 40, cellheight = 30, fontsize = 10, border_color = "grey60"),
           finally = dev.off())
}
z <- mk("Z_Mod"); p <- mk("P_Mod")
etq_z <- matrix(ifelse(is.na(z), "NA", sprintf("%.1f%s", z, ifelse(p <= 1 / (1 + N_REP), "*", ""))), nrow(z))
graficar(z, etq_z, file.path(DIRECTORIO_SALIDA, "Heatmap_Dano_vs_Azar_MismoModulo.pdf"),
  sprintf("CONTROL ESTRICTO: dano de los hubs vs. quitar el MISMO numero de genes al azar DEL MISMO MODULO (z de eficiencia global, %d rep; color saturado en +-10)\n%s\n(ROJO = los hubs danan la red MAS que genes cualquiera de su modulo | AZUL = menos | * = peor que las %d repeticiones | NA = se quita el modulo completo)",
          N_REP, sub_t, N_REP), lim = 10, centro = 0)

rz <- mk("Razon_Dano")
etq_r <- matrix(ifelse(is.na(rz), "NA", sprintf("%.1f", rz)), nrow(rz))
graficar(log2(rz), etq_r, file.path(DIRECTORIO_SALIDA, "Heatmap_Razon_Dano_Hubs_vs_MismoModulo.pdf"),
  sprintf("Razon de dano: perdida de eficiencia por quitar los hubs / perdida por quitar genes al azar del mismo modulo (%d rep)\n%s\n(numero = veces mas dano que el azar; 1 = igual | color en log2, saturado entre 1/8x y 8x)",
          N_REP, sub_t), lim = 3, centro = 0)
# el CSV que escribe graficar() queda en log2; se reescribe con la razon real en el mismo orden
f_r <- file.path(DIRECTORIO_SALIDA, "Heatmap_Razon_Dano_Hubs_vs_MismoModulo.csv")
orden_r <- sub(" \\(.*", "", rownames(read.csv(f_r, row.names = 1, check.names = FALSE)))
rz_out <- rz[orden_r, , drop = FALSE]
rownames(rz_out) <- sprintf("%s (%d)", orden_r, tam[tolower(orden_r)])
write.csv(rz_out, f_r)

# ---- Resumen ----
H <- grepl(" H", columnas)
resumen <- data.frame(Modulo = mods, Tamano = as.integer(tam[tolower(mods)]),
  Z_Global_Azar_1a15H = round(rowMeans(mk("Z_Eficiencia")[, H, drop = FALSE], na.rm = TRUE), 2),
  Z_MismoModulo_1a15H = round(rowMeans(z[, H, drop = FALSE], na.rm = TRUE), 2),
  Razon_Dano_mediana_1a15H = round(apply(rz[, H, drop = FALSE], 1, median, na.rm = TRUE), 2),
  Celdas_signif_MismoModulo = rowSums(p <= 1 / (1 + N_REP), na.rm = TRUE))
resumen <- resumen[order(-resumen$Z_MismoModulo_1a15H), ]
write.csv(resumen, file.path(DIRECTORIO_SALIDA, "Ranking_Control_MismoModulo.csv"), row.names = FALSE)
print(resumen, row.names = FALSE)
cat("LISTO\n")
