# ==============================================================================
# ROBUSTEZ GLOBAL: INTEGRIDAD DE LA RED ANTE ATAQUES DIRIGIDOS (red de 24 modulos)
# Grafo = pares de genes con TOM >= umbral (umbral fijo = top 1% de la red original).
# Tras cada ataque se recalcula la TOM sin los hubs y se aplica el MISMO umbral
# (asi se capturan aristas perdidas y aristas que se debilitan/fortalecen).
# Metricas (todas relativas a la red original, normalizadas por N original):
#   - Eficiencia_Rel = E_global(ataque) / E_global(original)  (caminos cortos entre genes)
#   - LCC_Rel        = tamano componente gigante (ataque) / (original)
#   - N_Componentes  = componentes con >1 gen
# Referencia: quitar el MISMO numero de genes AL AZAR (N_REP repeticiones).
#   - Z_Eficiencia = (media_azar - ataque) / sd_azar   (>0 = el ataque dana mas que el azar)
# ==============================================================================
.libPaths(c("/home/jmarquez/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
suppressMessages({ library(WGCNA); library(igraph); library(pheatmap); library(parallel) })
disableWGCNAThreads()
setwd("/STORAGE/genut/Manuel/directed_attack_network")
load("data/datExpr_colores_24Modulos.RData")

POWER_OPTIMO <- 12; TYPE_NETWORK <- "signed"
FRAC_ARISTAS <- 0.01
N_REP   <- as.integer(Sys.getenv("N_REP", "30"))
NUCLEOS <- as.integer(Sys.getenv("NUCLEOS", "16"))
PRUEBA  <- Sys.getenv("PRUEBA", "0") == "1"
DIRECTORIO_SALIDA <- if (PRUEBA) "results/prueba_integridad" else "results/integridad"
dir.create(DIRECTORIO_SALIDA, showWarnings = FALSE)
set.seed(2026)

genes <- colnames(datExpr_5000); N <- length(genes)
adj_base <- adjacency(datExpr_5000, power = POWER_OPTIMO, type = TYPE_NETWORK)
Alldegrees <- intramodularConnectivity(adj_base, colores_base)
TOM_base <- TOMsimilarity(adj_base, TOMType = TYPE_NETWORK, verbose = 0)
v <- TOM_base[upper.tri(TOM_base)]
UMBRAL <- sort(v, partial = length(v) - round(length(v) * FRAC_ARISTAS) + 1)[length(v) - round(length(v) * FRAC_ARISTAS) + 1]
rm(v)

metricas <- function(TOM) {
  n <- ncol(TOM)
  g <- graph_from_adjacency_matrix((TOM >= UMBRAL) * 1, mode = "undirected", diag = FALSE)
  cmp <- components(g)
  c(Eficiencia = global_efficiency(g) * n * (n - 1) / (N * (N - 1)),   # normalizada por N original
    LCC = max(cmp$csize), N_Componentes = sum(cmp$csize > 1), Aristas = gsize(g))
}
m_orig <- metricas(TOM_base); rm(TOM_base); gc()
cat(sprintf("Red original: umbral TOM=%.4f, E=%.4f, LCC=%d, componentes=%d, aristas=%d\n",
            UMBRAL, m_orig["Eficiencia"], m_orig["LCC"], m_orig["N_Componentes"], m_orig["Aristas"]))

red_sin <- function(rem) {
  keep <- which(!(genes %in% rem))
  metricas(TOMsimilarity(adj_base[keep, keep], TOMType = TYPE_NETWORK, verbose = 0))
}

modulos <- sort(setdiff(unique(colores_base), "grey"))
if (PRUEBA) modulos <- c("turquoise", "salmon")
tamanos_modulos <- table(colores_base)
columnas <- c(paste0(1:15, " H"), "5%", "10%", "20%", "90%")
hubs_de <- function(mod, col) {
  g <- names(colores_base[colores_base == mod])
  df <- Alldegrees[rownames(Alldegrees) %in% g, ]
  df <- df[order(df$kWithin, decreasing = TRUE), ]
  n <- if (grepl("%", col)) max(1, round(length(g) * as.numeric(sub("%", "", col)) / 100))
       else as.integer(sub(" H", "", col))
  rownames(df)[1:min(n, nrow(df))]
}
if (PRUEBA) columnas <- c("1 H", "15 H", "90%")

# ---- Trabajos: ataques dirigidos + controles al azar ----
ataques <- expand.grid(col = columnas, mod = modulos, stringsAsFactors = FALSE)
ataques$rem <- mapply(hubs_de, ataques$mod, ataques$col, SIMPLIFY = FALSE)
ataques$n <- lengths(ataques$rem)
tamanos_n <- sort(unique(ataques$n))
azar <- expand.grid(rep = seq_len(N_REP), n = tamanos_n)
azar$rem <- lapply(azar$n, function(n) sample(genes, n))
trabajos <- c(ataques$rem, azar$rem)
cat(sprintf("Ataques: %d | Controles al azar: %d (%d tamanos x %d rep)\n",
            nrow(ataques), nrow(azar), length(tamanos_n), N_REP))

t0 <- Sys.time()
res <- mclapply(trabajos, red_sin, mc.cores = NUCLEOS, mc.preschedule = FALSE)
fallos <- sapply(res, function(r) inherits(r, "try-error") || is.null(r))
if (any(fallos)) stop(sprintf("%d trabajos fallaron", sum(fallos)))
res <- do.call(rbind, res)
cat(sprintf("Tiempo: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

A <- cbind(ataques[, c("mod", "col", "n")], res[seq_len(nrow(ataques)), ])
R <- cbind(azar[, c("rep", "n")], res[-seq_len(nrow(ataques)), ])
A$Eficiencia_Rel <- A$Eficiencia / m_orig["Eficiencia"]
A$LCC_Rel <- A$LCC / m_orig["LCC"]
R$Eficiencia_Rel <- R$Eficiencia / m_orig["Eficiencia"]
R$LCC_Rel <- R$LCC / m_orig["LCC"]
ref <- aggregate(cbind(Eficiencia_Rel, LCC_Rel, N_Componentes) ~ n, R, function(x) c(media = mean(x), sd = sd(x)))
ref <- do.call(data.frame, ref)
A <- merge(A, ref, by = "n", suffixes = c("", "_azar"))
A$Z_Eficiencia <- (A$Eficiencia_Rel.media - A$Eficiencia_Rel) / A$Eficiencia_Rel.sd
A$Z_LCC        <- (A$LCC_Rel.media - A$LCC_Rel) / A$LCC_Rel.sd
A$P_Eficiencia <- sapply(seq_len(nrow(A)), function(i)        # p empirico: fraccion de azar con dano >= ataque
  (1 + sum(R$Eficiencia_Rel[R$n == A$n[i]] <= A$Eficiencia_Rel[i])) / (1 + N_REP))
A$Modulo <- toupper(A$mod); A$Ataque <- A$col
write.csv(A[, setdiff(names(A), c("mod", "col"))], file.path(DIRECTORIO_SALIDA, "Resultados_Integridad_Largo.csv"), row.names = FALSE)
write.csv(R, file.path(DIRECTORIO_SALIDA, "Controles_Azar.csv"), row.names = FALSE)
write.csv(data.frame(Umbral_TOM = UMBRAL, t(m_orig)), file.path(DIRECTORIO_SALIDA, "Red_Original_Metricas.csv"), row.names = FALSE)

# ---- Heatmaps ----
a_matriz <- function(var) {
  m <- matrix(NA, length(modulos), length(columnas), dimnames = list(toupper(modulos), columnas))
  for (r in seq_len(nrow(A))) m[A$Modulo[r], A$Ataque[r]] <- A[[var]][r]
  m
}
paleta <- colorRampPalette(c("#8B0000", "#CC0000", "#FF4500", "#FF8C00", "#FFA500",
  "#FFD700", "#FFFF00", "#ADFF2F", "#7CFC00", "#32CD32", "#00AA00", "#006400"))(100)
etiquetar <- function(m) sprintf("%s (%d)", rownames(m), tamanos_modulos[tolower(rownames(m))])

graficar <- function(m, archivo, titulo, formato, color, breaks = NA, orden_desc = FALSE, etiquetas = NULL) {
  o <- order(rowMeans(m), decreasing = orden_desc)
  m <- m[o, , drop = FALSE]
  if (!is.null(etiquetas)) etiquetas <- etiquetas[o, , drop = FALSE]
  rownames(m) <- etiquetar(m)
  write.csv(m, sub("\\.pdf$", ".csv", archivo))
  nums <- if (is.null(etiquetas)) TRUE else etiquetas
  pdf(archivo, width = max(16, (ncol(m) * 40 + 350) / 72), height = max(6, (nrow(m) * 30 + 220) / 72))
  tryCatch(pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = nums,
           number_format = formato, number_color = "black", fontsize_number = 7,
           color = color, breaks = breaks, main = titulo, angle_col = 45,
           cellwidth = 40, cellheight = 30, fontsize = 10, border_color = "grey60"),
           finally = dev.off())
}
sub_titulo <- sprintf("Grafo: aristas TOM >= %.4f (top 1%% de la red original; E=%.3f, componente gigante=%d genes)",
                      UMBRAL, m_orig["Eficiencia"], m_orig["LCC"])

graficar(a_matriz("Eficiencia_Rel"), file.path(DIRECTORIO_SALIDA, "Heatmap_Eficiencia_Global_Relativa.pdf"),
  paste0("Robustez global: eficiencia global de la red atacada / red original\n", sub_titulo,
         "\n(1 = la red sigue igual de comunicada | ROJO = la red pierde comunicacion entre genes)"),
  "%.3f", paleta)
graficar(a_matriz("LCC_Rel"), file.path(DIRECTORIO_SALIDA, "Heatmap_Componente_Gigante_Relativo.pdf"),
  paste0("Robustez global: tamano del componente gigante de la red atacada / red original\n", sub_titulo,
         "\n(1 = la red sigue unida | ROJO = la red se fragmenta)"),
  "%.3f", paleta)

z <- a_matriz("Z_Eficiencia"); p <- a_matriz("P_Eficiencia")
lim <- 10   # color saturado en +-10 para que no lo domine un solo modulo; los numeros son los reales
etq <- matrix(sprintf("%.1f%s", z, ifelse(p <= 1 / (1 + N_REP), "*", "")), nrow(z), dimnames = dimnames(z))
write.csv(z, file.path(DIRECTORIO_SALIDA, "Matriz_Z_Eficiencia.csv"))
graficar(pmin(pmax(z, -lim), lim), file.path(DIRECTORIO_SALIDA, "Heatmap_Dano_vs_Azar_Eficiencia.pdf"),
  paste0("Dano del ataque dirigido vs. quitar el MISMO numero de genes al azar (z de eficiencia global, ", N_REP, " repeticiones; color saturado en +-10)\n",
         sub_titulo, "\n(ROJO = el ataque dana la red MAS que el azar | AZUL = menos que el azar | * = peor que las ", N_REP, " repeticiones al azar)"),
  "%.1f", colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"))(100),
  breaks = seq(-lim, lim, length.out = 101), orden_desc = TRUE, etiquetas = etq)

# ---- Resumen por modulo ----
res_mod <- data.frame(Modulo = toupper(modulos), Tamano = as.integer(tamanos_modulos[modulos]),
  Eficiencia_15H = round(a_matriz("Eficiencia_Rel")[, if ("15 H" %in% columnas) "15 H" else ncol(z)], 3),
  Z_medio_1a15H = round(rowMeans(z[, grepl(" H", columnas), drop = FALSE]), 2),
  Z_medio_pct = round(rowMeans(z[, grepl("%", columnas), drop = FALSE]), 2),
  Celdas_signif = rowSums(p <= 1 / (1 + N_REP)))
res_mod <- res_mod[order(-res_mod$Z_medio_1a15H), ]
write.csv(res_mod, file.path(DIRECTORIO_SALIDA, "Ranking_Robustez_Global.csv"), row.names = FALSE)
print(res_mod, row.names = FALSE)
cat("LISTO\n")
