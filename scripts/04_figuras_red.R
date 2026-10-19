# ==============================================================================
# FIGURAS DE LA RED (red de 24 modulos)
#   Fig1_Dendrograma_Modulos.pdf : arbol jerarquico (1 - TOM) con los colores de modulo
#   Fig2_Red_Genes.pdf / .png    : grafo de genes (aristas TOM >= umbral top 1%, sin grey),
#                                  coloreado por modulo; hub principal de cada modulo resaltado
#   Fig3_Red_Modulos.pdf         : red de modulos; tamano = n genes, grosor de arista =
#                                  densidad de aristas fuertes entre los dos modulos
# Usa el mismo grafo que 01_robustez_global_integridad.R (umbral = top 1% de la TOM original).
# ==============================================================================
.libPaths(c("/home/jmarquez/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
suppressMessages({ library(WGCNA); library(igraph); library(AnnotationDbi); library(org.Hs.eg.db) })
enableWGCNAThreads(as.integer(Sys.getenv("NUCLEOS", "16")))
setwd("/STORAGE/genut/Manuel/directed_attack_network")
load("data/datExpr_colores_24Modulos.RData")
POWER_OPTIMO <- 12; TYPE_NETWORK <- "signed"; FRAC_ARISTAS <- 0.01
DIR_FIG <- "figures"; DIR_RES <- "results/red"
dir.create(DIR_FIG, showWarnings = FALSE); dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

colores_base <- colores_base[colnames(datExpr_5000)]
genes <- colnames(datExpr_5000)
adj <- adjacency(datExpr_5000, power = POWER_OPTIMO, type = TYPE_NETWORK)
Alldegrees <- intramodularConnectivity(adj, colores_base)
TOM <- TOMsimilarity(adj, TOMType = TYPE_NETWORK, verbose = 0)
dimnames(TOM) <- list(genes, genes)
v <- TOM[upper.tri(TOM)]; k <- round(length(v) * FRAC_ARISTAS)
UMBRAL <- sort(v, partial = length(v) - k + 1)[length(v) - k + 1]; rm(v)
modulos <- sort(setdiff(unique(colores_base), "grey"))
tam <- table(colores_base)

# ---- Hubs (top 5 por modulo) con simbolo ----
hubs <- do.call(rbind, lapply(modulos, function(m) {
  df <- Alldegrees[colores_base == m, ]
  df <- df[order(df$kWithin, decreasing = TRUE), ][1:min(5, nrow(df)), ]
  data.frame(Modulo = toupper(m), Rango = seq_len(nrow(df)), Ensembl = rownames(df),
             kWithin = round(df$kWithin, 2), kTotal = round(df$kTotal, 2))
}))
simb <- suppressMessages(mapIds(org.Hs.eg.db, keys = sub("\\..*", "", hubs$Ensembl),
                                column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
hubs$Simbolo <- ifelse(is.na(simb), hubs$Ensembl, simb)
write.csv(hubs, file.path(DIR_RES, "Hubs_Top5_Por_Modulo.csv"), row.names = FALSE)
cat(sprintf("Umbral TOM top 1%%: %.4f | hubs anotados: %d/%d\n", UMBRAL, sum(!is.na(simb)), nrow(hubs)))

# ---- Fig 1: dendrograma ----
arbol <- hclust(as.dist(1 - TOM), method = "average")
pdf(file.path(DIR_FIG, "Fig1_Dendrograma_Modulos.pdf"), width = 12, height = 6)
plotDendroAndColors(arbol, colores_base, "Modulo", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = sprintf("Dendrograma de genes (1 - TOM, power %d, signed) y asignacion a modulos: %d modulos + grey (%d genes sin modulo)",
                                   POWER_OPTIMO, length(modulos), tam["grey"]))
invisible(dev.off())

# ---- Fig 2: red de genes ----
gm <- genes[colores_base != "grey"]
A <- (TOM[gm, gm] >= UMBRAL) * 1; diag(A) <- 0
g <- graph_from_adjacency_matrix(A, mode = "undirected")
g <- delete_vertices(g, which(degree(g) == 0))
# componentes muy chicos (< 10 genes) quedan lejos en el layout y comprimen la figura
cmp <- igraph::components(g)
n_fuera <- sum(cmp$csize[cmp$csize < 10])
g <- delete_vertices(g, which(cmp$csize[cmp$membership] < 10))
V(g)$modulo <- colores_base[V(g)$name]
hub1 <- hubs$Ensembl[hubs$Rango == 1]
es_hub <- V(g)$name %in% hub1
lay <- layout_with_fr(g, niter = 1000)
etq <- ifelse(es_hub, hubs$Simbolo[match(V(g)$name, hubs$Ensembl)], NA)
dibujar_red_genes <- function() {
  par(mar = c(0, 0, 3, 0))
  plot(g, layout = lay, vertex.label = NA, vertex.size = ifelse(es_hub, 2.2, 0.8),
       vertex.color = V(g)$modulo, vertex.frame.color = ifelse(es_hub, "black", NA),
       edge.color = adjustcolor("grey40", alpha.f = 0.04), edge.width = 0.3,
       main = sprintf("Red de genes: %d genes en modulos, %d aristas con TOM >= %.4f (top 1%% de la red)\nColor = modulo | nodos grandes con borde = hub principal de cada modulo\n(se omiten %d genes en componentes de < 10 genes)",
                      vcount(g), ecount(g), UMBRAL, n_fuera))
  idx <- which(es_hub)
  lay_n <- norm_coords(lay)    # plot.igraph reescala el layout a [-1, 1]
  text(lay_n[idx, 1], lay_n[idx, 2], etq[idx], cex = 0.7, font = 2, pos = 3, offset = 0.4)
}
pdf(file.path(DIR_FIG, "Fig2_Red_Genes.pdf"), width = 15, height = 15); dibujar_red_genes(); invisible(dev.off())
png(file.path(DIR_FIG, "Fig2_Red_Genes.png"), width = 4200, height = 4200, res = 300, type = "cairo"); dibujar_red_genes(); invisible(dev.off())

# ---- Fig 3: red de modulos ----
D <- matrix(0, length(modulos), length(modulos), dimnames = list(modulos, modulos))
for (a in modulos) for (b in modulos) if (a < b) {
  ga <- genes[colores_base == a]; gb <- genes[colores_base == b]
  D[a, b] <- D[b, a] <- sum(TOM[ga, gb] >= UMBRAL) / (length(ga) * length(gb))
}
dens_int <- sapply(modulos, function(m) { gx <- genes[colores_base == m]; s <- TOM[gx, gx] >= UMBRAL
  (sum(s) - length(gx)) / (length(gx) * (length(gx) - 1)) })
write.csv(round(D, 5), file.path(DIR_RES, "Densidad_Aristas_Entre_Modulos.csv"))
write.csv(data.frame(Modulo = toupper(modulos), Tamano = as.integer(tam[modulos]), Densidad_Interna = round(dens_int, 4)),
          file.path(DIR_RES, "Densidad_Aristas_Interna.csv"), row.names = FALSE)
gmod <- graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE, diag = FALSE)
gmod <- delete_edges(gmod, which(E(gmod)$weight < quantile(E(gmod)$weight[E(gmod)$weight > 0], 0.5)))
lay_m <- layout_with_fr(gmod, weights = E(gmod)$weight / max(E(gmod)$weight))
pdf(file.path(DIR_FIG, "Fig3_Red_Modulos.pdf"), width = 11, height = 11)
par(mar = c(0, 0, 4, 0))
plot(gmod, layout = lay_m, vertex.color = modulos, vertex.size = 4 + 2.2 * sqrt(as.numeric(tam[modulos]) / 10),
     vertex.label = sprintf("%s\n(%d)", toupper(modulos), as.integer(tam[modulos])), vertex.label.cex = 0.7,
     vertex.label.color = "black", vertex.label.family = "sans",
     vertex.label.dist = 1.2 + 0.35 * sqrt(as.numeric(tam[modulos]) / 10), vertex.label.degree = -pi / 2, vertex.frame.color = "grey30",
     edge.width = 0.5 + 12 * E(gmod)$weight / max(E(gmod)$weight), edge.color = adjustcolor("grey30", 0.5),
     main = "Red de modulos: tamano del nodo = n genes; grosor = densidad de aristas fuertes (TOM top 1%) entre modulos\n(se muestran las conexiones por encima de la mediana)")
invisible(dev.off())
cat("LISTO\n")
