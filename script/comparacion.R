library(data.table)
library(dtwclust)
library(magrittr)
library(ggplot2)

# Probar las familias de ar_1, ar_2, ...ar_p con dtw
# Dejar la prueba de 0.5 y -0.5
# Dejar la prueba de 0.3, 0.8 y -0.3, -0.8

# Armar la grilla para cualquier caso.
# 

# Comparar familias ar de memoria corta y larga.

# ¿Por qué el jerárquico discrimina bien en positivos y negativos? Por la func.corr?
#


# Semilla
set.seed(123)

# Simulamos procesos autoregresivos
# N = 100
# phi = c(0.8)
# ar <- list(ar_1_pos = arima.sim(model = list(ar = phi), n = N),
#            ar_1_neg = arima.sim(model = list(ar = -phi), n = N))

generar_series <- function(phi, theta = NULL, N) {
    nombre <- ""
    if(!is.null(phi)) {
        signo_ar <- "p_"
        if(phi < 0) signo_ar <- "n_"
        nombre <- paste(nombre, "ar_", signo_ar, phi, sep = "")
    }
    if(!is.null(theta)) {
        signo_ma <- "p_"
        if(theta < 0) signo_ma <- "n_" 
        nombre <- paste(nombre, "ma_", signo_ma, theta, sep = "")
    }
    nombre <- gsub("\\.|-", "", nombre)
    arma <- data.table(modelo = nombre,
                       valores = as.vector(arima.sim(model = list(ar = phi, ma = theta), n = N)),
                       signo_ar = signo_ar,
                       signo_ma = signo_ma,
                       indices = 1:N)
    # setnames(arma, new = nombre)
    return(arma)
}

# generar_series(phi = phi, theta = 0.2, N = N)



# Parámetros
N = 100
phi = rep(c(0.3, 0.8), each = N)
# phi = rep(c(0.5), each = N)
theta = 0
k = 4
grilla <- setDT(expand.grid(N = N, phi = c(phi, -phi), theta = theta))

series_dt <- 
lapply(1:nrow(grilla), function(i) {
    N <- grilla[i, N]
    phi <- grilla[i, phi]
    theta <- grilla[i, theta]
    generar_series(phi = phi, theta = theta, N = N)
}) %>% 
    do.call("rbind", .)

series_dt$modelo <- paste(series_dt$modelo, rep(1:100, each = N), sep = "")
split(series_dt, by = "modelo", keep.by = FALSE) %>% 
    lapply(., function(x) x[["valores"]])

# series_dt[, modelo := paste(modelo, 1:.N, sep =)]

# Series simuladas:
ggplot(series_dt, aes(x = indices, y = valores, group = modelo, color = signo_ar)) + 
    geom_line(alpha = 0.1)

ggplot(series_dt, aes(x = indices, y = valores)) + 
    geom_line()

particional <- tsclust(
    series_list,
    type = "partitional",
    k = k,
    distance = "dtw_basic",
    centroid = "pam",
    seed = 3247L,
    trace = TRUE,
    args = tsclust_args(dist = list(window.size = 20L))
)

# Jerárquico
jerarquico <- tsclust(
    series_list,
    type = "hierarchical",
    k = k,
    distance = "sbd",
    trace = TRUE,
    control = hierarchical_control(method = "average")
)

# Agregar clusters
series_dt$clusters_dtw <- rep(particional@cluster, each = 2*N)
series_dt$clusters_jerarquico <- rep(jerarquico@cluster, each = 2*N)

# Labels
cluster_names <- c(
    `1` = "Cluster 1",
    `2` = "Cluster 2"
)

cluster_names <- c(
    `1` = "Cluster 1",
    `2` = "Cluster 2",
    `3` = "Cluster 3",
    `4` = "Cluster 4"
)

# plot(particional)

series_dt[grep("ar_p_.*", modelo), phi := gsub("ar_._(\\d+).*", "\\1", modelo)]
series_dt[!grep("ar_p_.*", modelo), phi := gsub("ar_._(\\d+).*", "-\\1", modelo)]

# Como quedan?
p1 <- ggplot(series_dt, aes(x = indices, y = valores, group = modelo, color = signo_ar)) + 
    geom_line(alpha = 0.1) +
    facet_wrap(clusters_dtw ~ ., labeller = as_labeller(cluster_names)) +
    labs(x = "Observaciones", y = "Valores simulados", title = "DTW")

# Opa!
p2 <- ggplot(series_dt, aes(x = indices, y = valores, group = modelo, color = signo_ar)) + 
    geom_line(alpha = 0.1) +
    facet_wrap(clusters_jerarquico ~ ., labeller = as_labeller(cluster_names)) +
    labs(x = "Observaciones", y = "Valores simulados", title = "Jerárquico")

# Color el phi utilizado
p1 <- ggplot(series_dt, aes(x = indices, y = valores, group = modelo, color = phi)) + 
    geom_line(alpha = 0.2) +
    facet_wrap(clusters_dtw ~ ., labeller = as_labeller(cluster_names)) +
    labs(x = "Observaciones", y = "Valores simulados", title = "DTW") +
    scale_color_discrete()

# Opa!
p2 <- ggplot(series_dt, aes(x = indices, y = valores, group = modelo, color = phi)) + 
    geom_line(alpha = 0.2) +
    facet_wrap(clusters_jerarquico ~ ., labeller = as_labeller(cluster_names)) +
    labs(x = "Observaciones", y = "Valores simulados", title = "Jerárquico") +
    scale_color_discrete()


library(patchwork)
p1+p2

# Cantidad de curvas por grupo dependiendo el signo
series_dt[, by = .(signo_ar, clusters_jerarquico), .N/N]
series_dt[, by = .(signo_ar, clusters_dtw), .N/N]

# Medidas de comparación entre clusters
particional@clusinfo
centroides <- particional@centroids
particional@cldist

particional

# Evaluación
k0 <- 2
k1 <- 10 #floor(nrow(grilla)/2)
data <- series_list
pc_k <- tsclust(data, k = k0:k1,
                distance = "dtw_basic", centroid = "pam",
                seed = 94L)
names(pc_k) <- paste0("k_", k0:k1)
sapply(pc_k, cvi, type = "internal")
    
indices <- as.data.table(sapply(pc_k, cvi, type = "valid"), keep.rownames = TRUE)
indices <- melt.data.table(data = indices, id.vars = "rn")    

indices[, variable := as.integer(gsub("k_", "", variable))]


ggplot(indices, aes(x = variable, y = value, color = rn, group = rn)) +
    geom_line() +
    geom_point() +
    facet_wrap(rn ~ ., scales = "free_y")
# Revisar el SF no tiene sentido.


# Comparación general -----------------------------------------------------

cfg <- compare_clusterings_configs(
    types = "partitional",
    k = k,
    controls = list(
        partitional = partitional_control(
            iter.max = 20L
        )
    ),
    distances = pdc_configs(
        "distance",
        partitional = list(
            dtw_basic = list(
                window.size = seq(from = 10L, to = 30L, by = 5L),
                norm = c("L1", "L2")
            )
        )
    ),
    centroids = pdc_configs(
        "centroid",
        share.config = c("p"),
        dba = list(
            window.size = seq(from = 10L, to = 30L, by = 5L),
            norm = c("L1", "L2")
        )
    ),
    no.expand = c(
        "window.size",
        "norm"
    )
)


evaluators <- cvi_evaluators("ARI", ground.truth = CharTrajLabels)
comparison <- compare_clusterings(series_list, types = "partitional",
                                  configs = cfg, seed = 8L,
                                  score.clus = evaluators$score,
                                  pick.clus = evaluators$pick)    
    


head(comparison$results$partitional[, c("distance",
                                        "centroid",
                                        "window.size_distance",
                                        "norm_distance",
                                        "ARI")])
