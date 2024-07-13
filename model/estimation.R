library(flexclust)
library(mclust)

estimate_init_model <- function(data, n_clusters) {
  return(
    with(
      kmeans(
        x = data, 
	centers = kcca(data, k=n_clusters, family=kccaFamily("kmeans"))@centers, # k-means++
	iter.max = 10, nstart = 1,
        algorithm = "Hartigan-Wong", trace = FALSE
      ),
      list(
        means = centers,
        variances = lapply(1:n_clusters, function(i) cov(data[cluster == i, ]))
      )
    )
  )
}


compute_params_model <- function(data, n_clusters) {
  with(
    Mclust(data,
      G = n_clusters,
      modelNames = "VVV",
      initialization = estimate_init_model(data, n_clusters)
    ), with(
      parameters,
      list(
        p = pro,
        mu = mean,
        sigma = variance$sigma,
	loglik = loglik,
	bic = bic
      )
    )
  )
}


