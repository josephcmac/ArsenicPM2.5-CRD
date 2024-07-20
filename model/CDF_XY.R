library(tidyverse)
library(mvtnorm)
library(Cairo)
library(RColorBrewer)

source("estimation.R")

filter_age_sex <- function(df, age, male) {
	if (male) {
		df %>% 
			filter(age_min == age, sex == "Homme") %>%
			select(X, Y)
	} else {
		df %>% 
			filter(age_min == age, sex == "Femme") %>%
			select(X, Y)
	}
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

produce_image <- function(n_clusters, age, male) {
    data <- filter_age_sex(df, age = age, male = male)

    X_min <- with(data, min(X))
    X_max <- with(data, max(X) + 0.01)
    Y_min <- with(data, min(Y))
    Y_max <- with(data, max(Y) + 0.01)

    X_len <- 50
    Y_len <- 50

    X_bin <- seq(X_min, X_max, length.out = X_len + 1)
    Y_bin <- seq(Y_min, Y_max, length.out = Y_len + 1)

    Observed <- sapply(1:X_len, function(i) {
        sapply(1:Y_len, function(j) {
            with(data, sum((X_bin[i] <= X) & (X < X_bin[i + 1]) & (Y_bin[j] <= Y) & (Y < Y_bin[j + 1])))
        })
    })

    params <- compute_params_model(data, n_clusters)
    Expected <- nrow(data) * sapply(1:X_len, function(i) {
        sapply(1:Y_len, function(j) {
            with(data, with(
                params,
                sapply(1:n_clusters, function(k) {
                    (Y_bin[j + 1] - Y_bin[j]) * (X_bin[i + 1] - X_bin[i]) * p[k] * mvtnorm::dmvnorm(
                        x = c(0.5 * (X_bin[i] + X_bin[i + 1]), 0.5 * (Y_bin[j] + Y_bin[j + 1])),
                        mean = mu[, k], sigma = sigma[, , k]
                    )
                }) %>% sum()
            ))
        })
    })

    plot_combined_matrices <- function(observed_matrix, expected_matrix, file_name, title1, title2, n_clusters) {
        CairoPNG(file_name, width = 1400, height = 700)
        layout(matrix(1:2, nrow = 1, ncol = 2))

        par(mar = c(5, 4, 4, 1) + 0.1)
        image(
            observed_matrix,
            col = gray.colors(100, start = 0, end = 1),
            axes = FALSE,
            main = title1
        )
        box()
        axis(1, at = seq(0, 1, length.out = ncol(observed_matrix)), labels = FALSE)
        axis(2, at = seq(0, 1, length.out = nrow(observed_matrix)), labels = FALSE)

        par(mar = c(5, 1, 4, 4) + 0.1)
        image(
            expected_matrix,
            col = gray.colors(100, start = 0, end = 1),
            axes = FALSE,
            main = title2
        )
        box()
        axis(1, at = seq(0, 1, length.out = ncol(expected_matrix)), labels = FALSE)
        axis(2, at = seq(0, 1, length.out = nrow(expected_matrix)), labels = FALSE)

       dev.off()
    }

    combined_file <- paste0(DATASETS, "/images/CDF_XY/", ifelse(male, "male", "female"),"/combined_CDF_XY-", n_clusters, "-", age, "-", ifelse(male, "male", "female"), ".png")

    plot_combined_matrices(Observed, Expected, combined_file, "Observed CDF XY", "Expected CDF XY", n_clusters)
}


for (n_clusters in 2:6) {
	for (male in c(TRUE, FALSE))
  produce_image(n_clusters=n_clusters, age=90, male=male)
}



