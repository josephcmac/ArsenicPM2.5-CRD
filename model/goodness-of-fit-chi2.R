library(tidyverse)
library(mvtnorm)

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


data <- filter_age_sex(df, age = 90, male = FALSE)

X_min <- with(data, min(X))
X_max <- with(data, max(X) + 0.01)
Y_min <- with(data, min(Y))
Y_max <- with(data, max(Y) + 0.01)

X_len <- 50
Y_len <- 50

X_bin <- seq(X_min, X_max, length.out = X_len + 1)
Y_bin <- seq(Y_min, Y_max, length.out = Y_len + 1)


X2 <- sapply(2:25, function(n_clusters) {
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



   sum((Observed - Expected)**2/Expected)
})

data.frame(n_clusters=2:25, chi2=X2) %>%
	ggplot(aes(n_clusters, chi2)) +
	geom_point() +
	geom_line() +
	scale_y_log10() +
	theme_classic()




