library(tidyverse)
library(mvtnorm)
library(gplm)

source("estimation.R")

filter_age_sex <- function(df, age, male) {
	if (male) {
		df %>% 
			filter(age_min == age) %>%
			select(ArsenicPM2.5LC_log, male_logit) %>%
			rename(X = ArsenicPM2.5LC_log, Y = male_logit)
	} else {
		df %>% 
			filter(age_min == age) %>%
			select(ArsenicPM2.5LC_log, female_logit) %>%
			rename(X = ArsenicPM2.5LC_log, Y = female_logit)
	}
}


calc_esperance_cond <- function(X, params, n_clusters) {
	pretty_params <- function(params, n_clusters) {
		pretty_params_aux <- function(params, n_clusters) {
			with(params,
			     list(
				  p = p,
				  mu = mu,
				  sigma = sigma,
				  moys_x = mu[1, ],
				  sds_x = sqrt(sigma[1, 1, 1:n_clusters]),
				  moys_y = mu[2, ],
				  sds_y = sqrt(sigma[2, 2, 1:n_clusters])
		     		)
			)
		}

		with(pretty_params_aux(params, n_clusters),
		     with(params,
		     list(
			  p = p,
			  mu = mu,
			  sigma = sigma,
			  moys_x = moys_x,
			  sds_x = sds_x,
			  moys_y = moys_y,
			  sds_y = sds_y,
			  rho_xy = sigma[1, 2, 1:n_clusters] / (sds_x * sds_y)
			  )
		     )
		)
	}


	calc_cond_probs <- function(X, moys_x, sds_x, p) {
		calc_cond_probs_aux <- function(A) {
			A / rowSums(A)
		}
		sapply(1:n_clusters, function(i) p[i] * dnorm(X, mean = moys_x[i], sd = sds_x[i])) %>% 
			calc_cond_probs_aux
		
	}

	calc_cond_means <- function(X, moys_x, sds_x, moys_y, sds_y, rho_xy) {
		sapply(1:n_clusters, function(i) moys_y[i] + rho_xy[i] * sds_y[i] * (X - moys_x[i]) / sds_x[i])
	}

	with(pretty_params(params, n_clusters), 
	     
		  rowSums(calc_cond_probs(X, moys_x, sds_x, p) * calc_cond_means(X, moys_x, sds_x, moys_y, sds_y, rho_xy))
	)
}



DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

# Number of clusters
n_clusters <- 4

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))
Y_min <- with(df, min(c(male_logit, female_logit)))
Y_max <- with(df, max(c(male_logit, female_logit)))

data <- filter_age_sex(df, age=90, male=TRUE)

params <- compute_params_model(data, n_clusters) 

esperance_cond <- with(data, calc_esperance_cond(X, params, n_clusters))
kernel_regression <- with(data, kreg(X, Y))

with(data, plot(X, Y, type = "p", col = "black"))
with(data, points(X, esperance_cond, col = "red"))
with(kernel_regression, points(x, y, col = "green"))
title("Esperance conditionnelle de Y sachant X=x en fonction de x")



