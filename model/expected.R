library(tidyverse)
library(gplm)
library(Cairo)

set.seed(1729)

source("estimation.R")

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
		sapply(1:n_clusters, function(i) p[i] * dnorm(X, mean = moys_x[i], 
		                                              sd = sds_x[i])) %>% 
			calc_cond_probs_aux
		
	}

	calc_cond_means <- function(X, moys_x, sds_x, moys_y, sds_y, rho_xy) {
		sapply(1:n_clusters, function(i) 
		  moys_y[i] + rho_xy[i] * sds_y[i] * (X - moys_x[i]) / sds_x[i])
	}

	with(pretty_params(params, n_clusters),  
		  rowSums(calc_cond_probs(X, moys_x, sds_x, p) * 
		            calc_cond_means(X, moys_x, sds_x, moys_y, sds_y, rho_xy))
	)
}


plot_expected <- function(DATASETS, data, age, male, n_clusters, Y_min, Y_max) {
	params <- compute_params_model(data, n_clusters) 
	kernel_regression <- with(data, kreg(X, Y, bandwidth=0.3, kernel="triweight"))
	esperance_cond <- with(kernel_regression, 
	                       calc_esperance_cond(x, params, n_clusters))
	
	CairoPNG(paste0(DATASETS,"/images/expected/", 
	                ifelse(male, "male", "female"),"/expected_Y-",age,"-",
	                ifelse(male, "male", "female"),".png"))
	with(data, plot(X, Y, type = "p", col = rgb(0,0,0,.1), pch=19, 
	     main="Espérance conditionnelle de Y sachant X=x en fonction de x",
	     sub=paste0("Tranche d'Àge : ",age, "-", 
	                age+4, "; Sexe : ", ifelse(male, "Homme", "Femme")),
	     ylim=c(Y_min, Y_max)

	)
	)
	with(kernel_regression, 
	     lines(x, y, col = "blue", pch=19, lty="dashed", lwd=3))
	with(kernel_regression, 
	     lines(x, esperance_cond, col = "red", pch=19, lty="solid", lwd=3))
	dev.off()

}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))
dr <- readRDS(paste0(DATASETS, "/tables/clusters.rds"))

Y_min <- with(df, min(Y))
Y_max <- with(df, max(Y))

sapply(seq(0,90,5), function(age) sapply(c(TRUE,FALSE), function(male)
	plot_expected(DATASETS=DATASETS,
	      data=filter_age_sex(df, age=age, male=male), 
	      age=age, 
	      male=male, 
	      n_clusters=dr %>% 
	        filter(
	          age_min == age, 
	          sex == ifelse(male, "Homme", "Femme")
	        ) %>% 
	        getElement("clusters"), 
	      Y_min = Y_min,
	      Y_max = Y_max
)
					  ))



