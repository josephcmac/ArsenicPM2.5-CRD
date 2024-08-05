library(tidyverse)
library(mvtnorm)
library(Cairo)

source("estimation.R")

density_model <- function(params, x) {
	with(params,
	sapply(1:length(p), function(i) p[i]*pnorm(x, mean=params$mu[1,][i], sd=sqrt(sigma[,,i][1,1])) ) %>%
		sum
	)
}


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



plot_CDF <- function(DATASETS, data, age, male, n_clusters, h_step) {
	params <- compute_params_model(data, n_clusters)
	CairoPNG(paste0(DATASETS,"/images/CDF_X/",ifelse(male, "male", "female"),"/CDF_X","-", n_clusters,"-",age,"-",ifelse(male, "male", "female")))
	print(
	      plot(ecdf(data$X), 
		col="red", 
		main = "Fonction de Répartition", 
	 	xlab = "x", ylab = "P(X < x)", lwd=8,
		sub=paste0("Tranche d'Àge : ",age, "-", age+4, "; Sexe : ", ifelse(male, "Homme", "Femme"), "; Clusters : ", n_clusters)))

	X_min <- with(data, min(X))
	X_max <- with(data, max(X))

	print(lines( seq(X_min, X_max, h_step), sapply(seq(X_min, X_max, h_step), function(x) density_model(params, x)), col="blue", type="l", lty=2, lwd=3))
	dev.off()
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

age <- 90
male <- FALSE

for (age in seq(0, 90, 5)) {
  for (male in c(TRUE, FALSE)) {
	plot_CDF(DATASETS, 
		 data = filter_age_sex(df = df, age = age, male=male),
		 age = age,
		 male = male,
		 n_clusters = 6, 
		 h_step = 0.01
		)
  }
}




