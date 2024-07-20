library(tidyverse)
library(mvtnorm)
library(Cairo)

source("estimation.R")

density_model <- function(params, x) {
	with(params,
	sapply(1:length(p), function(i) p[i]*pnorm(x, mean=params$mu[2,][i], sd=sqrt(sigma[,,i][2,2])) ) %>%
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
	CairoPNG(paste0(DATASETS,"/images/CDF_Y/CDF_Y","-", n_clusters,"-",age,"-",ifelse(male, "male", "female")))
	print(
	      plot(ecdf(data$Y), 
		col="red", 
		main = "Fonction de Répartition", 
	 	xlab = "y", ylab = "P(Y < y)", lwd=8,
		sub=paste0("Tranche d'Àge : ",age, "-", age+4, "; Sexe : ", ifelse(male, "Homme", "Femme"), "; Clusters : ", n_clusters)))

	Y_min <- with(data, min(Y))
	Y_max <- with(data, max(Y))

	print(lines( seq(Y_min, Y_max, h_step), sapply(seq(Y_min, Y_max, h_step), function(x) density_model(params, x)), col="blue", type="l", lty=2, lwd=3))
	dev.off()
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

age <- 90
male <- FALSE

plot_CDF(DATASETS, 
	 data = filter_age_sex(df = df, age = age, male=male),
	 age = age,
	 male = male,
	 n_clusters = 2, 
	 h_step = 0.01
	)

plot_CDF(DATASETS, 
	 data = filter_age_sex(df = df, age = age, male=male),
	 age = age,
	 male = male,
	 n_clusters = 3, 
	 h_step = 0.01
	)

plot_CDF(DATASETS, 
	 data = filter_age_sex(df = df, age = age, male=male),
	 age = age,
	 male = male,
	 n_clusters = 4, 
	 h_step = 0.01
	)

plot_CDF(DATASETS, 
	 data = filter_age_sex(df = df, age = age, male=male),
	 age = age,
	 male = male,
	 n_clusters = 5, 
	 h_step = 0.01
	)

plot_CDF(DATASETS, 
	 data = filter_age_sex(df = df, age = age, male=male),
	 age = age,
	 male = male,
	 n_clusters = 6, 
	 h_step = 0.01
	)




