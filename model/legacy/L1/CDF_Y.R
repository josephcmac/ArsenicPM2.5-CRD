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


plot_CDF <- function(DATASETS, data, age, male, n_clusters, h_step, Y_min, Y_max) {
	params <- compute_params_model(data, n_clusters)
	CairoPNG(paste0(DATASETS,"/images/CDF_Y/CDF_Y","-", n_clusters,"-",age,"-",ifelse(male, "male", "female")))
	print(
	      plot(ecdf(data$Y), 
		col="red", 
		xlim=c(Y_min, Y_max), 
		main = "Fonction de Répartition", 
	 	xlab = "y", ylab = "P(Y < y)", 
		sub=paste0("Tranche d'Àge : ",age, "-", age+4, "; Sexe : ", ifelse(male, "Homme", "Femme"))))

	print(lines( seq(Y_min, Y_max, h_step), sapply(seq(Y_min, Y_max, h_step), function(x) density_model(params, x)), col="blue", type="l", lty=2, lwd=3))
	dev.off()
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))
Y_min <- with(df, min(c(male_logit, female_logit)))
Y_max <- with(df, max(c(male_logit, female_logit)))

lapply(seq(0,90,5), function(age) lapply(c(FALSE, TRUE), function(male) 
	plot_CDF(DATASETS, 
	 data = filter_age_sex(df = df, age = age, male=male),
	 age = age,
	 male = male,
	 n_clusters = 2, 
	 h_step = 0.01,
	 Y_min = Y_min,
	 Y_max = Y_max
	)
))


