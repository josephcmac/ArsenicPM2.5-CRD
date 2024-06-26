library(tidyverse)
library(mvtnorm)
library(Cairo)

source("estimation.R")

simulate_model <- function(params, n) {
	sim_format <- function(simulation) {
		data.frame(X = simulation[,1], Y = simulation[,2]) 
	}

  with(params, sapply(1:n, function(rep) {
			      sample(1:length(p), size = 1, prob = p) %>%
				      as.integer %>% 
				      { rmvnorm(1, mean = mu[, .], sigma = sigma[, , .]) }
	
  })) %>% t %>% sim_format
}



filter_age_sex <- function(age, male) {
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


plot_QQ <- function(DATASETS, data, age, male, n_clusters) {
	params <- compute_params_model(data, n_clusters)
	simulation <- simulate_model(params, nrow(data))
	p <- ks.test(data$Y, simulation$Y) %>% getElement("p.value")
	CairoPNG(paste0(DATASETS,"/images/QQ/QQ","-", n_clusters,"-",age,"-",ifelse(male, "male", "female")))
	qqplot(data$Y, simulation$Y, xlab="Données Réelles", ylab="Données Simulées",
		main = "Graphique QQ",
		sub=paste0("Tranche d'Àge : ",age, "-", age+4, "; Sexe : ", ifelse(male, "Homme", "Femme"), "; p (Kolmogorov-Smirnow) = ", round(p, 2)))
	dev.off()	
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

data <- filter_age_sex(age=90, male=TRUE)

n_clusters <- 4



lapply(seq(0,90,5), function(age) lapply(c(FALSE, TRUE), function(male) 
	plot_QQ(
		DATASETS=DATASETS, 
		data=data, 
		age=age, 
		male=male, 
		n_clusters=n_clusters)
))








