library(tidyverse)

set.seed(123)
source("estimation.R")

custom_cdf <- function(params, x) {
  with(params,
       sapply(1:length(p), function(i) p[i] * pnorm(x, mean = params$mu[2,][i], sd = sqrt(sigma[,,i][2,2])) ) %>%
         sum
  )
}

filter_age_sex <- function(df, male, age) {
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

chi_sex_age <- function(df, male, age, num_bins, n_clusters) {
	bootstrap <- function(data) {
		data[sample(1:nrow(data), replace=TRUE, size = nrow(data)),]

	}


	create_table_sex_age_data <- function(data, male, age, num_bins, n_clusters) {
		normalize_add <- function(A) {
			A/sum(A)
		}	

		compute_p <- function(params, num_bins) {
			compute_p_aux <- function(params, bins) {
				sapply(1:(length(bins)-1), function(i) {
					custom_cdf(params, bins[i+1]) - custom_cdf(params, bins[i])
		  		}) %>% normalize_add
			}

			compute_p_aux(params=params, 
				      bins=seq(min(data$Y), max(data$Y), length.out = num_bins + 1))
		}


	with(
	     chisq.test(
			x = with(data, with(hist(Y, 
					       breaks = with(data, seq(min(Y), max(Y), length.out = num_bins + 1)), 
					       plot = FALSE), counts)), 
			p = compute_p(params = compute_params_model(data, n_clusters),
				      num_bins = num_bins)),
	     list(
		  p = p.value
		)
	)


	} 

	create_table_sex_age_data(data = filter_age_sex(df = df, age = age, male = male) %>% bootstrap,
				  male, age, num_bins, n_clusters)
}


calculate_mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}


DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

N <- 1034
num_bins <- round( 2*N**(2/5) ) # https://www.itl.nist.gov/div898/handbook/prc/section2/prc211.htm
rm(N)


r <- sapply(1:100, function(i)
c(
			  chi_sex_age(df, male=TRUE, age=90, num_bins=num_bins, n_clusters=2) %>% as.numeric,
			  chi_sex_age(df, male=TRUE, age=90, num_bins=num_bins, n_clusters=3) %>% as.numeric,
			  chi_sex_age(df, male=TRUE, age=90, num_bins=num_bins, n_clusters=4) %>% as.numeric,
			  chi_sex_age(df, male=TRUE, age=90, num_bins=num_bins, n_clusters=5) %>% as.numeric 
			  ) %>% rank %>% (function(x) x+1)
)

r[4,] %>% plot 

r[4,] %>% calculate_mode





