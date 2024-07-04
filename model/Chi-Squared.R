library(tidyverse)
library(xtable)

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

create_table_sex_age <- function(df, male, age, num_bins, n_clusters) {
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
	     		age_min = age,
			age_max = age+4,
			sex = ifelse(male, "Homme", "Femme"),
			p = p.value
		)
	)


	} 

	create_table_sex_age_data(data = filter_age_sex(df = df, age = age, male = male),
				  male, age, num_bins, n_clusters)
}




create_table_sex <- function(df, male, num_bins, n_clusters) {
	sapply(seq(0, 90, 5), function(age) create_table_sex_age(df, male, age, num_bins, n_clusters)
	) %>% t %>% as.data.frame %>% 
	mutate(
	       age_min = as.integer(age_min),
	       age_max = as.integer(age_max),
	       sex = as.character(sex),
	       p = as.numeric(p)
	)
}

create_table <- function(df, num_bins, n_clusters) {
	rbind(create_table_sex(df, TRUE, num_bins, n_clusters), create_table_sex(df, FALSE, num_bins, n_clusters))
}



DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

N <- 1034
num_bins <- round( 2*N**(2/5) ) # https://www.itl.nist.gov/div898/handbook/prc/section2/prc211.htm
rm(N)

d <- create_table(df, num_bins=num_bins , n_clusters=2) %>% 
	rename(p2 = p)

d_next <- create_table(df, num_bins=num_bins , n_clusters=3) %>% 
	rename(p3 = p)

d <- merge(d, d_next)

d_next <- create_table(df, num_bins=num_bins , n_clusters=4) %>% 
	rename(p4 = p)

d <- merge(d, d_next)

d_next <- create_table(df, num_bins=num_bins , n_clusters=5) %>% 
	rename(p5 = p)

d <- merge(d, d_next)
rm(d_next)

index <- sort(d$age_min, index.return=TRUE)

table_latex <- xtable(d[index$ix,])
print(table_latex, type = "latex", include.rownames = FALSE, booktabs = TRUE)




