library(tidyverse)
library(xtable)
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

mae <- function(x,y) {
	(x-y) %>% abs %>% mean
}


create_table_age_sex <- function(df, age, male, alpha, n_clusters, error) {
	split_data <- function(data, alpha) {
		sample(1:nrow(data), size = round(alpha*nrow(data))) %>%
		list(
		     train = data[.,],
		     test = data[-.,]
		)
	}

	compute_error <- function(train, test) {
		with(test,{ 
		 Y_pred <- calc_esperance_cond(X=X, 
				params=compute_params_model(train, n_clusters),
				n_clusters=n_clusters)
		 list(
		      age_min = age,
		      age_max = age + 4,
		      sex = ifelse(male, "Homme", "Femme"),
		      mae = error(Y, Y_pred)) 
		}
		)
	}

 
	with(split_data(data=filter_age_sex(df=df, age=age, male=male), alpha=alpha),
	     compute_error(train, test)
	)
}


create_table_sex <- function(df, male, alpha, n_clusters, error) {
	sapply(seq(0,90,5), function(age)
	       create_table_age_sex(df=df, age=age, male=male, alpha = alpha, n_clusters = n_clusters, error = error)
	) %>% t %>% as.data.frame
}


create_table <- function(df, alpha, n_clusters, error) {
	rbind(
		create_table_sex(df, male=TRUE, alpha = alpha, n_clusters, error = error),
		create_table_sex(df, male=FALSE, alpha = alpha, n_clusters, error = error)
	) %>%
	mutate(
	       age_min = as.integer(age_min),
	       age_max = as.integer(age_max),
	       sex = as.character(sex),
	       mae = as.numeric(mae)
	)
}

set.seed(123)

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

d <- create_table(df, alpha = 0.8, n_clusters = 2, error = mae) %>% 
	rename(mae2 = mae)

d_next <- create_table(df, alpha = 0.8, n_clusters = 3, error = mae) %>% 
	rename(mae3 = mae)

d <- merge(d, d_next)

d_next <- create_table(df, alpha = 0.8, n_clusters = 4, error = mae) %>% 
	rename(mae4 = mae)

d <- merge(d, d_next)

d_next <- create_table(df, alpha = 0.8, n_clusters = 5, error = mae) %>% 
	rename(mae5 = mae)

d <- merge(d, d_next)
rm(d_next)

index <- sort(d$age_min, index.return=TRUE)

table_latex <- xtable(d[index$ix,])
print(table_latex, type = "latex", include.rownames = FALSE, booktabs = TRUE)


