library(tidyverse)
source("estimation.R")

compute_params_model_list <- function(params, n_clusters) {
  cov_format <- function(params, n_clusters) {
    p_list <- with(params, 
		   setNames(as.list(p[1:n_clusters]), paste0("p", 1:n_clusters))
    )
    
    mu_list <- with(params, as.list(as.vector(mu)))
    names(mu_list) <- sapply(1:n_clusters, 
			     function(i) c(paste0("muX", i), paste0("muY", i))) %>%
	    as.list %>% unlist


    sigma_list <- list()
    for (i in 1:n_clusters) {
      sigma_list[[paste0("sigmaX", i)]] <- sqrt(params$sigma[1, 1, i])
      sigma_list[[paste0("sigmaY", i)]] <- sqrt(params$sigma[2, 2, i])
      sigma_list[[paste0("cov", i)]] <- params$sigma[1, 2, i]
    }
    
    c(p_list, mu_list, sigma_list, list(loglik = params$loglik))
  }

  pearson_format <- function(params, n_clusters) {
    rho_list <- list()
    for (i in 1:n_clusters) {
      rho_list[[paste0("rho", i)]] <- params[[paste0("cov", i)]] / (params[[paste0("sigmaX", i)]] * params[[paste0("sigmaY", i)]])
    }
    rho_list
  }

  delete_cov <- function(params, n_clusters) {
    params[!names(params) %in% paste0("cov", 1:n_clusters)]
  }

  pretty_format <- function(cov_params, n_clusters) {
    c(delete_cov(cov_params, n_clusters), pearson_format(cov_params, n_clusters))
  }

  cov_params <- cov_format(params, n_clusters)
  pretty_format(cov_params, n_clusters)
}

########### Table Parameters

filter_age_sex <- function(df, age_min0, sex_male) {
  filter_age_sex_male <- function(data, age_min0) {
    data %>%
      select(-female_logit) %>%
      rename(X = ArsenicPM2.5LC_log, Y = male_logit)
  }

  filter_age_sex_female <- function(data, age_min0) {
    data %>%
      select(-male_logit) %>%
      rename(X = ArsenicPM2.5LC_log, Y = female_logit)
  }

  if (sex_male) {
    filter_age_sex_male(
      data = df %>%
        filter(age_min == age_min0) %>%
        select(-age_min, -age_max),
      age_min0
    )
  } else {
    filter_age_sex_female(
      data = df %>%
        filter(age_min == age_min0) %>%
        select(-age_min, -age_max),
      age_min0
    )
  }
}

create_table_age_sex <- function(df, age_min, sex_male, compute_params, n_clusters) {
  c(
    list(
      age_min = age_min, age_max = age_min + 4,
      sex = ifelse(sex_male, "male", "female")
    ),
    compute_params(data = filter_age_sex(
      df = df, age_min = age_min,
      sex_male = sex_male
    ))
  ) %>%
    t() %>%
    as.data.frame()
}

create_table_age <- function(df, age_min, compute_params, n_clusters) {
  rbind(
    create_table_age_sex(df, age_min, sex_male = TRUE, compute_params, n_clusters),
    create_table_age_sex(df, age_min, sex_male = FALSE, compute_params, n_clusters)
  )
}

create_table <- function(df, compute_params, n_clusters) {
  map_df(seq(0, 90, 5), function(age_min) {
    create_table_age(df, age_min = age_min, compute_params, n_clusters)
  }) %>% mutate_at(vars(-sex), as.numeric) %>% mutate(sex = as.character(sex)) 
}



###################### Start


# Example usage
n_clusters <- 4


DATASETS <- "../../../datasets/ArsenicPM2.5-CRD/tables"

create_table(
	df = read.csv(paste0(DATASETS, "/combined.csv")), 
	compute_params = function(data) compute_params_model_list(compute_params_model(data, n_clusters), n_clusters),
	n_clusters = n_clusters
) %>%
	write.csv(paste0(DATASETS, paste0("/parameters",n_clusters,".csv")), row.names=FALSE)



