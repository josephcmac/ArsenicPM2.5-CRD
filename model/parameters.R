library(tidyverse)
library(mclust)
library(flexclust)
library(mvtnorm)

set.seed(123)

################ Model

simulate_model <- function(params, n) {
  with(params, sapply(1:n, function(rep) {
    switch(EXPR=sample(1:3, size = 1, prob = p),
      rmvnorm(1, mean = mu[, 1], sigma = sigma[, , 1]),
      rmvnorm(1, mean = mu[, 2], sigma = sigma[, , 2]),
      rmvnorm(1, mean = mu[, 3], sigma = sigma[, , 3])
    )
  }))
}

predict_model <- function(params, x) {
  pretty_params_cov <- function(params) {
    with(
      params,
      list(
        p = p,
        mean_x = mu[1, ],
        sd_x = sqrt(sigma[1, 1, 1:3]),
        mean_y = mu[2, ],
        sd_y = sqrt(sigma[2, 2, 1:3]),
        cov_xy = sigma[1, 2, 1:3]
      )
    )
  }

  pretty_params_without_d <- function(params) {
    with(
      pretty_params_cov(params),
      list(
        p = p,
        mean_x = mean_x,
        sd_x = sd_x,
        mean_y = mean_y,
        sd_y = sd_y,
        rho_xy = cov_xy / (sd_x * sd_y),
        d1 = p[1] * dnorm(x, mean = mean_x[1], sd = sd_x[1]),
        d2 = p[2] * dnorm(x, mean = mean_x[2], sd = sd_x[2]),
        d3 = p[3] * dnorm(x, mean = mean_x[3], sd = sd_x[3])
      )
    )
  }

  pretty_params_with_d <- function(params) {
    with(pretty_params_without_d(params), list(
      p = p,
      mean_x = mean_x,
      sd_x = sd_x,
      mean_y = mean_y,
      sd_y = sd_y,
      rho_xy = rho_xy,
      d1 = d1,
      d2 = d2,
      d3 = d3,
      d = d1 + d2 + d3
    ))
  }


  p_and_m <- function(params) {
    with(pretty_params_with_d(params), {
      list(
        p1_cond = d1 / d,
        p2_cond = d2 / d,
        p3_cond = d3 / d,
        m1_cond = mean_y[1] + rho_xy[1] * sd_y[1] * (x - mean_x[1]) / sd_x[1],
        m2_cond = mean_y[2] + rho_xy[2] * sd_y[2] * (x - mean_x[2]) / sd_x[2],
        m3_cond = mean_y[3] + rho_xy[3] * sd_y[3] * (x - mean_x[3]) / sd_x[3]
      )
    })
  }

  with(p_and_m(params),
    p1_cond*m1_cond + p2_cond*m2_cond + p3_cond*m3_cond)
}

estimate_init_model <- function(data, n_clusters = 3) {
  return(
    with(
      kmeans(
        x = data, 
	centers = kcca(data, k=n_clusters, family=kccaFamily("kmeans"))@centers, # k-means++
	iter.max = 10, nstart = 1,
        algorithm = "Hartigan-Wong", trace = FALSE
      ),
      list(
        means = centers,
        variances = lapply(1:n_clusters, function(i) cov(data[cluster == i, ]))
      )
    )
  )
}

compute_params_model <- function(data, n_clusters = 3) {
  with(
    Mclust(data,
      G = n_clusters,
      modelNames = "VVV",
      initialization = estimate_init_model(data, n_clusters)
    ), with(
      parameters,
      list(
        p = pro,
        mu = mean,
        sigma = variance$sigma,
	loglik = loglik
      )
    )
  )
}



compute_params_model_list <- function(data) {
  cov_format <- function(params) {
    with(
      params,
      list(
        p1 = p[1],
        p2 = p[2],
        p3 = p[3],
        muX1 = mu[1, 1],
        muX2 = mu[1, 2],
        muX3 = mu[1, 3],
        muY1 = mu[2, 1],
        muY2 = mu[2, 2],
        muY3 = mu[2, 3],
        sigmaX1 = sqrt(sigma[1, 1, 1]),
        sigmaX2 = sqrt(sigma[1, 1, 2]),
        sigmaX3 = sqrt(sigma[1, 1, 3]),
        sigmaY1 = sqrt(sigma[2, 2, 1]),
        sigmaY2 = sqrt(sigma[2, 2, 2]),
        sigmaY3 = sqrt(sigma[2, 2, 3]),
        cov1 = sigma[1, 2, 1],
        cov2 = sigma[1, 2, 2],
        cov3 = sigma[1, 2, 3],
	loglik = loglik
      )
    )
  }

  pearson_format <- function(params) {
    with(params, {
      list(
        rho1 = cov1 / (sigmaX1 * sigmaY1),
        rho2 = cov2 / (sigmaX2 * sigmaY2),
        rho3 = cov3 / (sigmaX3 * sigmaY3)
      )
    })
  }

  delete_cov <- function(params) {
    params[c(-16, -17, -18)]
  }

  pretty_format <- function(cov_params) {
    c(delete_cov(cov_params), pearson_format(cov_params))
  }

  data %>%
    compute_params_model() %>%
    cov_format() %>%
    pretty_format()
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


create_table_age_sex <- function(df, age_min, sex_male, compute_params) {
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

create_table_age <- function(df, age_min, compute_params) {
  rbind(
    create_table_age_sex(df, age_min, sex_male = TRUE, compute_params),
    create_table_age_sex(df, age_min, sex_male = FALSE, compute_params)
  )
}

create_table <- function(df, compute_params) {
  map_df(seq(0, 90, 5), function(age_min) {
    create_table_age(df, age_min = age_min, compute_params)
  }) %>% mutate_at(vars(-sex), as.numeric) %>% mutate(sex = as.character(sex)) 
}

############## Start


DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

create_table(
	df = read.csv(paste0(DATASETS, "/combined.csv")), 
	compute_params = compute_params_model_list
) %>%
	write.csv(paste0(DATASETS, "/parameters.csv"), row.names=FALSE)


