library(tidyverse)

# Source the required file
source("estimation.R")

# Function to filter data based on age and sex
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

# Function to calculate conditional expectation
calc_esperance_cond <- function(X, params, n_clusters) {
  pretty_params <- function(params, n_clusters) {
    with(params, {
      moys_x <- mu[1, ]
      sds_x <- sqrt(sigma[1, 1, 1:n_clusters])
      moys_y <- mu[2, ]
      sds_y <- sqrt(sigma[2, 2, 1:n_clusters])
      rho_xy <- sigma[1, 2, 1:n_clusters] / (sds_x * sds_y)
      
      list(
        p = p,
        moys_x = moys_x,
        sds_x = sds_x,
        moys_y = moys_y,
        sds_y = sds_y,
        rho_xy = rho_xy
      )
    })
  }
  
  calc_cond_probs <- function(X, moys_x, sds_x, p) {
    sapply(1:n_clusters, function(i) p[i] * dnorm(X, mean = moys_x[i], sd = sds_x[i])) %>%
      { . / rowSums(.) }
  }
  
  calc_cond_means <- function(X, moys_x, sds_x, moys_y, sds_y, rho_xy) {
    sapply(1:n_clusters, function(i) moys_y[i] + rho_xy[i] * sds_y[i] * (X - moys_x[i]) / sds_x[i])
  }
  
  params_list <- pretty_params(params, n_clusters)
  
  with(params_list, {
    cond_probs <- calc_cond_probs(X, moys_x, sds_x, p)
    cond_means <- calc_cond_means(X, moys_x, sds_x, moys_y, sds_y, rho_xy)
    rowSums(cond_probs * cond_means)
  })
}

# Function to calculate mean absolute error
mae <- function(x, y) {
  mean(abs((x - y)))
}

# Function to calculate mean squared error
mse <- function(x, y) {
  mean((x - y)^2)
}

# Function to compute error based on age, sex, alpha, and number of clusters
error_age_sex <- function(df, age, male, alpha, n_clusters, error_func) {
  
  bootstrap <- function(data) {
    data[sample(1:nrow(data), replace = TRUE), ]
  }
  
  split_data <- function(data, alpha) {
    idx <- sample(1:nrow(data), size = round(alpha * nrow(data)))
    list(train = data[idx, ], test = data[-idx, ])
  }
  
  compute_error <- function(train, test) {
    params <- compute_params_model(train, n_clusters)
    Y_pred <- calc_esperance_cond(X = test$X, params = params, n_clusters = n_clusters)
    error_func(test$Y, Y_pred)
  }
  
  data_filtered <- filter_age_sex(df, age, male) %>% bootstrap()
  data_split <- split_data(data_filtered, alpha)
  compute_error(data_split$train, data_split$test)
}

# Function to calculate the mode
calculate_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

set.seed(123)

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"
df <- read.csv(file.path(DATASETS, "tables/combined.csv"))

alpha <- 0.3

# Add debugging information
print("Starting error calculations")

# Running the error calculation
r <- sapply(1:100, function(i) {
  cat("Iteration:", i, "\n")
  c(
    error_age_sex(df, age = 90, male = TRUE, alpha = alpha, n_clusters = 2, error_func = mae),
    error_age_sex(df, age = 90, male = TRUE, alpha = alpha, n_clusters = 3, error_func = mae),
    error_age_sex(df, age = 90, male = TRUE, alpha = alpha, n_clusters = 4, error_func = mae),
    error_age_sex(df, age = 90, male = TRUE, alpha = alpha, n_clusters = 5, error_func = mae)
  ) %>% rank %>% { . + 1 }
})

print("First error calculation completed")

r[1, ] %>% plot

r[1, ] %>% calculate_mode

#--------------------------

print("Second error calculation starting")

r <- sapply(1:100, function(i) {
  cat("Iteration:", i, "\n")
  c(
    error_age_sex(df, age = 90, male = TRUE, alpha = alpha, n_clusters = 2, error_func = mse),
    error_age_sex(df, age = 90, male = TRUE, alpha = alpha, n_clusters = 3, error_func = mse),
    error_age_sex(df, age = 90, male = TRUE, alpha = alpha, n_clusters = 4, error_func = mse),
    error_age_sex(df, age = 90, male = TRUE, alpha = alpha, n_clusters = 5, error_func = mse)
  ) %>% rank %>% { . + 1 }
})

print("Second error calculation completed")


r[1, ] %>% plot

r[1, ] %>% calculate_mode

