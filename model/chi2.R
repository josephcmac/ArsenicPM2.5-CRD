library(tidyverse)
library(mvtnorm)
library(xtable)

set.seed(121)
source("estimation.R")

new_row <- function(data, n_clusters, age, male) {
  N <- round( 2*nrow(data)**(2/5) ) # https://www.itl.nist.gov/div898/handbook/prc/section2/prc211.htm
  X_min <- with(data, min(X))
  X_max <- with(data, max(X) + 0.01)
  Y_min <- with(data, min(Y))
  Y_max <- with(data, max(Y) + 0.01)
  X_bin <- seq(X_min, X_max, length.out = N + 1)
  Y_bin <- seq(Y_min, Y_max, length.out = N + 1)
  
  params <- compute_params_model(data, n_clusters)
  
  ObservedX <- with(data, sapply(1:N, function(i) {
    sum( (X_bin[i] <= X) & (X < X_bin[i+1]) )
  }))
  
  
  ProbX <-  with(data, sapply(1:N, function(i) {
    with(params, 
         sapply(1:n_clusters, function(j)
           p[j]*(X_bin[i+1] - X_bin[i])*dnorm( 0.5*(X_bin[i] + X_bin[i+1]), 
                               mean = params$mu[1,j], sd = sqrt(sigma[1,1,j])) 
         ) %>% sum
    )
  }))
  
  ObservedY <- with(data, sapply(1:N, function(i) {
    sum( (Y_bin[i] <= Y) & (Y < Y_bin[i+1]) )
  }))
  
  
  ProbY <-  with(data, sapply(1:N, function(i) {
    with(params, 
         sapply(1:n_clusters, function(j)
           p[j]*(Y_bin[i+1] - Y_bin[i])*dnorm( 0.5*(Y_bin[i] + Y_bin[i+1]), 
                            mean = params$mu[2,j], sd = sqrt(sigma[2,2,j])) 
         ) %>% sum
    )
  }))
  
  data.frame(
    age_min = age,
    age_max = age+4,
    sex = ifelse(male, "Homme", "Femme"),
    bins = N,
    size = nrow(data),
    chi2X = chisq.test(x=ObservedX, p=ProbX, rescale.p=TRUE, 
                       simulate.p.value = TRUE) %>% getElement("p.value"),
    chi2Y = chisq.test(x=ObservedY, p=ProbY, rescale.p=TRUE, 
                       simulate.p.value = TRUE) %>% getElement("p.value")
  )
}

create_table <- function(df, dr) {
  d <- data.frame()
  for (age in seq(0,90,5)) {
    for (male in c(TRUE, FALSE)) {
      d <- rbind(d, new_row(data = filter_age_sex(df, age = age, male = male),
                            n_clusters = dr %>% 
                              filter(
                                age_min == age, 
                                sex == ifelse(male, "Homme", "Femme")
                              ) %>% 
                              getElement("clusters"),
                            age = age, 
                            male = male)
                 )
    }
  }
  d
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

d <- create_table(df = read.csv(paste0(DATASETS, "/tables/combined.csv")),
             dr = readRDS(paste0(DATASETS, "/tables/clusters.rds"))) %>% 
  mutate(age_min = as.integer(age_min), age_max = as.integer(age_max),
                  bins = as.integer(bins))

d %>% 
  filter(sex == "Femme") %>% 
  select(-sex) %>% 
  xtable(digits = 5) %>% 
  print(include.rownames = FALSE)


d %>% 
  filter(sex == "Homme") %>% 
  select(-sex) %>% 
  xtable(digits = 5) %>% 
  print(include.rownames = FALSE)


