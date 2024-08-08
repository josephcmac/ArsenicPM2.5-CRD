library(tidyverse)
library(mvtnorm)
library(caret)
library(mclust)
library(xtable)

set.seed(700)
source("estimation.R")

new_row <- function(data, n_clusters, age, male, prob=0.8) {
  trainIndex <- createDataPartition(data$Y, p = prob, list = FALSE)
  trainData <- data[trainIndex, ]
  testData <- data[-trainIndex, ]
  
  model <- Mclust(trainData, G = n_clusters)
  
  predicted_clusters <- predict(model, newdata = testData)$classification
  predicted_means <- sapply(predicted_clusters, function(cluster) model$parameters$mean[cluster])
  
  rmse_value <- RMSE(predicted_means, testData$Y)
  mae_value <- MAE(predicted_means, testData$Y)
  
  train_size <- nrow(trainData)
  test_size <- nrow(testData)
  
  data.frame(
    age_min = age,
    age_max = age + 4,
    sex = ifelse(male, "Homme", "Femme"),
    train_size = train_size,
    test_size = test_size,
    RMSE = rmse_value,
    MAE = mae_value
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
                              pull(clusters),
                            age = age, 
                            male = male,
                            prob=0.8)
      )
    }
  }
  d
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

d <- create_table(df = read.csv(paste0(DATASETS, "/tables/combined.csv")),
                  dr = readRDS(paste0(DATASETS, "/tables/clusters.rds"))) %>% 
  mutate(age_min = as.integer(age_min), age_max = as.integer(age_max))

d %>% 
  filter(sex == "Femme") %>% 
  select(-sex) %>% 
  xtable(digits = 2) %>% 
  print(include.rownames = FALSE)


d %>% 
  filter(sex == "Homme") %>% 
  select(-sex) %>% 
  xtable(digits = 2) %>% 
  print(include.rownames = FALSE)
