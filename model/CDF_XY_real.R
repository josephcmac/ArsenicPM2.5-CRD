library(tidyverse)
library(mvtnorm)
library(Cairo)
library(RColorBrewer)

set.seed(9809)


source("estimation.R")

produce_image <- function(data, n_clusters, age, male, 
                          X_bin, Y_bin, X_lab, Y_lab) {
  
  Observed <- sapply(1:X_len, function(i) {
    sapply(1:Y_len, function(j) {
      with(data, sum((X_bin[i] <= X) & (X < X_bin[i + 1]) & 
                       (Y_bin[j] <= Y) & (Y < Y_bin[j + 1])))
    })
  })
  
  params <- compute_params_model(data, n_clusters)
  Expected <- nrow(data) * sapply(1:X_len, function(i) {
    sapply(1:Y_len, function(j) {
      with(data, with(
        params,
        sapply(1:n_clusters, function(k) {
          (Y_bin[j + 1] - Y_bin[j]) * (X_bin[i + 1] - 
                                         X_bin[i]) * p[k] * mvtnorm::dmvnorm(
                                        x = c(0.5 * (X_bin[i] + X_bin[i + 1]), 
                                             0.5 * (Y_bin[j] + Y_bin[j + 1])),
                                           mean = mu[, k], sigma = sigma[, , k]
                                         )
        }) %>% sum()
      ))
    })
  })
  
  plot_combined_matrices <- function(observed_matrix, expected_matrix, 
                                     file_name, title1, title2, 
                                     X_bin, Y_bin, X_lab, Y_lab, age, male) {
    
    CairoPNG(file_name, width = 1400, height = 700)
    par(mfrow = c(1, 2))  
    par(mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 2, 0))  
    n_labels_x <- 10  
    n_labels_y <- 10
    
    X_labels_to_show <- seq(1, length(X_bin), length.out = n_labels_x)
    Y_labels_to_show <- seq(1, length(Y_bin), length.out = n_labels_y)
    
    # Plot observed matrix
    image(
      x = X_bin, y = Y_bin, z = t(observed_matrix),
      col = gray.colors(100, start = 0, end = 1),
      axes = FALSE,
      main = title1,
      xlab = "Arsenic PM2.5 (microgrammes/mètre cube (LC))",
      ylab = "Taux d'Incidence (nouvaux cas pour 100 000 personnes)"
    )
    axis(1, at = X_bin[X_labels_to_show], labels = X_lab[X_labels_to_show], cex.axis = 0.7)
    axis(2, at = Y_bin[Y_labels_to_show], labels = Y_lab[Y_labels_to_show], cex.axis = 0.7)
    box()
    
    # Plot expected matrix
    image(
      x = X_bin, y = Y_bin, z = t(expected_matrix),
      col = gray.colors(100, start = 0, end = 1),
      axes = FALSE,
      main = title2,
      xlab = "Arsenic PM2.5 (microgrammes/mètre cube (LC))",
      ylab = "Taux d'Incidence (nouvaux cas pour 100 000 personnes)"
    )
    axis(1, at = X_bin[X_labels_to_show], labels = X_lab[X_labels_to_show], cex.axis = 0.7)
    axis(2, at = Y_bin[Y_labels_to_show], labels = Y_lab[Y_labels_to_show], cex.axis = 0.7)
    box()
    
    mtext(paste0("Tranche d'Âge : ", age, "-", age + 4, ", Sexe : ", 
                 ifelse(male, "Homme", "Femme")), 
          side = 3, outer = TRUE, cex = 1.5)
    dev.off()
  }
  
  
  combined_file <- paste0(DATASETS, "/images/CDF_XY_real/", 
                          ifelse(male, "male", "female"),"/combined_CDF_XY-", 
                          age, "-", ifelse(male, "male", "female"), ".png")
  
  plot_combined_matrices(Observed, Expected, combined_file, 
                         "Observées", "Espérées", X_bin, Y_bin, X_lab, Y_lab,
                         age, male)
}

expint <- function(y) {
  v <- exp(y)
  100000*v/(v+1)
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

X_len <- 50
Y_len <- 50


df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

dr <- readRDS(paste0(DATASETS, "/tables/clusters.rds"))

# Calculate global X and Y limits across the entire dataset
X_global_min <- min(df$X)
X_global_max <- max(df$X)
Y_global_min <- min(df$Y)
Y_global_max <- max(df$Y)

X_bin <- seq(X_global_min, X_global_max, length.out = X_len + 1)
Y_bin <- seq(Y_global_min, Y_global_max, length.out = Y_len + 1)

X_lab <- sapply(X_bin, exp) %>% round(6)
Y_lab <- sapply(Y_bin, expint) %>% round(0)

for (age in seq(0, 90, 5)) {
  for (male in c(TRUE, FALSE))
    produce_image(data=filter_age_sex(df, age = age, male = male), 
                  n_clusters=dr %>% 
                    filter(
                      age_min == age, 
                      sex == ifelse(male, "Homme", "Femme")
                    ) %>% 
                    getElement("clusters"), 
                  age=age, male=male,
                  X_bin=X_bin, Y_bin=Y_bin, 
                  X_lab=X_lab, Y_lab=Y_lab)
}


