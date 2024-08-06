library(tidyverse)
library(Cairo)
source("estimation.R")

compute_params_model_list <- function(params, n_clusters) {
  cov_format <- function(params, n_clusters) {
    p_list <- list()
    p_list[["p"]] <- params$p
    
    mu_list <- list()
    mu_list[["muX"]] <- params$mu[1,]
    mu_list[["muY"]] <- params$mu[2,]
    
    sigma_list <- list()
    for (i in 1:n_clusters) {
      sigma_list[["sigmaX"]][i] <- sqrt(params$sigma[1, 1, i])
      sigma_list[["sigmaY"]][i] <- sqrt(params$sigma[2, 2, i])
      sigma_list[["cov"]][i] <- params$sigma[1, 2, i]
    }
    
    c(p_list, mu_list, sigma_list)
  }
  
  pearson_format <- function(params, n_clusters) {
    rho_list <- list()
    for (i in 1:n_clusters) {
      rho_list[["rho"]][i] <- params[["cov"]][i] / 
        (params[["sigmaX"]][i] * params[["sigmaY"]][i] )
    }
    rho_list
  }
  
  delete_cov <- function(params, n_clusters) {
    params[!names(params) %in% c("cov")]
  }
  
  pretty_format <- function(cov_params, n_clusters) {
    c(delete_cov(cov_params, n_clusters), 
      pearson_format(cov_params, n_clusters))
  }
  
  pretty_format(cov_format(params, n_clusters), n_clusters)
}


produce_graphics <- function(data, n_clusters, age_min, male) {
  produce_graphics_aux <- function(pretty_params, n_clusters, age_min, male) {
    list(
      g1=ggplot(pretty_params) +
        geom_label(aes(muX, sigmaX, label=1:n_clusters)) + 
        xlim(-9,-5) +
        ylim(0,2) +
        labs(x = "Moyenne X", y = "Écart type X",
             title="Distribution de X",
             subtitle=paste0("Tranche d'Âge : ",age_min, "-", age_min+4, " ; ", 
                             "Sexe : ", ifelse(male, "Homme", "Femme"))) +
        theme_classic() +
        theme(axis.title = element_text(size = 14), 
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 16, face = "bold")),
      
      g2=ggplot(pretty_params) +
        geom_label(aes(muY, sigmaY, label=1:n_clusters)) + 
        xlim(-1, -7) +
        ylim(0,1) +
        labs(x = "Moyenne Y", y = "Écart type Y", 
             title="Distribution de Y",
             subtitle=paste0("Tranche d'Âge : ",age_min, "-", age_min+4, " ; ", 
                             "Sexe : ", ifelse(male, "Homme", "Femme"))) +
        theme_classic() +
        theme(axis.title = element_text(size = 14), 
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 16, face = "bold")),
      
      g3=ggplot(pretty_params, aes(x = rho, y = y)) +
        geom_point(size = 10, color="white") +
        annotate("segment", x = -1, xend = 1, y = 0, yend = 0, size = 2) +
        annotate("segment", x = -1, xend = -1, y = -0.1, yend = 0.1, size = 2) +
        annotate("segment", x = 1, xend = 1, y = -0.1, yend = 0.1, size = 2) +
        geom_text(aes(label = 1:n_clusters), col = "black") +
        annotate("text", x = -1, y = -0.3, label = "-1", color = "black", size = 5, 
                 vjust = -1) +
        annotate("text", x = 1, y = -0.3, label = "+1", color = "black", size = 5, 
                 vjust = -1) +
        geom_segment(aes(x = rho, xend = rho, 
                         y = y, yend = 0), 
                     color = "gray", linetype = "dashed") +
        geom_point(aes(x = rho, y = 0), color = "red", size = 5) +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
        scale_y_continuous(limits = c(-1, 1)) +
        labs(title="Corrélation de Pearson entre X et Y",
             subtitle=paste0("Tranche d'Âge : ",age_min, "-", age_min+4, " ; ", 
                             "Sexe : ", ifelse(male, "Homme", "Femme"))) +
        theme(panel.background = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(size = 16, face = "bold")),
      
      g4=ggplot(pretty_params, aes(x = "", y = p, fill = factor(1:nrow(pretty_params)))) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar(theta = "y") +
        geom_text(aes(label = 1:n_clusters), 
                  position = position_stack(vjust = 0.5), 
                  size = 5, color = "white") +
        labs(title = "Poids du Mélange",
             subtitle=paste0("Tranche d'Âge : ",age_min, "-", age_min+4, " ; ", 
                             "Sexe : ", ifelse(male, "Homme", "Femme")),
             fill = "Cluster") +
        theme_void() + 
        theme(legend.position = "",
              plot.title = element_text(size = 16, face = "bold"))
    )
  }
  
  clear_params2 <- function(pretty_params, n_clusters) {
    if(n_clusters %% 2 == 0) {
      pretty_params$y <- seq(-0.5, 0.5, length.out = n_clusters)
    } else {
      pretty_params$y <- setdiff(seq(-0.5, 0.5, length.out = n_clusters+1), -0.5)
    }
    pretty_params
  }
  
  clean_params1 <- function(data, n_clusters) {
    compute_params_model(data, n_clusters) %>% 
      compute_params_model_list(n_clusters) %>%
      as.data.frame
    
  }
  
  clean_params1(data, n_clusters) %>% 
    clear_params2(n_clusters) %>%
    produce_graphics_aux(n_clusters, age_min, male)
}

save_pictures <- function(df, dr, DATASETS) {
  print_images <- function(g, age, male, DATASETS) {
    with(g, {
      CairoPNG(paste0(DATASETS,"/images/params/", 
                      ifelse(male, "male", "female"),"/g1-",age,"-",
                      ifelse(male, "male", "female"),".png"))
      print(g1)
      dev.off()
      
      CairoPNG(paste0(DATASETS,"/images/params/", 
                      ifelse(male, "male", "female"),"/g2-",age,"-",
                      ifelse(male, "male", "female"),".png"))
      print(g2)
      dev.off()
      
      CairoPNG(paste0(DATASETS,"/images/params/", 
                      ifelse(male, "male", "female"),"/g3-",age,"-",
                      ifelse(male, "male", "female"),".png"))
      print(g3)
      dev.off()
      
      CairoPNG(paste0(DATASETS,"/images/params/", 
                      ifelse(male, "male", "female"),"/g4-",age,"-",
                      ifelse(male, "male", "female"),".png"))
      print(g4)
      dev.off()
      
    })
  }
  
  for (age in seq(0, 90, 5)) {
    for (male in c(TRUE, FALSE)) {
      produce_graphics(
        data = filter_age_sex(df, age=age, male=male),
        n_clusters = dr %>% 
          filter(
            age_min == age, 
            sex == ifelse(male, "Homme", "Femme")
          ) %>% 
          getElement("clusters"),
        age_min=age,
        male=male
      ) %>% print_images(age, male, DATASETS)
    }
  }
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))
dr <- readRDS(paste0(DATASETS, "/tables/clusters.rds"))

save_pictures(df, dr, DATASETS)

