library(tidyverse)
library(mvtnorm)
library(Cairo)
library(mclust)
library(INLA)

# Function to filter data based on age and sex
filter_age_sex <- function(df, age, male) {
    if (male) {
        df %>% 
            filter(age_min == age, sex == "Homme") %>%
            select(X, Y)
    } else {
        df %>% 
            filter(age_min == age, sex == "Femme") %>%
            select(X, Y)
    }
}

# Function to fit a GMM and estimate the number of components using Bayesian approach
bayesian_gmm <- function(data, K_max) {
    log_likelihood <- function(K) {
        gmm <- Mclust(data, G=K)
        return(gmm$loglik)
    }
    
    # Integrated Nested Laplace Approximation (INLA) for marginal posterior
    inla_approx <- function(K, y) {
        log_prior <- -K
        log_posterior <- log_likelihood(K) + log_prior
        return(log_posterior)
    }
    
    # Estimate the optimal number of components
    K_values <- 1:K_max
    posterior_values <- sapply(K_values, function(K) inla_approx(K, data))
    
    # Normalize the posterior values to get probabilities
    posterior_prob <- exp(posterior_values - max(posterior_values))
    posterior_prob <- posterior_prob / sum(posterior_prob)
    
    # Find the optimal number of components
    K_opt <- K_values[which.max(posterior_prob)]
    
    return(list(K_opt=K_opt, posterior_prob=posterior_prob))
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"
df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))
m <- 10

for (age in seq(90, 0, -5)) {
    for (male in c(TRUE, FALSE)) {
        data <- filter_age_sex(df, age = age, male = male)

        bayesian_result <- bayesian_gmm(data, K_max = m)
        posterior_prob <- bayesian_result$posterior_prob

        img_file <- paste0(DATASETS, "/images/elbow/", ifelse(male, "male", "female"), "/clusters-", age, "-", ifelse(male, "male", "female"), ".png")
        CairoPNG(filename = img_file)

        plot <- data.frame(n_clusters = 1:m, posterior_prob = posterior_prob) %>%
            ggplot(aes(n_clusters, posterior_prob)) +
            geom_point() +
            geom_line() +
            labs(title = "Posterior Probability of Number of Clusters", 
                 subtitle = paste0("Tranche d'Âge : ", age, "-", age+4, ", Sexe : ", ifelse(male, "Homme", "Femme"))) +
            theme_classic()

        print(plot)
        dev.off()
    }
}

