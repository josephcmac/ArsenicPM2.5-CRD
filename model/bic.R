library(tidyverse)
library(Cairo)
library(xtable)

set.seed(123)
source("estimation.R")

custom_cdf <- function(params, x) {
  with(params,
       sapply(1:length(p), function(i) p[i] * pnorm(x, mean = params$mu[2,][i], 
                                              sd = sqrt(sigma[,,i][2,2])) ) %>%
         sum
  )
}

create_table_sex_age <- function(df, male, age, n_clusters) {
	create_table_sex_age_data <- function(data, male, age, n_clusters) {
		params <- compute_params_model(data, n_clusters)
		with(params,
			list(
			     age_min = age,
			     age_max = age+4,
			     sex = ifelse(male, "Homme", "Femme"),
			     bic = bic
			)
		)
	}

	create_table_sex_age_data(
	  data = filter_age_sex(df = df, age = age, male = male),
				  male, age, n_clusters)
}




create_table_sex <- function(df, male, n_clusters) {
	sapply(seq(0, 90, 5), function(age) 
	  create_table_sex_age(df, male, age, n_clusters)) %>% 
		t %>% as.data.frame %>% 
	mutate(
	       age_min = as.integer(age_min),
	       age_max = as.integer(age_max),
	       sex = as.character(sex),
	       bic = as.numeric(bic)
	)
}

create_table <- function(df, n_clusters) {
	rbind(create_table_sex(df, TRUE, n_clusters), 
	      create_table_sex(df, FALSE, n_clusters))
}

DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

K_max <- 10

df <- read.csv(paste0(DATASETS, "/tables/combined.csv"))

# Initialize the first table with n_clusters=2
d <- create_table(df, n_clusters=2) %>%
  rename(bic2 = bic)

# Loop over the remaining n_clusters values
for (k in 3:K_max) {
  d_next <- create_table(df, n_clusters=k) %>%
    rename(!!paste0("bic", k) := bic)
  d <- merge(d, d_next)
}

rm(d_next)

index <- sort(d$age_min, index.return=TRUE)

d <- d[index$ix,]
rm(index)


d_long <- d %>%
	filter(sex == "Homme") %>%
	select(-sex, -age_max) %>%
	reshape2::melt(id.vars = c("age_min"), variable.name = "k", value.name = "p_value")

CairoPNG(paste0(DATASETS,"/images/BIC/BIC-male.png"))
ggplot(d_long, aes(x = age_min, y = p_value, color = k)) +
  geom_point() +
  geom_line() +
  labs(title = "BIC par tranche d'âge (hommes)",
       x = "Tranche d'Âge",
       y = "BIC",
       color = "k") +
  theme_minimal()
dev.off()
# ---------------------

d_long <- d %>%
	filter(sex == "Femme") %>%
	select(-sex, -age_max) %>%
	reshape2::melt(id.vars = c("age_min"), variable.name = "k", value.name = "p_value")

CairoPNG(paste0(DATASETS,"/images/BIC/BIC-female.png"))
ggplot(d_long, aes(x = age_min, y = p_value, color = k)) +
  geom_point() +
  geom_line() +
  labs(title = "BIC par tranche d'âge (femmes)",
       x = "Tranche d'Âge",
       y = "BIC",
       color = "k") +
  theme_minimal()
dev.off()




dr <- d %>%
  rowwise() %>%
  mutate(clusters = 1+which.max(c_across(bic2:bic10))) %>%
  ungroup() %>%
  mutate(clusters = as.integer(clusters)) %>%
  select(age_min, age_max, sex, clusters)

dr %>% filter(sex=="Homme") %>% select(-sex) %>% xtable() %>% 
  print(include.rownames=FALSE)

dr %>% filter(sex=="Femme") %>% select(-sex) %>% xtable() %>% 
  print(include.rownames=FALSE)

saveRDS(dr, file=paste0(DATASETS, "/tables/clusters.rds"))

