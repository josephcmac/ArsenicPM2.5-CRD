library(tidyverse)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(Cairo)

DATASETS <- "../../../../datasets"


read_year <- function(year, parameter, sample_duration) {
        read.csv(paste0(DATASETS, "/daily_HAPS/daily_HAPS_",year,".csv")) %>%
                filter(Parameter.Name == parameter, Sample.Duration == sample_duration) %>%
                filter(!(State.Name %in% c("Country Of Mexico", "Puerto Rico", "Virgin Islands"))) %>%
                mutate(State.Name = ifelse(State.Name == "District Of Columbia",
                        "District of Columbia", State.Name)) %>%
                select(Date.Local, State.Name, X1st.Max.Value) %>%
                mutate(Date.Local = as.Date(Date.Local)) %>%
                group_by(Date.Local, State.Name) %>%
                summarise(value = ifelse(length(X1st.Max.Value) > 0, max(X1st.Max.Value), NA), .groups = "drop") %>%
                mutate(value = sapply(value, function(x) max(x,0))) %>%
                rename(date = Date.Local, location = State.Name)
}

geom_mean_nonzero <- function(x) {
  x_nonzero <- x[x > 0]
  if (length(x_nonzero) == 0) {
    return(NA)
  }
  x_nonzero %>% log %>% mean %>% exp
}

exp_label <- function(x) {
  parse(text = sprintf("10^%s", as.character(log10(x))))
}

read_file_health <- function(filename, age_group) {
  read.csv(paste0(DATASETS, "/IHME_GHDx/Incidence/", filename,"/", filename, ".csv")) %>%
	  select(location, sex, age, year, val) %>%
	  rename(incidence = val) %>% 
	  filter(age == age_group) %>%
	  select(-age) %>%
	  mutate(sex = as.factor(sex), location = as.factor(location))
}

exp_label <- function(x) {
  parse(text = sprintf("10^%s", as.character(log10(x))))
}

male_female_separation <- function(df) {
  df_Male <- df %>% filter(sex=="Male") %>% select(-sex) %>% rename(Male = incidence)
  df_Female <- df %>% filter(sex=="Female") %>% select(-sex) %>% rename(Female = incidence)
  return(merge(df_Male, df_Female))
}

make_combination <- function(df_env_yearly, df_health, latency) {
  df_env_yearly$year <- df_env_yearly$year + latency 
  df_health <- male_female_separation(df_health) 
  return(merge(df_env_yearly, df_health))
}

comp_p <- function(df) {
  df0 <- df %>% na.omit()
  cor.test(x=df0$geom_mean_nonzero, y=df0$Female, method="kendall", alternative="greater")$p.value
}


read_file <- function(filename) {
  read.csv(paste0(DATASETS,"/IHME_GHDx/Incidence/", filename,"/", filename, ".csv")) %>%
	  select(location, sex, age, year, val) %>%
	  rename(incidence = val)
}

# age group n ranges from 5n to 5n+4
fix_age_single <- function(x) {
  switch(x, 
    "<5 years"="0",
    "5-9 years"="1",
    "10-14 years"="2",
    "15-19 years"="3",
    "20-24 years"="4",
    "25-29 years"="5",
    "30-34 years"="6",
    "35-39 years"="7",
    "40-44 years"="8",
    "45-49 years"="9",
    "50-54 years"="10",
    "55-59 years"="11",
    "60-64 years"="12",
    "65-69 years"="13",
    "70-74 years"="14",
    "75-79 years"="15",
    "80-84"="16",
    "85-89"="17",
    "90-94"="18",
    x
  )
}


fix_age_single0 <- function(x) {
  switch(x, 
    "<20 years"="<20 ans",
    "20-54 years"="20-54 ans",
    "55+ years"="55+ ans",
   x
  )
}

fix_age <- function(x, fix_age_single) {
  sapply(x, function(y) fix_age_single(y))
}

read_files <- function(filenames) {
map_df(filenames, ~ read_file(.x)) %>%
  filter( !(age %in% c("All ages", "Age-standardized") ) ) %>%
  mutate(age = fix_age(age, fix_age_single) %>% as.integer, sex = as.factor(sex), location = as.factor(location))
}

read_files0 <- function(filenames) {
map_df(filenames, ~ read_file(.x)) %>%
  filter( !(age %in% c("All ages", "Age-standardized") ) ) %>%
  mutate(age = as.factor(fix_age(age,fix_age_single0)), sex = as.factor(sex), location = as.factor(location))
}

df_env <- map_df(1988:2019, ~ read_year(.x, parameter="Arsenic PM2.5 LC", sample_duration="24 HOUR")) %>% 
        mutate(location = as.factor(location))

df_env_yearly <- df_env %>%
    mutate(year = floor_date(date, "year") %>% year) %>% 
    group_by(year, location) %>%
    summarise(geom_mean_nonzero=geom_mean_nonzero(value), .groups = "drop") 
rm(df_env)
gc()


filenames <- c("IHME-GBD_2019_DATA-5fcc42dd-1", "IHME-GBD_2019_DATA-28e7f8dc-1", "IHME-GBD_2019_DATA-77ebed24-1", "IHME-GBD_2019_DATA-610d44b5-1", "IHME-GBD_2019_DATA-a0e15397-1", "IHME-GBD_2019_DATA-c8d94d50-1", "IHME-GBD_2019_DATA-d0cb166f-1")

latency <- 6

df_health <- read_files(filenames)

df2 <- map_df(0:18, function(age0) {
  df <- make_combination(df_env_yearly, 
          df_health %>% filter(age == age0) %>% select(-age), latency) %>%
    na.omit()
  df1 <- data.frame(geom_mean = df$geom_mean_nonzero, Male = df$Male, Female = df$Female)
  df1$age <- age0
  return(df1)
})


variable_labeller <- function(a,i) {
  5*a[i]
}


CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/bivariate/BioGradiantArsenicHealthFemale.png"), width=1200, height=1000)
ggplot(df2, aes(x=geom_mean, y=Female)) + 
  geom_point(color="gray", alpha=0.1) +
  geom_smooth(color="black", method="gam", formula=y ~ s(x, bs = "cs")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = exp_label) +
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
  theme_classic() +
  facet_wrap(~age, ncol=5, labeller= variable_labeller) +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(title="Gradient Biologique chez les Femmes", subtitle="Latence : 6 ans", caption="Source : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency.") +
  xlab("Moyenne Géométrique Positive (µg/m³ (25°C))") +
  ylab("Incidence chez les Femmes (malades pour 100,000)")
dev.off()


CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/bivariate/BioGradiantArsenicHealthMale.png"), width=1200, height=1000)
ggplot(df2, aes(x=geom_mean, y=Male)) + 
  geom_point(color="gray", alpha=0.1) +
  geom_smooth(color="black", method="gam", formula=y ~ s(x, bs = "cs")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = exp_label) +
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
  theme_classic() +
  facet_wrap(~age, ncol=5, labeller= variable_labeller) +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(title="Gradient Biologique chez les Hommes", subtitle="Latence : 6 ans", caption="Source : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency.") +
  xlab("Moyenne Géométrique Positive (µg/m³ (25°C))") +
  ylab("Incidence chez les Femmes (malades pour 100,000)")
dev.off()


rm(df_health, df2)

df_health <- read_files0("IHME-GBD_2019_DATA-f6a7f43a-1")

df_health %>% glimpse 

df_health$age %>% unique

df2 <- map_df(c("<20 ans", "20-54 ans", "55+ ans"), function(age0) {
  df <- make_combination(df_env_yearly, 
          df_health %>% filter(age == age0) %>% select(-age), latency) %>%
    na.omit()
  df1 <- data.frame(geom_mean = df$geom_mean_nonzero, Male = df$Male, Female = df$Female, location = df$location)
  df1$age <- age0
  return(df1)
})

CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/bivariate/BioGradiantArsenicHealthAgeGroupFemale.png"), width=1200, height=1000)
ggplot(df2, aes(x=geom_mean, y=Female, color=age)) + 
  geom_point(alpha=0.1) +
  geom_smooth(method="gam", formula=y ~ s(x, bs = "cs")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = exp_label) +
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(title="Gradient Biologique chez les Femmes", subtitle="Latence : 6 ans", caption="Source : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency") +
  xlab("Moyenne Géométrique Positive (µg/m³ (25°C))") +
  ylab("Incidence chez les Femmes (malades pour 100,000)")
dev.off()


CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/bivariate/BioGradiantArsenicHealthAgeGroupMale.png"), width=1200, height=1000)
ggplot(df2, aes(x=geom_mean, y=Male, color=age)) + 
  geom_point(alpha=0.1) +
  geom_smooth(method="gam", formula=y ~ s(x, bs = "cs")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = exp_label) +
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(title="Gradient Biologique chez les Hommes", subtitle="Latence : 6 ans", caption="Source : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency") +
  xlab("Moyenne Géométrique Positive (µg/m³ (25°C))") +
  ylab("Incidence chez les Hommes (malades pour 100,000)")
dev.off()


glimpse(df2)



#########################
# Geospatial Statistics #
#########################

#Load US states map
states <- ne_states(country = "United States of America", returnclass = "sf")

# Choose a specific year
year_of_interest <- 2019


###############

model <- MASS::rlm(Female ~ geom_mean + location - 1, 
          data = df2 %>% filter(age == "20-54 ans") %>% select(-age))

df <- as.data.frame(summary(model)$coefficients)

df <- df[2:nrow(df),] %>%
  select("Value")

rownames(df) <- sapply(rownames(df), function(s) gsub("location", "", s))

df <- cbind(location = rownames(df), df)
rownames(df) <- 1:nrow(df)

# Merge the geometric means data with the states map
geo_data <- merge(states, df, by.x = "name", by.y = "location")

# Plotting
CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/bivariate/BivariateMap20-54.png"), width=1200, height=1000)
ggplot(data = geo_data) +
  geom_sf(aes(fill = Value), color = "black") +
  scale_fill_viridis_c(option = "plasma", name = "Contribution\nMalades pour 100,000   ") +
  labs(title = "Contribution régionale pour le taux d'incidence des maladies respiratoires chroniques chez les femmes (20-54 ans) en fonction de la concentration d'arsenic PM2,5",
       subtitle = "Sources : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency") +
  coord_sf(xlim = c(-125, -66.5), ylim = c(24, 49.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(3, "line"),
        legend.title = element_text(margin = margin(b = 20)),
   legend.justification.bottom = "right" 
 )
dev.off()

#############

model <- MASS::rlm(Female ~ geom_mean + location - 1, 
          data = df2 %>% filter(age == "<20 ans") %>% select(-age))

df <- as.data.frame(summary(model)$coefficients)

df <- df[2:nrow(df),] %>%
  select("Value")

rownames(df) <- sapply(rownames(df), function(s) gsub("location", "", s))

df <- cbind(location = rownames(df), df)
rownames(df) <- 1:nrow(df)

# Merge the geometric means data with the states map
geo_data <- merge(states, df, by.x = "name", by.y = "location")

# Plotting
CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/bivariate/BivariateMap20Minus.png"), width=1200, height=1000)
ggplot(data = geo_data) +
  geom_sf(aes(fill = Value), color = "black") +
  scale_fill_viridis_c(option = "plasma", name = "Contribution\nMalades pour 100,000   ") +
  labs(title = "Contribution régionale pour le taux d'incidence des maladies respiratoires chroniques chez les femmes (<20 ans) en fonction de la concentration d'arsenic PM2,5",
       subtitle = "Sources : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency") +
  coord_sf(xlim = c(-125, -66.5), ylim = c(24, 49.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(3, "line"),
        legend.title = element_text(margin = margin(b = 20)),
   legend.justification.bottom = "right" 
 )
dev.off()

#############

model <- MASS::rlm(Female ~ geom_mean + location - 1, 
          data = df2 %>% filter(age == "55+ ans") %>% select(-age))

df <- as.data.frame(summary(model)$coefficients)

df <- df[2:nrow(df),] %>%
  select("Value")

rownames(df) <- sapply(rownames(df), function(s) gsub("location", "", s))

df <- cbind(location = rownames(df), df)
rownames(df) <- 1:nrow(df)

# Merge the geometric means data with the states map
geo_data <- merge(states, df, by.x = "name", by.y = "location")

# Plotting
CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/bivariate/BivariateMap55Plus.png"), width=1200, height=1000)
ggplot(data = geo_data) +
  geom_sf(aes(fill = Value), color = "black") +
  scale_fill_viridis_c(option = "plasma", name = "Contribution\nMalades pour 100,000   ") +
  labs(title = "Contribution régionale pour le taux d'incidence des maladies respiratoires chroniques chez les femmes (55+ ans) en fonction de la concentration d'arsenic PM2,5",
       subtitle = "Sources : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency") +
  coord_sf(xlim = c(-125, -66.5), ylim = c(24, 49.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(3, "line"),
        legend.title = element_text(margin = margin(b = 20)),
   legend.justification.bottom = "right" 
 )
dev.off()

