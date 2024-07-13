library(tidyverse)
library(Cairo)

DATASETS <- "../../../../datasets"


read_year <- function(year, parameter, sample_duration) {
        read.csv(paste0(DATASETS,"/daily_HAPS/daily_HAPS_",year,".csv")) %>%
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
  read.csv(paste0(DATASET, "/IHME_GHDx/Incidence/", filename,"/", filename, ".csv")) %>%
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

comp_p_Female <- function(df) {
  df0 <- df %>% na.omit()
  cor.test(x=df0$geom_mean_nonzero, y=df0$Female, method="kendall", alternative="greater")$p.value
}

comp_p_Male <- function(df) {
  df0 <- df %>% na.omit()
  cor.test(x=df0$geom_mean_nonzero, y=df0$Male, method="kendall", alternative="greater")$p.value
}



read_file <- function(filename) {
  read.csv(paste0(DATASETS, "/IHME_GHDx/Incidence/", filename,"/", filename, ".csv")) %>%
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

fix_age <- function(x) {
  sapply(x, function(y) fix_age_single(y))
}

read_files <- function(filenames) {
map_df(filenames, ~ read_file(.x)) %>%
  filter( !(age %in% c("All ages", "Age-standardized") ) ) %>%
  mutate(age = age %>% fix_age %>% as.integer, sex = as.factor(sex), location = as.factor(location))
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

CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/latency/LatencyArsenicHealthFemale.png"), width=1200, height=1000)
plot(c(), 
     xlim=c(0,10), ylim=c(-150, 0), 
     main="Série chronologique des p-valeurs pour l'hypothèse d'une 
     corrélation de Kendall non-nulle (Femmes)",
     xlab="Latence (Année)",
     ylab="Log-cotes des p-valeurs",
     sub="Source : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency",
     cex.sub=0.7, font.sub=3
)

df_health <- read_files(filenames)

for (age0 in 0:18) {
  x <- 0:10
  y <- sapply(x, function(i) comp_p_Female(make_combination(df_env_yearly, 
          df_health %>% filter(age == age0) %>% select(-age), 
        i)))

  points(x,log10(y/(1-y)))
  lines(x,log10(y/(1-y)))
}

dev.off()


CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/latency/LatencyArsenicHealthMale.png"), width=1200, height=1000)
plot(c(), 
     xlim=c(0,10), ylim=c(-150, 0), 
     main="Série chronologique des p-valeurs pour l'hypothèse d'une 
     corrélation de Kendall non-nulle (Hommes)",
     xlab="Latence (Année)",
     ylab="Log-cotes des p-valeurs",
     sub="Source : Institute for Health Metrics and Evaluation; U.S. Environmental Protection Agency",
     cex.sub=0.7, font.sub=3
)

df_health <- read_files(filenames)

for (age0 in 0:18) {
  x <- 0:10
  y <- sapply(x, function(i) comp_p_Male(make_combination(df_env_yearly, 
          df_health %>% filter(age == age0) %>% select(-age), 
        i)))

  points(x,log10(y/(1-y)))
  lines(x,log10(y/(1-y)))
}

dev.off()


