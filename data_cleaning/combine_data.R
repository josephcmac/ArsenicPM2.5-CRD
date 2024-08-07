library(tidyverse)
library(Rcpp)

Rcpp::cppFunction("
  double logit100000(double x) {
    return log(x) - log(100000 - x);
  }
")

Rcpp::cppFunction("
  double log_geom_mean_nonzero(NumericVector x) {
    NumericVector x_nonzero = x[x > 0];
    if (x_nonzero.length() == 0) {
      return NA_REAL;
    }
    double sum = 0;
    for (int i = 0; i < x_nonzero.length(); i++) {
      sum += log(x_nonzero[i]);
    }
    return sum / x_nonzero.length();
  }
")

Rcpp::cppFunction('
  int fix_age_single(std::string x) {
    if (x == "<5 years") {
      return 0;
    } else if (x == "5-9 years") {
      return 1;
    } else if (x == "10-14 years") {
      return 2;
    } else if (x == "15-19 years") {
      return 3;
    } else if (x == "20-24 years") {
      return 4;
    } else if (x == "25-29 years") {
      return 5;
    } else if (x == "30-34 years") {
      return 6;
    } else if (x == "35-39 years") {
      return 7;
    } else if (x == "40-44 years") {
      return 8;
    } else if (x == "45-49 years") {
      return 9;
    } else if (x == "50-54 years") {
      return 10;
    } else if (x == "55-59 years") {
      return 11;
    } else if (x == "60-64 years") {
      return 12;
    } else if (x == "65-69 years") {
      return 13;
    } else if (x == "70-74 years") {
      return 14;
    } else if (x == "75-79 years") {
      return 15;
    } else if (x == "80-84") {
      return 16;
    } else if (x == "85-89") {
      return 17;
    } else if (x == "90-94") {
      return 18;
    } else {
      return -1;
    }
  }
')

convert_age_to_group <- function(x) {
  sapply(x, function(y) fix_age_single(y))
}

separate_by_sex <- function(df) {
  male_data <- df %>%
    filter(sex == "Male") %>%
    select(-sex) %>%
    mutate(Male = sapply(incidence, logit100000)) %>%
    select(-incidence)
  female_data <- df %>%
    filter(sex == "Female") %>%
    select(-sex) %>%
    mutate(Female = sapply(incidence, logit100000)) %>%
    select(-incidence)
  return(merge(male_data, female_data))
}

combine_environment_and_health_data <- function(df_env_yearly, df_health, latency) {
  df_env_yearly$year <- df_env_yearly$year + latency
  df_health <- separate_by_sex(df_health)
  return(merge(df_env_yearly, df_health))
}

read_incidence_file <- function(filename, DIR) {
  read.csv(paste0(DIR, "/IHME_GHDx/Incidence/", filename, "/", filename, ".csv")) %>%
    select(location, sex, age, year, val) %>%
    rename(incidence = val)
}


read_incidence_files <- function(filenames, DIR) {
  map_df(filenames, ~ read_incidence_file(.x, DIR)) %>%
    filter(!(age %in% c("All ages", "Age-standardized"))) %>%
    mutate(age = convert_age_to_group(age) %>% as.integer(), sex = as.factor(sex), location = as.factor(location))
}


read_env_year <- function(year, parameter, sample_duration, DIR) {
  read.csv(paste0(DIR,"/daily_HAPS/daily_HAPS_", year, ".csv")) %>%
    filter(Parameter.Name == parameter, Sample.Duration == sample_duration) %>%
    filter(!(State.Name %in% c("Country Of Mexico", "Puerto Rico", "Virgin Islands"))) %>%
    mutate(State.Name = ifelse(State.Name == "District Of Columbia", "District of Columbia", State.Name)) %>%
    select(Date.Local, State.Name, X1st.Max.Value) %>%
    mutate(Date.Local = as.Date(Date.Local)) %>%
    group_by(Date.Local, State.Name) %>%
    summarise(value = ifelse(length(X1st.Max.Value) > 0, max(X1st.Max.Value), NA), .groups = "drop") %>%
    mutate(value = sapply(value, function(x) max(x, 0))) %>%
    rename(date = Date.Local, location = State.Name)
}


read_data <- function(filenames, latency, year0, year1, DIR_env, DIR_health) {
}

DATASETS <- "../../../datasets"

filenames <- c("IHME-GBD_2019_DATA-5fcc42dd-1", "IHME-GBD_2019_DATA-28e7f8dc-1", "IHME-GBD_2019_DATA-77ebed24-1", "IHME-GBD_2019_DATA-610d44b5-1", "IHME-GBD_2019_DATA-a0e15397-1", "IHME-GBD_2019_DATA-c8d94d50-1", "IHME-GBD_2019_DATA-d0cb166f-1")

year0 <- 1988
year1 <- 2019
 

df_env <- map_df(year0:year1, ~ read_env_year(.x, parameter = "Arsenic PM2.5 LC", sample_duration = "24 HOUR", DIR = DATASETS)) %>%
    mutate(location = as.factor(location))

df_env_yearly <- df_env %>%
    mutate(year = floor_date(date, "year") %>% year()) %>%
    group_by(year, location) %>%
    summarise(ArsenicPM2.5LC = log_geom_mean_nonzero(value), .groups = "drop")
  rm(df_env)
  gc()

df_health <- read_incidence_files(filenames, DIR = DATASETS)


combine_df <- function(age0, male, latency) {

  df0 <- combine_environment_and_health_data(
      df_env_yearly,
      df_health %>% filter(age == age0) %>% select(-age), latency
    ) %>%
      na.omit()

  if (male) {
	  df1 <- data.frame(X = df0$ArsenicPM2.5LC, Y = df0$Male)
  } else {
	  df1 <- data.frame(X = df0$ArsenicPM2.5LC, Y = df0$Female)
  }

    df1$age_min <- 5 * age0
    df1$age_max <- 5 * age0 + 4
    df1$sex <- ifelse(male, "Homme", "Femme")

  gc()

  df1 %>% select(age_min, age_max, sex, X, Y) 
}


df <- combine_df(age0=0, male=FALSE, latency=7)
df <- rbind(df, combine_df(age0=1, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=2, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=3, male=FALSE, latency=1))
df <- rbind(df, combine_df(age0=4, male=FALSE, latency=1))
df <- rbind(df, combine_df(age0=5, male=FALSE, latency=1))
df <- rbind(df, combine_df(age0=6, male=FALSE, latency=5))
df <- rbind(df, combine_df(age0=7, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=8, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=9, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=10, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=11, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=12, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=13, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=14, male=FALSE, latency=7))
df <- rbind(df, combine_df(age0=15, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=16, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=17, male=FALSE, latency=6))
df <- rbind(df, combine_df(age0=18, male=FALSE, latency=6))

df <- rbind(df, combine_df(age0=0, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=1, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=2, male=TRUE, latency=3))
df <- rbind(df, combine_df(age0=3, male=TRUE, latency=4))
df <- rbind(df, combine_df(age0=4, male=TRUE, latency=5))
df <- rbind(df, combine_df(age0=5, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=6, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=7, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=8, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=9, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=10, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=11, male=TRUE, latency=7))
df <- rbind(df, combine_df(age0=12, male=TRUE, latency=7))
df <- rbind(df, combine_df(age0=13, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=14, male=TRUE, latency=6))
df <- rbind(df, combine_df(age0=15, male=TRUE, latency=5))
df <- rbind(df, combine_df(age0=16, male=TRUE, latency=5))
df <- rbind(df, combine_df(age0=17, male=TRUE, latency=5))
df <- rbind(df, combine_df(age0=18, male=TRUE, latency=1))


write.csv(df, paste0(DATASETS, "/ArsenicPM2.5-CRD/tables/combined.csv"), row.names = FALSE)

