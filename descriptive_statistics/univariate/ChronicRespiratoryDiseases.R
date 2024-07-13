library(tidyverse)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(Cairo)

DATASETS <- "../../../../datasets"


read_file_standarized <- function(filename) {
read.csv(paste0(DATASETS, "/IHME_GHDx/Incidence/", filename,"/", filename, ".csv")) %>%
	select(location, sex, age, year, val) %>%
	rename(incidence = val) %>% 
	filter(age =="Age-standardized") %>%
	select(-age) %>%
	mutate(sex = as.factor(sex), location = as.factor(location))
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

exp_label <- function(x) {
parse(text = sprintf("10^%s", as.character(log10(x))))
}

male_female_separation <- function(df) {
df_Male <- df %>% filter(sex=="Male") %>% select(-sex) %>% rename(Male = incidence)
df_Female <- df %>% filter(sex=="Female") %>% select(-sex) %>% rename(Female = incidence)
return(merge(df_Male, df_Female))
}


filenames <- c("IHME-GBD_2019_DATA-5fcc42dd-1", "IHME-GBD_2019_DATA-28e7f8dc-1", "IHME-GBD_2019_DATA-77ebed24-1", "IHME-GBD_2019_DATA-610d44b5-1", "IHME-GBD_2019_DATA-a0e15397-1", "IHME-GBD_2019_DATA-c8d94d50-1", "IHME-GBD_2019_DATA-d0cb166f-1")

df <- read_files(filenames)

x_min <- df$incidence %>% min(na.rm=T)
x_max <- df$incidence %>% max(na.rm=T)

df <- df %>% male_female_separation()

df <- df %>%
  mutate(age_group = cut(year, breaks = sapply(seq(from = floor(min(year)), 
                                            to = ceiling(max(year)) + 1,
                                            by = 6), function(x) min(x,2019)),
                          include.lowest = TRUE, 
                          right = FALSE))

CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/boxplotHealthFemale.png"), width=1200, height=1000)
ggplot(df, aes(x=as.factor(age_group), y=Female)) +
  geom_boxplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
  labs(title = "Série chronologique des diagrammes de moustache",
       subtitle = "Chaque boîte contient les données de toutes les tranches d'âge chez les femmes",
       x = "Groupe d'âge",
       y = "Taux d'incidence (malades pour 100,000)",
       caption = "Source : Institute for Health Metrics and Evaluation") +
  theme_classic() +
  facet_wrap(~location, ncol=6) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1),
        axis.ticks.x=element_blank())
dev.off()


CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/boxplotHealthMale.png"), width=1200, height=1000)
ggplot(df, aes(x=as.factor(age_group), y=Male)) +
  geom_boxplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
  labs(title = "Série chronologique des diagrammes de moustache",
       subtitle = "Chaque boîte contient les données de toutes les tranches d'âge chez les hommes",
       x = "Groupe d'âge",
       y = "Taux d'incidence (malades pour 100,000)",
       caption = "Source : Institute for Health Metrics and Evaluation") +
  theme_classic() +
  facet_wrap(~location, ncol=6) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle = 90, hjust = 1),
        axis.ticks.x=element_blank())
dev.off()


variable_labeller <- function(a,i) {
5*a[i]
}


CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/SexComparisonHealth.png"), width=1200, height=1000)
ggplot(df) +
  geom_point(aes(Male, Female), alpha=0.1) +
  geom_smooth(aes(Male, Female), method=MASS::rlm, formula=y ~ x, color="lightgray")+
  coord_fixed(ratio=1, xlim=c(x_min, x_max), ylim=c(x_min, x_max)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = exp_label) +
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
  labs(title="Comparaison du taux d'incidence entre les hommes et les femmes",
     subtitle="Stratification par tranche d'âge",
     caption = "Source : Institute for Health Metrics and Evaluation") +
  xlab("Taux d'incidence chez les hommes (malades pour 100,000)") +
  ylab("Taux d'incidence chez les femmes (malades pour 100,000)") +
  theme_classic() +
  facet_wrap(~age, ncol=5, labeller= variable_labeller) +
  theme(strip.background = element_blank(), strip.placement = "outside")
dev.off()



#########################
# Geospatial Statistics #
#########################

#Load US states map
states <- ne_states(country = "United States of America", returnclass = "sf")

# Choose a specific year
year_of_interest <- 2019

# Filter the dataset for the selected year
df_filtered <- df %>%
  filter(year == year_of_interest, age == 18) %>%
  select(-age)

# Merge the geometric means data with the states map
geo_data <- merge(states, df_filtered, by.x = "name", by.y = "location")

# Plotting
CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/ChronicRespiratoryDiseasesMap.png"), width=1200, height=1000)
ggplot(data = geo_data) +
  geom_sf(aes(fill = Female), color = "black") +
  scale_fill_viridis_c(option = "plasma", name = "Taux d'incidence\nmalades pour 100,000   ") +
  labs(title = paste("Taux d'incidence des maladies respiratoires chroniques chez les femmes (90-94 ans) par État pendant l'année", year_of_interest),
       subtitle = "Source : Institute for Health Metrics and Evaluation") +
  coord_sf(xlim = c(-125, -66.5), ylim = c(24, 49.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(3, "line"),
        legend.title = element_text(margin = margin(b = 20)),
   legend.justification.bottom = "right"
        
 )
dev.off()



