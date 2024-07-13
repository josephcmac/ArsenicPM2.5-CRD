library(tidyverse)
library(scales)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(Cairo)

DATASETS <- "../../../../datasets"

# Preliminaries
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


########################
# Microscopic Analysis #
########################

# Reading
df <- map_df(1988:2023, ~ read_year(.x, parameter="Arsenic PM2.5 LC", sample_duration="24 HOUR")) %>% 
	mutate(location = as.factor(location))

# Overview
glimpse(df)

nrow(df)

length(unique(df$location))

length(unique(df$date))

summary(df)

(mean( df$value > 0 )*100) %>% round(0)


# Longitudinal study 
CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/SmoothArsenicDaily.png"), width=1200, height=1000)
df %>% 
	ggplot(aes(x=date, y=value))   +
	geom_smooth(color="black", method = 'gam', formula = y ~ s(x, bs = "cs")) +
	labs(title = "Concentration (lissée) d'arsenic PM2,5 LC",
	     x = "Année",
	     y = "Concentration (µg/m³ (25°C))",
	     caption = "Source : U.S. Environmental Protection Agency") +
	theme_classic() +
	facet_wrap(~location, ncol=6)
dev.off()

CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/SmoothArsenicDailyPositive.png"), width=1200, height=1000)
df %>% 
	filter(value > 0) %>%
	ggplot(aes(x=date, y=value))   +
	geom_smooth(color="black", method = 'gam', formula = y ~ s(x, bs = "cs")) +
	labs(title = "Concentration (lissée) des valeurs positives d'arsenic PM2,5 LC",
	     x = "Année",
	     y = "Concentration Positives (µg/m³ (25°C))",
	     caption = "Source : U.S. Environmental Protection Agency") +
	scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
	theme_classic() +
	facet_wrap(~location, ncol=6)
dev.off()

df <- df %>%
  mutate(year = floor_date(date, "year") %>% year,  
	 age_group = cut(year, 
                         breaks = sapply(seq(from = floor(min(year)), 
                                            to = ceiling(max(year)) + 1, 
                                            by = 6), 
                                         function(x) min(x, 2023)),
                         include.lowest = TRUE, 
                         right = FALSE))

CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/BoxplotArsenicDailyPositve.png"), width=1200, height=1000)
ggplot(df %>% filter(value > 0), aes(x=as.factor(age_group), y=value)) +
  geom_boxplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
  labs(title = "Série chronologique des diagrammes de moustache",
    x = "Année",
    y = "Moyenne géométrique positive (µg/m³ (25°C))",
    caption = "Source : U.S. Environmental Protection Agency") +
  theme_classic() +
  facet_wrap(~location, ncol=6) +
  theme(axis.title.x=element_blank(),
  axis.text.x=element_text(angle = 90, hjust = 1),
  axis.ticks.x=element_blank())
dev.off()

########################
# Macroscopic Analysis #
########################
df_yearly <- df %>%
  mutate(year = floor_date(date, "year") %>% year) %>% 
  group_by(year, location) %>%
  summarise(geom_mean_nonzero=geom_mean_nonzero(value), .groups = "drop") 
rm(df)
gc(df)


CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/ArsenicGeomMean.png"), width=1200, height=1000)
ggplot(df_yearly %>% na.omit) + 
	geom_point(color="gray", aes(year, geom_mean_nonzero)) +
       	geom_line(color ="black", aes(year, geom_mean_nonzero)) +
	scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y), labels = exp_label) +
	labs(title = "Série chronologique moyenne géométrique positive",
		x = "Année",
		y = "Moyenne géométrique positive (µg/m³ (25°C))",
		caption = "Source : U.S. Environmental Protection Agency") +
	theme_classic() +
	facet_wrap(~location, ncol = 6)
dev.off()


#########################
# Geospatial Statistics #
#########################

#Load US states map
states <- ne_states(country = "United States of America", returnclass = "sf")

# Choose a specific year
year_of_interest <- 2019 - 6 

# Filter the dataset for the selected year
df_filtered <- df_yearly %>%
  filter(year == year_of_interest)

# Merge the geometric means data with the states map
geo_data <- merge(states, df_filtered, by.x = "name", by.y = "location")

# Plotting
CairoPNG(paste0(DATASETS, "/ArsenicPM2.5-CRD/images/descriptive_statistics/univariate/ArsenicMap.png"), width=1200, height=1000)
ggplot(data = geo_data) +
  geom_sf(aes(fill = geom_mean_nonzero), color = "black") +
  scale_fill_viridis_c(option = "plasma", name = "Arsenic PM2,5\nµg/m³ (25°C)   ") +
  labs(title = paste("Moyenne géométrique positive par État pendant l'année", year_of_interest),
       subtitle = "Source: U.S. Environmental Protection Agency") +
  coord_sf(xlim = c(-125, -66.5), ylim = c(24, 49.5)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(3, "line"),
        legend.title = element_text(margin = margin(b = 20)),
   legend.justification.bottom = "right"
        
 )
dev.off()




