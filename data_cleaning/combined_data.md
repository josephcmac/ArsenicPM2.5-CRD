# R Script Documentation

## Name
Environmental and Health Data Processing

## Synopsis
This script processes environmental data related to Arsenic PM2.5 levels and health incidence data to create a combined dataset for analyzing the relationship between environmental pollutants and health outcomes.

## Introduction
This documentation provides a detailed description of the functions and their usage in the R script that processes environmental and health data. The script utilizes both R and C++ functions for specific computations to optimize performance.

## Requirements
- `tidyverse` package
- `Rcpp` package

## Functions

### logit100000
`logit100000[x]`
- **Description**: Computes the logit transformation for a value `x` based on a fixed total of 100,000.
- **Usage**: 
  ```R
  logit100000(x)
  ```
- **Arguments**:
  - `x`: A numeric value.
- **Returns**: The logit transformation of `x`.

### log_geom_mean_nonzero
`log_geom_mean_nonzero[x]`
- **Description**: Calculates the logarithm of the geometric mean of a numeric vector, excluding zeros.
- **Usage**:
  ```R
  log_geom_mean_nonzero(x)
  ```
- **Arguments**:
  - `x`: A numeric vector.
- **Returns**: The logarithm of the geometric mean of the non-zero elements in `x`.

### fix_age_single
`fix_age_single[x]`
- **Description**: Converts age groups from string representations to integer codes.
- **Usage**:
  ```R
  fix_age_single(x)
  ```
- **Arguments**:
  - `x`: A string representing an age group.
- **Returns**: An integer code corresponding to the age group.

### convert_age_to_group
`convert_age_to_group[x]`
- **Description**: Applies the `fix_age_single` function to a vector of age groups, converting them to integer codes.
- **Usage**:
  ```R
  convert_age_to_group(x)
  ```
- **Arguments**:
  - `x`: A vector of strings representing age groups.
- **Returns**: A vector of integer codes corresponding to the age groups.

### separate_by_sex
`separate_by_sex[df]`
- **Description**: Separates a data frame by sex and applies the `logit100000` transformation to the incidence values.
- **Usage**:
  ```R
  separate_by_sex(df)
  ```
- **Arguments**:
  - `df`: A data frame containing health incidence data with columns for `sex` and `incidence`.
- **Returns**: A data frame with separate columns for male and female logit-transformed incidence values.

### combine_environment_and_health_data
`combine_environment_and_health_data[df_env_yearly, df_health, latency]`
- **Description**: Combines environmental data and health data with a specified latency period.
- **Usage**:
  ```R
  combine_environment_and_health_data(df_env_yearly, df_health, latency)
  ```
- **Arguments**:
  - `df_env_yearly`: A data frame containing yearly environmental data.
  - `df_health`: A data frame containing health incidence data.
  - `latency`: An integer representing the latency period.
- **Returns**: A combined data frame with adjusted environmental data and separated health data by sex.

### read_incidence_file
`read_incidence_file[filename, DIR]`
- **Description**: Reads a single incidence file and selects specific columns for further processing.
- **Usage**:
  ```R
  read_incidence_file(filename, DIR)
  ```
- **Arguments**:
  - `filename`: A string representing the file name of the incidence data.
  - `DIR`: A string representing the directory containing the incidence data files.
- **Returns**: A data frame containing selected columns from the incidence file.

### read_incidence_files
`read_incidence_files[filenames, DIR]`
- **Description**: Reads multiple incidence files and processes them.
- **Usage**:
  ```R
  read_incidence_files(filenames, DIR)
  ```
- **Arguments**:
  - `filenames`: A vector of strings representing the file names of the incidence data.
  - `DIR`: A string representing the directory containing the incidence data files.
- **Returns**: A combined data frame containing processed incidence data from all specified files.

### read_env_year
`read_env_year[year, parameter, sample_duration, DIR]`
- **Description**: Reads environmental data for a specific year, parameter, and sample duration.
- **Usage**:
  ```R
  read_env_year(year, parameter, sample_duration, DIR)
  ```
- **Arguments**:
  - `year`: An integer representing the year of the environmental data.
  - `parameter`: A string representing the parameter of interest (e.g., "Arsenic PM2.5 LC").
  - `sample_duration`: A string representing the sample duration (e.g., "24 HOUR").
  - `DIR`: A string representing the directory containing the environmental data files.
- **Returns**: A processed data frame containing environmental data for the specified year.

### read_data
`read_data[filenames, latency, year0, year1, DIR_env, DIR_health]`
- **Description**: Reads and processes environmental and health data over a range of years, and returns a combined dataset.
- **Usage**:
  ```R
  read_data(filenames, latency, year0, year1, DIR_env, DIR_health)
  ```
- **Arguments**:
  - `filenames`: A vector of strings representing the file names of the health incidence data.
  - `latency`: An integer representing the latency period.
  - `year0`: An integer representing the starting year for the environmental data.
  - `year1`: An integer representing the ending year for the environmental data.
  - `DIR_env`: A string representing the directory containing the environmental data files.
  - `DIR_health`: A string representing the directory containing the health incidence data files.
- **Returns**: A combined data frame containing processed environmental and health data.

## Examples
### Example 1: Basic Usage
```R
# Load necessary libraries
library(tidyverse)
library(Rcpp)

# Call read_data function with appropriate parameters
combined_data <- read_data(
  filenames = c("file1", "file2"),
  latency = 6,
  year0 = 1988,
  year1 = 2019,
  DIR_env = "path/to/env/data",
  DIR_health = "path/to/health/data"
)

# Write the combined data to a CSV file
write.csv(combined_data, "path/to/output/combined.csv", row.names = FALSE)
```

## Properties & Relations
- The `logit100000` function is used within `separate_by_sex` to transform incidence values.
- The `fix_age_single` function is used within `convert_age_to_group` to convert age groups.
- The `log_geom_mean_nonzero` function is used within `read_env_year` to process environmental data.
- The `combine_environment_and_health_data` function integrates the processed environmental and health data considering latency.

## Applications
- Analyzing the impact of environmental pollutants on health outcomes.
- Studying the relationship between Arsenic PM2.5 levels and health incidence rates across different age groups and sexes.

## See Also
- `tidyverse` package documentation
- `Rcpp` package documentation

## Tutorials
- Using C++ with R: Introduction to Rcpp
- Data Manipulation with dplyr and tidyr

## Implementation Notes
- Ensure that the necessary packages are installed and loaded.
- Update directory paths and file names as per your dataset structure before running the script.
