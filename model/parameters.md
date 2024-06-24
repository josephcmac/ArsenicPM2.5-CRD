
# Documentation for Model Simulation and Parameter Estimation in R

## Overview

This R script provides a comprehensive framework for simulating data based on a mixture model, predicting outcomes, estimating initial model parameters using k-means clustering, computing model parameters using the Mclust package, and creating tables of parameters for different age and sex groups. It relies on multiple libraries, including `tidyverse`, `mclust`, `flexclust`, and `mvtnorm`.

## Libraries

The following libraries are used in this script:
- `tidyverse`: For data manipulation and visualization.
- `mclust`: For model-based clustering.
- `flexclust`: For flexible clustering algorithms.
- `mvtnorm`: For multivariate normal distributions.

## Functions

### Model Simulation

#### `simulate_model(params, n)`

Simulates data from a mixture model.

- **Parameters**:
  - `params`: A list containing the parameters of the mixture model (`p`, `mu`, `sigma`).
  - `n`: Number of samples to simulate.

- **Returns**:
  - A matrix of simulated data points.

### Model Prediction

#### `predict_model(params, x)`

Predicts the outcome based on the mixture model parameters and input `x`.

- **Parameters**:
  - `params`: A list of mixture model parameters.
  - `x`: The input variable.

- **Returns**:
  - Predicted value based on the mixture model.

### Initial Parameter Estimation

#### `estimate_init_model(data, n_clusters = 3)`

Estimates initial model parameters using k-means clustering.

- **Parameters**:
  - `data`: The dataset.
  - `n_clusters`: Number of clusters (default is 3).

- **Returns**:
  - A list containing initial means and variances for each cluster.

### Model Parameter Computation

#### `compute_params_model(data, n_clusters = 3)`

Computes the model parameters using the `Mclust` package.

- **Parameters**:
  - `data`: The dataset.
  - `n_clusters`: Number of clusters (default is 3).

- **Returns**:
  - A list of computed model parameters (`p`, `mu`, `sigma`, `loglik`).

#### `compute_params_model_list(data)`

Formats the computed parameters into a user-friendly list.

- **Parameters**:
  - `data`: The dataset.

- **Returns**:
  - A list of formatted model parameters.

### Table Creation

#### `filter_age_sex(df, age_min0, sex_male)`

Filters the dataset based on age and sex.

- **Parameters**:
  - `df`: The dataset.
  - `age_min0`: The minimum age for filtering.
  - `sex_male`: Boolean indicating if the data should be filtered for males.

- **Returns**:
  - Filtered dataset.

#### `create_table_age_sex(df, age_min, sex_male, compute_params)`

Creates a table of parameters for a specific age group and sex.

- **Parameters**:
  - `df`: The dataset.
  - `age_min`: The minimum age for the age group.
  - `sex_male`: Boolean indicating if the data should be filtered for males.
  - `compute_params`: Function to compute parameters.

- **Returns**:
  - A data frame with parameters for the specified age group and sex.

#### `create_table_age(df, age_min, compute_params)`

Creates a table of parameters for a specific age group.

- **Parameters**:
  - `df`: The dataset.
  - `age_min`: The minimum age for the age group.
  - `compute_params`: Function to compute parameters.

- **Returns**:
  - A data frame with parameters for the specified age group.

#### `create_table(df, compute_params)`

Creates a comprehensive table of parameters for all age groups and both sexes.

- **Parameters**:
  - `df`: The dataset.
  - `compute_params`: Function to compute parameters.

- **Returns**:
  - A data frame with parameters for all age groups and both sexes.

## Execution

The script reads a dataset, computes the parameters for each age and sex group, and writes the resulting table to a CSV file.

```r
DATASETS <- "../../../datasets/ArsenicPM2.5-CRD"

create_table(
  df = read.csv(paste0(DATASETS, "/combined.csv")), 
  compute_params = compute_params_model_list
) %>%
  write.csv(paste0(DATASETS, "/parameters.csv"), row.names=FALSE)
```

## Data Source

The dataset used is located at `../../../datasets/ArsenicPM2.5-CRD/combined.csv`, and the output is saved to `../../../datasets/ArsenicPM2.5-CRD/parameters.csv`.
