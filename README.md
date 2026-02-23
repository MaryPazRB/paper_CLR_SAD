# Coffee Leaf Rust SAD Analysis: Reproducibility Guide

This repository contains the full analysis workflow for evaluating Standard Area Diagrams (SADs) used to assess Coffee Leaf Rust (CLR) severity. The analysis compares an "Old" SAD versus a "New" SAD, evaluating rater performance, assessment time, and sampling efficiency.

## 📋 Project Overview

The analysis is divided into four main stages:
1.  **Lab Validation**: Performance assessment using standardized leaf images in a controlled setting.
2.  **Field Validation**: Accuracy and precision comparison against a "Gold Standard" (GS) in field conditions.
3.  **Time Analysis**: Modeling the time-to-assess per branch using different SADs.
4.  **Sampling Efficiency**: Determining optimal sampling designs and simulating seasonal time savings.

## ⚙️ Dependencies

Ensure you have R installed along with the following packages:

```r
# Data Wrangling & Utilities
install.packages(c("gsheet", "dplyr", "tidyr", "readr", "purrr", "forcats", "writexl"))

# Statistical Modeling
install.packages(c("lme4", "lmerTest", "emmeans", "DHARMa", "epiR", "irr", "psych", "pbkrtest"))

# Visualization
install.packages(c("ggplot2", "patchwork", "cowplot", "scales", "viridis", "ggridges", "ggthemes"))
```

## 🚀 Execution Sequence

To reproduce the findings, run the scripts in the following order:

### 1. Laboratory Validation (`lab_validation.R`)
*   **Purpose**: Analyzes rater concordance (CCC) and coverage probability (using a hybrid tolerance rule) for laboratory-based ratings.
*   **Data Source**: Fetched via `gsheet` from the Lab Ratings dataset.
*   **Key Metrics**: Agreement, Precision, Bias coefficient, and Interrater Reliability (ICC).
*   **Main Outputs**: `panel_2x2.png`, `panel_est_actual.png`.

### 2. Field Validation (`field_validation.R`)
*   **Purpose**: Validates rater accuracy against a Gold Standard for field-collected samples.
*   **Data Source**: Fethed via `gsheet` from the Field Ratings dataset.
*   **Process**: Filters human raters, removes problematic images, and calculates accuracy components.
*   **Intermediate Data**: Produces the `df_ratings2` dataset used in subsequent sampling analyses. (Ensure this is exported as `df_ratings2.csv` for the final step).
*   **Main Outputs**: `plots_field.png`.

### 3. Assessment Time Analysis (`field_time.R`)
*   **Purpose**: compares the speed of assessment between Old and New SADs across three different fields.
*   **Process**: Fits a Gamma Generalized Linear Mixed Model (GLMM) with a log link to account for hierarchical data (Field > Evaluator > Plant > Branch).
*   **Main Outputs**: 
    *   `field_all.csv`: Cleaned time data used for subsequent simulations.
    *   `figs/fig_time_branch.png` and `figs/figura_plots.png`.

### 4. Precision & Seasonal Simulations (`field_time_minimum_sample.R`)
*   **Purpose**: Integrates performance and time data to optimize field protocols.
*   **Dependencies**: Requires `field_all.csv` and `df_ratings2.csv`.
*   **Analysis**:
    *   Calculates variance components (Plant, Branch, Leaf).
    *   Determines the minimum number of plants sampled per field to reach a target precision (95% CI ± 2.0 p.p.).
    *   Simulates "Seasonal Time Savings" across 500 fields under Mild, Typical, and Severe disease scenarios.
*   **Main Outputs**: 
    *   `minimaall.csv`: Summary of minimum sampling requirements.
    *   `figs/fig_min_plants_time.png`: Visualization of sampling designs and predicted assessment times.

## 📂 Data Sources

Data is hosted on Google Sheets and accessed programmatically:
- **Lab Ratings**: [Google Sheet ID: 141vQED6DLbzC0tQfvz16r_qR-UeQ0Hd7VOU46fJ_Z-Q]
- **Field Ratings**: [Google Sheet ID: 1rVFraYtTIUxoIfk5w7F4LdoniNPRXu6R]
- **Field Time**: [Google Sheet ID: 1_fO1nLXZzxPsKPBDb4JY3d76Ddthaf3bbTo7Xw76sVY]

## 📁 Output Directory
Most figures are saved within the `figs/` directory. Ensure this directory exists before running the scripts or the scripts will attempt to create it.


## 📷 About the datasets.

Datasets are disponibilized for usage in Kagglehub in the following link:

https://www.kaggle.com/datasets/marypazrb/paper-clr-sad

Annotated datasets for golden standard is available in roboflow universe:

 - **40 Images used for Lab validation**

https://universe.roboflow.com/clr-zky50/40_img_lab/dataset/1

 - **100 Images used for field validation**

https://universe.roboflow.com/clr-zky50/imgtest-fvn9j/dataset/1
