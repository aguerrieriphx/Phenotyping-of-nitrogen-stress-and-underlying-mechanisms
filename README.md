# Power Analysis for Digital Biomass, NDVI, and NPCI Phenotyping

This repository contains the official dataset and R code used for the power analysis presented in our paper: **"Efficient and accurate quantification of the effects of nitrogen stress on spectral and morphological traits, and their coordination with photosynthetic traits for underlying mechanisms in maize"**. 

The analysis uses a simulation-based approach powered by Quantile Regression to model biological variance and calculate:
1. **Sensitivity / Minimum Detectable Effect (MDE):** The percentage difference detectable with a fixed sample size (n = 4).
2. **Sample Size (N) Requirements:** The number of plant replicates needed to detect small-to-moderate biological changes (10%, 15%, and 20%) at 80% statistical power.

---

## Repository Structure

* `Merged_dataset.csv` – The clean, curated phenotypic dataset used for modeling.
* `Publication script for power analysis v2.R` – The execution script that runs the noise modeling, Monte Carlo power simulations, and generates plots.
* `README.md` – This documentation file.

---

## Prerequisites & Installation

To run the analysis script, you will need R (version 4.0 or higher recommended) and a few essential packages. You can install all dependencies by running the following command in your R console:

```R
install.packages(c("dplyr", "tidyr", "ggplot2", "gridExtra", "quantreg"))
```

## Summary of Package Roles

* **quantreg**: Fits the 10th, 50th (median), and 84th percentiles to extract robust noise standard deviations.
* **dplyr & tidyr**: Data cleaning and structural manipulation.
* **ggplot2 & gridExtra**: Generating and aligning publication-ready multi-panel plots.

---

## How to Run the Analysis

1. **Clone or Download** this repository to your local machine.
2. Ensure both `Merged_dataset.csv` and `Publication script for power analysis v2.R` are placed in the same folder (working directory).
3. Open your R environment or RStudio, set your working directory to this folder, and run the script:

```R
source("Publication script for power analysis v2.R")
```

> **Note on Runtime:** The script runs 400 Monte Carlo simulations across 30 iterative loops per trait to ensure highly accurate, bounded estimations. Consequently, execution may take a few minutes to complete.

---

## Methodology Overview

Instead of assuming constant variance (homoscedasticity), this pipeline uses **Quantile Regression** to dynamically model noise based on biological variance across trait values:

* **Noise Modeling:** Fits `tau = 0.5` (median), `tau = 0.16` (lower bounds), and `tau = 0.84` (upper bounds) to capture the exact distribution of noise relative to the reference measurements (`biomass_g` and `avg_Chla`).
* **Power Simulation:** Simulates experimental control and shifted groups over a range of sample sizes and effect sizes using a standard t-test (`alpha = 0.05`).
* **Outlier Filtering:** Automatically drops sample `NN-4` during the initial data-cleaning step to preserve model robustness.

---

## Generated Outputs

Once execution is complete, two new validation files will appear in your folder:

* **`Power_Analysis_Results.csv`**: A summary table containing the Mean, Minimum, and Maximum values across all 30 simulation runs for:
  * **`MDE_Mean`**: Minimum detectable effect size at `n = 4`.
  * **`N10 / N15 / N20`**: Exact number of plants required to confidently detect a 10%, 15%, or 20% change.
* **`Power_Analysis_Plot.png`**: A high-resolution (300 DPI) 4-panel publication-ready visualization plotting the sensitivity margins and sample size thresholds for Digital Biomass, NDVI, and NPCI.
