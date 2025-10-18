# opImputation: a model‑agnostic framework for dataset‑specific missing value imputation

**opImputation** is an R package for comparing, benchmarking, and applying missing value imputation strategies to numerical tabular data.  
It is primarily developed for biomedical and clinical research but is broadly applicable to any numerical dataset with missing values. This repository contains the source code and accompanies:

> Lötsch J, Ultsch A. _A model‑agnostic framework for dataset‑specific selection of missing value imputation methods in pain‑related numerical data._ Can J Pain (2025, in minor revision)

---

## Features

- **Model-agnostic benchmarking:** Compare imputation methods from diverse theoretical families  
- **Dataset-specific selection:** Identify the best-performing method for your actual data  
- **Diagnostic imputations:** Includes reference models for objective bias/error evaluation  
- **Parallel processing:** Efficient computation for large datasets  
- **Reproducible results:** Standardized workflows and seedable randomization  
- **Extensible integration:** Easily add your own imputation methods or reference results

---

## Installation
```r 
if (!requireNamespace("devtools")) install.packages("devtools") devtools::install_github("JornLotsch/opImputation")
``` 

**Requirements:**  
- R ≥ 3.5.0  
- Imports: parallel, Rfit, methods, stats, caret, ABCanalysis, ggplot2, pbmcapply,
  missForest, utils, mice, miceRanger, multiUS, Amelia, mi, reshape2,
  DataVisualizations, cowplot, twosamples, ggh4x, ggrepel, doParallel, foreach,
  tools, Rcpp

---

## Usage

### Basic example
```r 
library(opImputation)
# Compare imputation methods for your dataset
# Always set a seed for deterministic, fast results (strongly recommended!)
results <- compare_imputation_methods(
  Data = iris[,1:4],
  ImputationMethods = c("rf_missForest", "median", "plus"),
  nIter = 5,
  Seed = 42
)
# Impute missing data using the selected best method
imputed <- imputeData(
  Data = iris[,1:4],
  ImputationMethod = "rf_missForest",
  Seed = 42
)
``` 

---

## Main functions

### compare_imputation_methods

| Argument | Description |
|-----------|--------------|
| `Data` | Numeric data frame or matrix. All columns must be numeric. |
| `ImputationMethods` | Vector of method names to test. Use `all_imputation_methods` for a full list. |
| `ImputationRepetitions` | Number of repeats per method (default: 20). |
| `Seed` | Integer seed for reproducibility.<br>**For fast and reproducible results, always explicitly set `Seed` (e.g. `Seed = 42`).**<br>If not set, attempts slow RNG state recovery. |
| `nIter` | Number of iterations for diagnostic missing value insertions (default: 20). |
| `nProc` | Number of CPU cores for parallel processing (default: `getOption("mc.cores", 2L)`). |
| `probMissing` | Proportion (0–1) of values to remove for diagnostics (default: 0.1). |
| `PValueThresholdForMetrics` | Threshold for considering bias/error non‑significant (default: 0.1). |
| `pfctMtdsInABC` | Include calibration methods in ABC summaries (logical, default: FALSE). |
| `mnarity` | Fraction of MNAR missingness in diagnostics (default: 0). |
| `lowOnly` | Restrict diagnostic missings to lowest values (default: FALSE). |
| `mnarshape` | Shape parameter for MNAR curve (default: 1). |
| `test_only_variables_with_missings` | Analyze only columns with missings (default: FALSE). |
| `PlotIt` | Generate diagnostic summary plots (default: TRUE). |
| `overallBestzDelta` | Compare results to overall best method (default: FALSE). |

**Returns:**  
- `RepeatedSampleImputations` — all individual results  
- `zDeltas` — error and deviation metrics  
- `MethodsResults` — summary statistics and rankings  
- `BestMethodPerDataset` — optimal method per dataset  
- `Fig_zDeltaDistributions_bestMethods` — visual comparison of top methods  
- `Fig_compare_imputation_methods` — overall summary plot  

---

### imputeData

| Argument | Description |
|-----------|--------------|
| `Data` | Numeric data frame or matrix to be imputed. |
| `ImputationMethod` | Method name (e.g., `"median"`, `"rf_missForest"`). |
| `ImputationRepetitions` | Number of repeated imputations (default: 20 for repeated methods). |
| `Seed` | Integer seed for reproducibility.<br>**For fast and reproducible results, always explicitly set `Seed` (e.g. `Seed = 42`).**<br>If not set, attempts slow RNG state recovery. |
| `nProc` | Number of CPU cores for parallel operations (default: `getOption("mc.cores", 2L)`). |
**Returns:**  
A data frame or matrix of identical dimensions with missing values imputed.

---

## Output and diagnostics

The main product of **opImputation** is a detailed summary of imputation benchmarking, including:

- **Standardized mean ranks** for all candidate methods (ABC category overlays)
- **Mean absolute standardized errors** from diagnostics
- **Per-variable error distributions** across all imputations

These summaries enable transparent, reproducible, and data‑driven imputation selection.

### Summary plot

*Example diagnostic summary from the Iris dataset: mean rank bar plots, ABC curves, and error diagnostics for multiple imputation models.*

<img src="./Iris_TestOutput_annotated.svg" width="750">

A: Bar graph of standardized mean ranks for all imputation methods, with cABC category coloring and an overlaid ABC curve.  
B: Mean absolute standardized errors for different models using inserted diagnostic missing values.  
C: Mean absolute standardized error per variable and model across all iterations.  
Sets A, B, and C represent best, next-best, and discouraged models. See manuscript for further details.

---

## When to use opImputation

- Biomedical/clinical studies with complex or substantial missing data
- Preprocessing for multivariate analysis or machine learning
- Transparent benchmarking and method selection
- Teaching, methods development, or software validation in imputation

---

## Citation

If you use **opImputation**, please cite:

> Lötsch J, Ultsch A. (2025). _A model‑agnostic framework for dataset‑specific selection of missing value imputation methods in pain‑related numerical data._ Can J Pain (in minor revision)

---

## Authors and license

- Jorn Lotsch (author, maintainer)  
- Alfred Ultsch (author)  
- License: GPL‑3

---

## About this project

**opImputation** provides a robust, reproducible, and extensible framework for dataset-specific comparison, selection, and benchmarking of missing value imputation strategies.  
It enables transparent analysis, integration of new methods, and facilitates reliable, data-driven data preprocessing.
```
