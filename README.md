# opImputation: A model-agnostic framework for dataset-specific selection of missing value imputation methods in pain-related numerical data 

**opImputation** is an R package for comparing, benchmarking, and applying missing‑value imputation strategies to numerical tabular data. 
It is developed for biomedical and clinical research but is broadly applicable to any numerical dataset containing missing values.

---

## Features

- **Model‑agnostic benchmarking:** Compare imputation methods from diverse algorithmic families  
- **Dataset‑specific selection:** Automatically identify the best‑performing method for your data  
- **Automated imputation:** Optionally produce a final imputed dataset using the top method  
- **Parallel processing:** Efficient computation using the *future* framework (`future.apply`, `progressr`)  
- **Reproducible analysis:** Seedable, standardized workflows  
- **Extensible integration:** Add new methods or external benchmarking data easily  

---

## Installation
```
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("JornLotsch/opImputation")
```

**Package metadata:**  
| Field | Value |
|-------|--------|
| **Type** | R Package |
| **Title** | Optimal Selection of Imputation Methods for Bio‑Medical Data |
| **Version** | 0.4 |
| **Depends** | R (≥ 3.5.0) |
| **Imports** | parallel, Rfit, methods, stats, caret, ABCanalysis, ggplot2, future.apply, progressr, missForest, utils, mice, miceRanger, multiUS, Amelia, mi, reshape2, DataVisualizations, cowplot, twosamples, ggh4x, ggrepel, tools, Rcpp (≥ 1.0.0) |
| **LinkingTo** | Rcpp |
| **License** | GPL‑3 |
| **Authors** | Jörn Lotsch ([ORCID 0000‑0002‑5818‑6958]), Alfred Ultsch ([ORCID 0000‑0002‑7845‑3283]) |
| **Maintainer** | Jörn Lotsch <j.lotsch@em.uni‑frankfurt.de> |
| **Repository** | [https://github.com/JornLotsch/opImputation](https://github.com/JornLotsch/opImputation) |
| **Date** | 2025‑05‑03 |
| **NeedsCompilation** | yes |

---

## Usage

### Basic example
```
library(opImputation)

# Load example dataset (numeric columns only)
data <- iris[, 1:4]

# Introduce 5% random missing values, avoiding empty rows
set.seed(42)
n_total <- prod(dim(data))
n_missing <- round(0.05 * n_total)
repeat {
  tmp <- data
  pos <- sample(n_total, n_missing)
  tmp[pos] <- NA
  if (all(rowSums(is.na(tmp)) < ncol(tmp))) break
}
data <- tmp

# Step 1: Compare imputation methods
results <- compare_imputation_methods(
  data = data,
  imputation_methods = c("rf_mice", "median", "knn5"),
  imputation_repetitions = 20,
  n_iterations = 10,
  percent_missing = 0.05,
  seed = 123,
  n_proc = 2,
  plot_results = TRUE
)

# Step 2: Retrieve automatically generated final imputation
imputed_data <- results$imputed_data
print(results$method_used_for_imputation)
```
---

## Main functions

### compare_imputation_methods

| Argument | Description |
|-----------|-------------|
| `data` | Numeric data frame or matrix. May contain existing missing values. |
| `imputation_methods` | Character vector of imputation method names to compare. Default: `all_imputation_methods`. Must include at least two non‑calibrating methods. |
| `imputation_repetitions` | Integer. Number of repeated imputations for each method and iteration (default = 20). |
| `perfect_methods_in_ABC` | Logical. If `TRUE`, calibration methods are included in the final categorization of methods (default value = FALSE). For testing purposes only; do not set to TRUE in real test environments.|
| `n_iterations` | Number of missing data patterns to test (default = 20). |
| `n_proc` | Number of CPU cores for parallel processing (default: `getOption("mc.cores", 2L)`). |
| `percent_missing` | Numeric. Proportion of data to randomly set missing (0‑1; default = 0.1). |
| `seed` | Integer. Random seed for reproducibility (recommended). |
| `mnar_shape` | Shape parameter for the *Missing Not At Random* (MNAR) mechanism (default = 1). |
| `mnar_ity` | Degree of MNAR dependency (0–1; default = 0 → MCAR). |
| `low_only` | Logical. If `TRUE`, insert missings only in lower‑valued observations. |
| `fixed_seed_for_inserted_missings` | Logical. Repeat identical random pattern across iterations. |
| `max_attempts` | Maximum attempts to avoid creating empty rows (default = 1000). |
| `plot_results` | Logical. If `TRUE`, create summary plots (default = TRUE). |
| `overall_best_z_delta` | Logical. Compare to global best or category best method (default = FALSE). |
| `produce_final_imputations` | Logical. If `TRUE`, generates final imputed dataset using the best‑ranked valid method (default = TRUE). |

**Returns:**  
- `repeated_sample_imputations` — list of all iteration results  
- `z_deltas` — standardized z‑delta metrics (raw values, medians, and row medians)  
- `methods_results` — ABC analysis metrics and method rankings  
- `best_method_per_dataset` — name of the overall best method  
- `best_univariate_method`, `best_multivariate_method`, `best_uni_multivariate_method`, `best_poisoned_method` — category bests  
- `df_abc_results` — results of ABC analysis (category and score)  
- `fig_z_delta_distributions_best_methods` — comparison of z‑delta distributions  
- `fig_comparison_summary` — combined summary figure (ABC + z‑delta plots)  
- `imputed_data` — final imputed dataset (if `produce_final_imputations = TRUE`)  
- `method_used_for_imputation` — name of the actual method used  

---

### impute_missings

| Argument | Description |
|-----------|-------------|
| `x` | Numeric data frame or matrix with missing values. |
| `method` | Imputation method name (default = `"rf_missForest"`). |
| `ImputationRepetitions` | Number of repetitions for methods ending with `_repeated` (default = 10). |
| `seed` | Random seed for reproducibility (recommended). |
| `x_orig` | Original dataset required for “poisoned” or “calibrating” methods. |

**Returns:**  
A numeric data frame of the same dimensions and column names, with all missing values imputed.  

---

## Output and diagnostics

Performance evaluation is based on the standardized **Δz (z‑delta)** metric—  
a robust measure of the absolute deviation between true and imputed values.  
**ABC (Activity‑Based Classification)** categorizes imputation methods by their relative performance,  
highlighting “A‑class” models as top performers.

Example output table from `res_abc$df_abc_results[,1:3]`(generic dataset):

| abc_score | abc_category | method |
|-----------:|:-------------:|:-------|
| 36.5755 | A | plusminus |
| 21.8880 | A | cart_repeated |
| 17.5513 | A | pmm_repeated |
| 16.8750 | A | rf_mice_repeated |
| 16.0810 | A | miceRanger |
| 15.3061 | A | miceRanger_repeated |
| 9.7959 | A | cart |
| 9.0947 | B | pmm |
| 7.2345 | B | rf_missForest |
| 7.0602 | B | amelia_repeated |
| 4.8430 | B | miImp |
| 4.8430 | B | rf_mice |
| 3.2741 | C | plus |
| 3.2153 | C | rf_missForest_repeated |
| 1.9199 | C | amelia |
| 1.4161 | C | knn3 |
| 1.3021 | C | linear |
| 0.0000 | C | bag |
| 0.0000 | C | bag_repeated |
| 0.0000 | C | factor |
| 0.0000 | C | knn10 |
| 0.0000 | C | knn5 |
| 0.0000 | C | knn7 |
| 0.0000 | C | knn9 |
| 0.0000 | C | mean |
| 0.0000 | C | median |
| 0.0000 | C | mode |
| 0.0000 | C | rSample |

**Legend:**  
- `abc_score`: zDelta values. Quantitative measure of imputation performance (higher = better).  
- `abc_category`: ABC‑derived ranking class (“A” = top, “B” = medium, “C” = low).  
- `method`: Name of the evaluated imputation algorithm.  
“A‑class” methods (top seven in this example) represent the highest‑performing algorithms for the tested dataset.  
Lower tiers correspond to progressively weaker or calibration‑only approaches.


---

## Example summary plot

*Diagnostic summary from the Iris dataset: ABC curves and variable‑specific Δz distributions (different dataset than that used in above table).*

<img src="./Iris_TestOutput_annotated.svg" width="750">

A: Standardized mean ranks for all imputation methods with ABC category coloring.  
B: Mean standardized Δz deviations for diagnostic missings.  
C: Variable‑level Δz distributions across methods.  

---

## When to use opImputation

- Biomedical or clinical datasets with incomplete numerical data  
- Multivariate analysis or machine‑learning preprocessing  
- Benchmarking and transparent method selection  
- Fully automated, reproducible imputation pipelines  

---

## Citation

If you use **opImputation**, please cite:

> Lötsch J, Ultsch A. (2025).  
> *A model‑agnostic framework for dataset‑specific selection of missing value imputation methods in pain‑related numerical data.*  
> *Can J Pain* (in minor revision)

---

## Authors and license

- **Jörn Lötsch** (author, creator, maintainer)  
- **Alfred Ultsch** (author)  
- License: GPL‑3  

---

## About this project

**opImputation** provides an automated, transparent, and reproducible framework for dataset‑specific benchmarking and optimal selection of missing‑value imputation methods.  
The framework incorporates *Activity‑Based Classification (ABC)* and *computed ABC (cABC)* analyses to identify statistically top‑performing algorithms and to optionally generate a fully imputed dataset automatically.  

For theoretical background, see:  
- Ultsch A, Lötsch J. *Computed ABC Analysis for Rational Selection of Most Informative Variables in Multivariate Data.* *PLoS One.* 2015; 10(6): e0129767.  
- Lötsch J, Ultsch A. *Recursive computed ABC (cABC) Analysis for Reducing Machine‑Learning Feature Sets to Their Minimum Informative Size.* *Sci Rep.* 2023; 13(1): 5470.
```
