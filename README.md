# MSc Thesis: Decomposing health inequalities through regression trees

2026-06-16

## Overview

Health inequalities are often studied using linear regression–based
decompositions, but these impose restrictive functional-form and
additivity assumptions that may obscure important heterogeneity. This
project develops and evaluates a non-parametric, regression-tree–based
decomposition method that pairs classification and regression trees with
standard inequality indices.

The goal is to:

1.  Minimise within-group and/or between-group disparities.
2.  Reveal how determinants jointly segment the population into
    subgroups relevant to inequality.

## Research Questions

1.  Does tree-based partitioning reveal more policy-relevant subgroup
    structures than linear decompositions under non-linear or
    interactive relationships?
2.  How do decomposition results vary across different outcome
    definitions (e.g., concentration index–based measures)?
3.  What is the marginal contribution of each determinant in a
    recursively partitioned model?

## Data

The method will be applied to demographic surveillance datasets that
include socioeconomic status, health outcomes, and a range of
covariates.

Methodological considerations include algorithm implementation and
accounting for survey design features such as weights and
stratification.

## Repository Structure

- [**R /**](./R) - R scripts and functions
- [**Report/**](./report) - methodology, results, and discussion
  write-up
- [**Simulation/**](./simulation) - code and results for simulation
  studies
- [**MSE
  Report/**](./report/MSE_Thesis_Decomposition_of_CI_Trees_MM.qmd) - QMD
  used to generate thesis report
- [**DRC Results
  Generator/**](./report/Generate_DRC_Results_Objects.qmd) - QMD used
  to generate the main analysis objects
- [**Design-Aware
  Validation/**](./report/design_aware_validation.qmd) - QMD used for
  PSU-level bootstrap stability checks
- [**Inequality Tree
  Implementation**](https://github.com/m-mburu/ineqTrees) - R package
  for concentration index tree and forest implementation

## Reproducible Analysis

Open the repository from the project root, either with the
[`thesis-mse.Rproj`](./thesis-mse.Rproj) file in RStudio or as the
project folder in VS Code. The analysis uses `here()`, so running code
from the project root keeps file paths consistent.

Most dependencies can be installed from CRAN. Two important packages
should be installed from GitHub:

``` r
install.packages("remotes")
remotes::install_github("m-mburu/ineqTrees")
remotes::install_github("bgreenwell/fastshap")
```

`ineqTrees` contains the concentration-index tree and forest
implementation used in this thesis. `fastshap` is used for the
predictive-forest SHAP decomposition and has been removed from CRAN, so
it should be installed from GitHub.

To reproduce the thesis results, run the analysis files in this order:

1.  Render [**Simulation Fits**](./simulation/simulation_fits.qmd). This
    builds the null and poor-subgroup simulation examples and saves the
    simulation tree plots used in the report.
2.  Render [**DRC Results
    Generator**](./report/Generate_DRC_Results_Objects.qmd). This
    creates `data/drc_report_results_objects.rda`, which is the main
    object consumed by the thesis report. It includes the DRC descriptive
    summaries, survey-weighted regression and `rineq` decompositions,
    concentration-index trees for `CI`, `CIg`, `CIc`, and `L`,
    lecturer/rpart comparison trees, concentration-index forests and
    surrogate trees, the predictive `ranger` forest, `fastshap` SHAP
    decompositions, variable-importance objects, and the method-comparison
    tables.
3.  Render [**Design-Aware
    Validation**](./report/design_aware_validation.qmd). This uses the
    generated DRC results object and runs PSU-level bootstrap stability
    checks for the selected concentration-index trees.
4.  Render [**MSE
    Report**](./report/MSE_Thesis_Decomposition_of_CI_Trees_MM.qmd). This
    is the final thesis report QMD. It reads the saved DRC results,
    simulation plots, and bootstrap validation objects produced by the
    previous steps.

From a terminal, the same order is:

``` sh
quarto render simulation/simulation_fits.qmd
quarto render report/Generate_DRC_Results_Objects.qmd
quarto render report/design_aware_validation.qmd
quarto render report/MSE_Thesis_Decomposition_of_CI_Trees_MM.qmd
```

## Data Availability Note

The source data are not shared in this GitHub repository. Public GitHub
users can inspect the code and report structure, but the full analysis
will not run unless the private data files are supplied separately and
placed in the expected project folders. My supervisors will be able to
reproduce the workflow after I share the data with them separately.
