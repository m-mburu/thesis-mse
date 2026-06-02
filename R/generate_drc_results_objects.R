params <-
list(run_cv = FALSE, run_forest = FALSE, run_forest_cv = FALSE, 
    run_shap = FALSE, shap_n = 300L, shap_nsim = 8L, ranger_death_target_share = 0.2)

knitr::opts_chunk$set(
  fig.width = 9,
  fig.height = 5.5,
  dpi = 500,
  dev = "png"
)

library(data.table)
library(ggplot2)
library(grid)
library(knitr)
library(ineqTrees)
library(partykit)
library(here)
source(here("R", "report_helpers.R"))
source(here("R", "ci_known_subgroup_demo.R"))
source(here("R", "ci_rpart_comparison_methods.R"))

if (requireNamespace("kableExtra", quietly = TRUE)) {
  kableExtra::use_latex_packages()
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  tolower(value) %in% c("1", "true", "yes", "y")
}

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = NA_character_)))
  if (is.na(value)) default else value
}

env_num <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = NA_character_)))
  if (is.na(value)) default else value
}

default_params <- list(
  run_cv = FALSE,
  run_forest = FALSE,
  run_forest_cv = FALSE,
  run_shap = FALSE,
  shap_n = 300L,
  shap_nsim = 8L,
  ranger_death_target_share = 0.20
)

runtime_params <- if (!exists("params", inherits = FALSE)) {
  default_params
} else {
  utils::modifyList(default_params, as.list(params))
}

runtime_params$run_cv <- env_flag("RUN_CV", runtime_params$run_cv)
runtime_params$run_forest <- env_flag("RUN_FOREST", runtime_params$run_forest)
runtime_params$run_forest_cv <- env_flag("RUN_FOREST_CV", runtime_params$run_forest_cv)
runtime_params$run_shap <- env_flag("RUN_SHAP", runtime_params$run_shap)
runtime_params$shap_n <- env_int("SHAP_N", runtime_params$shap_n)
runtime_params$shap_nsim <- env_int("SHAP_NSIM", runtime_params$shap_nsim)
runtime_params$ranger_death_target_share <- env_num(
  "RANGER_DEATH_TARGET_SHARE",
  runtime_params$ranger_death_target_share
)

criterion_types <- c("CI", "CIg", "CIc", "L")
rineq_ci_types <- c("CI", "CIg", "CIc", "CIw")
results_object_file <- here("data", "drc_report_results_objects.rda")
outcome_labels <- c(
  alive = "Alive at age 5",
  died = "Died before age 5"
)

main_tree_control <-ci_tree_control(
  minsplit = 500L,
  minbucket = 250L,
  minprob = 0.02,
  maxdepth = 5L,
  min_gain = 0.00001,
  min_relative_gain = 0.20
)

map_model_terms_to_predictors <- function(term_names, predictors) {
  vapply(term_names, function(term) {
    hit <- predictors[term == predictors | startsWith(term, predictors)]
    if (length(hit)) hit[which.max(nchar(hit))] else term
  }, character(1L))
}

canonicalize_model_terms <- function(term_names, reference_terms) {
  if (!length(term_names) || !length(reference_terms)) {
    return(term_names)
  }

  lookup <- data.table(
    variable = reference_terms,
    syntactic = make.names(reference_terms, unique = FALSE)
  )
  lookup <- lookup[
    !duplicated(syntactic) & !duplicated(syntactic, fromLast = TRUE)
  ]

  out <- term_names
  exact <- match(out, reference_terms)
  out[!is.na(exact)] <- reference_terms[exact[!is.na(exact)]]

  needs_lookup <- which(is.na(exact))
  syntactic <- match(out[needs_lookup], lookup$syntactic)
  out[needs_lookup[!is.na(syntactic)]] <- lookup$variable[syntactic[!is.na(syntactic)]]

  out
}

rank_from_named_vector <- function(x, method_name) {
  if (is.null(x) || !length(x)) {
    return(data.table(variable = character(), method = method_name, rank = integer()))
  }
  data.table(variable = names(x), score = as.numeric(x))[
    is.finite(score) & score != 0
  ][
    order(-abs(score))
  ][
    ,
    .(variable, method = method_name, rank = seq_len(.N), score)
  ]
}

load(here("data", "congo_model_data.rda"))

raw_congo_model_dt <- as.data.table(congo_model_dt)
raw_congo_var_labels <- congo_var_labels
raw_congo_predictors <- names(raw_congo_var_labels)

predictor_level_counts <- vapply(
  raw_congo_model_dt[, raw_congo_predictors, with = FALSE],
  function(x) length(levels(as.factor(x))),
  integer(1L)
)
nominal_predictors <- raw_congo_predictors[
  vapply(
    raw_congo_model_dt[, raw_congo_predictors, with = FALSE],
    function(x) is.factor(x) || is.character(x),
    logical(1L)
  )
]
multi_level_predictors <- intersect(
  names(predictor_level_counts[predictor_level_counts > 2L]),
  nominal_predictors
)
kept_predictors <- setdiff(raw_congo_predictors, multi_level_predictors)

make_indicator_name <- function(variable, level) {
  make.names(paste(variable, level, sep = "_"), unique = FALSE)
}

indicator_level_map <- rbindlist(lapply(multi_level_predictors, function(variable) {
  levels_i <- levels(as.factor(raw_congo_model_dt[[variable]]))
  data.table(
    raw_variable = variable,
    level = levels_i,
    variable = make_indicator_name(variable, levels_i),
    variable_label = paste0(unname(raw_congo_var_labels[variable]), ": ", levels_i)
  )
}), fill = TRUE)

if (anyDuplicated(indicator_level_map$variable)) {
  stop("Indicator dummy names are not unique. Check predictor level names.", call. = FALSE)
}

model_predictor_dt <- copy(raw_congo_model_dt[, kept_predictors, with = FALSE])
for (i in seq_len(nrow(indicator_level_map))) {
  source_variable <- indicator_level_map$raw_variable[i]
  source_level <- indicator_level_map$level[i]
  target_variable <- indicator_level_map$variable[i]
  model_predictor_dt[
    ,
    (target_variable) := factor(
      fifelse(as.character(raw_congo_model_dt[[source_variable]]) == source_level, "Yes", "No"),
      levels = c("No", "Yes")
    )
  ]
}

congo_model_dt <- data.table(
  raw_congo_model_dt[, .(wealth, deadu5_num, sample_weight)],
  model_predictor_dt
)
congo_predictors <- names(model_predictor_dt)
congo_var_labels <- setNames(congo_predictors, congo_predictors)
congo_var_labels[intersect(kept_predictors, names(raw_congo_var_labels))] <-
  raw_congo_var_labels[intersect(kept_predictors, names(raw_congo_var_labels))]
congo_var_labels[indicator_level_map$variable] <- indicator_level_map$variable_label

congo_ci_formula <- stats::as.formula(
  paste("cbind(wealth, deadu5_num) ~", paste(congo_predictors, collapse = " + "))
)

canonicalize_regression_terms <- function(term_names) {
  vapply(term_names, function(term) {
    parent <- map_model_terms_to_predictors(term, raw_congo_predictors)
    if (parent %in% multi_level_predictors && parent != term) {
      source_level <- sub(paste0("^", parent), "", term)
      hit <- indicator_level_map[raw_variable == parent & level == source_level]
      if (nrow(hit)) {
        return(hit$variable[1L])
      }
    }
    if (parent %in% raw_congo_predictors && parent != term) {
      return(parent)
    }
    term
  }, character(1L))
}

format_variable_label <- function(x) {
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(value)) {
      return(NA_character_)
    }
    if (value %in% names(congo_var_labels)) {
      return(unname(congo_var_labels[value]))
    }
    if (exists("raw_congo_var_labels", inherits = TRUE) &&
        value %in% names(raw_congo_var_labels)) {
      return(unname(raw_congo_var_labels[value]))
    }

    parent <- map_model_terms_to_predictors(value, raw_congo_predictors)
    if (parent %in% names(raw_congo_var_labels) && parent != value) {
      suffix <- trimws(sub(paste0("^", parent), "", value))
      suffix <- gsub("^[:._]+", "", suffix)
      if (nzchar(suffix)) {
        return(paste0(unname(raw_congo_var_labels[parent]), ": ", suffix))
      }
      return(unname(raw_congo_var_labels[parent]))
    }

    value
  }, character(1L))
}



outcome_plot_dt <- congo_model_dt[
  ,
  .(weighted_n = sum(sample_weight)),
  by = .(deadu5_num)
][
  ,
  outcome_label := fifelse(deadu5_num == 1, outcome_labels[["died"]], outcome_labels[["alive"]])
][
  ,
  weighted_percent := 100 * weighted_n / sum(weighted_n)
]

candidate_predictor_table <- rbindlist(lapply(raw_congo_predictors, function(var) {
  x <- raw_congo_model_dt[[var]]
  level_dt <- raw_congo_model_dt[, .(weighted_n = sum(sample_weight)), by = var][order(-weighted_n)]
  level_dt[, weighted_percent := 100 * weighted_n / sum(weighted_n)]
  data.table(
    variable = var,
    label = unname(raw_congo_var_labels[var]),
    class = paste(class(x), collapse = "/"),
    n_levels = uniqueN(x),
    most_common_level = as.character(level_dt[[var]][1L]),
    most_common_percent = level_dt$weighted_percent[1L]
  )
}), fill = TRUE)

drc_root_ci_table <- rbindlist(lapply(criterion_types, function(type) {
  ci_fun <-ci_factory(type)
  data.table(
    criterion = type,
    estimate = ci_fun(
      cbind(congo_model_dt$wealth, congo_model_dt$deadu5_num),
      congo_model_dt$sample_weight
    )
  )
}))

drc_data_summary_table <- data.table(
  item = c(
    "Analysis rows",
    "Weighted sample size",
    "Weighted under-five mortality",
    "Ranking variable",
    "Candidate predictors"
  ),
  value = c(
    format(nrow(congo_model_dt), big.mark = ",", scientific = FALSE),
    format(round(sum(congo_model_dt$sample_weight)), big.mark = ",", scientific = FALSE),
    sprintf("%.2f%%", 100 * weighted_mean_safe(congo_model_dt$deadu5_num, congo_model_dt$sample_weight)),
    "DHS wealth index",
    paste(unname(raw_congo_var_labels), collapse = "; ")
  )
)





library(survey)
library(rineq)
linear_decomp_formula <- stats::as.formula(
  paste("deadu5_num ~", paste(raw_congo_predictors, collapse = " + "))
)
linear_decomp_design <- survey::svydesign(
  ids = ~1,
  weights = ~sample_weight,
  data = raw_congo_model_dt
)
linear_decomp_fit <- survey::svyglm(
  linear_decomp_formula,
  design = linear_decomp_design,
  family = stats::quasibinomial()
)



linear_x <- stats::model.matrix(linear_decomp_formula, data = raw_congo_model_dt)
linear_beta <- stats::coef(linear_decomp_fit)
linear_terms <- setdiff(colnames(linear_x), "(Intercept)")

linear_term_info <- data.table(
  term = linear_terms,
  variable = canonicalize_regression_terms(linear_terms),
  linear_coefficient = unname(linear_beta[linear_terms]),
  mean_of_variable = vapply(
    linear_terms,
    function(term) weighted_mean_safe(linear_x[, term], raw_congo_model_dt$sample_weight),
    numeric(1L)
  )
)

if (requireNamespace("broom", quietly = TRUE)) {
  regression_model_summary <- as.data.table(
    as.data.frame(broom::tidy(linear_decomp_fit, conf.int = TRUE))
  )
  regression_model_summary[, variable := fifelse(
    term == "(Intercept)",
    term,
    canonicalize_regression_terms(term)
  )]
  regression_model_summary[, term_label := fifelse(
    term == "(Intercept)",
    "Intercept",
    format_variable_label(variable)
  )]
  data.table::setcolorder(
    regression_model_summary,
    c("term", "term_label", "variable", setdiff(names(regression_model_summary), c("term", "term_label", "variable")))
  )
} else {
  regression_model_summary <- data.table(
    note = "Install broom to build the regression model summary."
  )
}

decomposition_to_table <- function(object, digits = 4L) {
  data.table(
    term = names(object$rel_contribution),
    `Contribution (%)` = round(object$rel_contribution, digits),
    `Contribution (Abs)` = round(object$ci_contribution, digits),
    Elasticity = round(c(0, object$elasticities), digits),
    `Concentration Index` = c(NA_real_, round(object$partial_cis, digits)),
    `lower 5%` = c(NA_real_, round(object$confints[1, ], digits)),
    `upper 5%` = c(NA_real_, round(object$confints[2, ], digits)),
    Corrected = c("", ifelse(object$corrected_coefficients, "yes", "no"))
  )
}

summarise_rineq_decomposition <- function(object, method_name) {
  out <- merge(
    linear_term_info,
    decomposition_to_table(object, digits = 6L)[
      term != "residual",
      .(
        term,
        elasticity = Elasticity,
        concentration_index_of_variable = `Concentration Index`,
        regression_contribution = `Contribution (Abs)`
      )
    ],
    by = "term",
    all.y = TRUE
  )
  out[, variable := fifelse(!is.na(variable) & nzchar(variable), variable, term)]
  out[, `:=`(
    criterion = method_name,
    variable_label = format_variable_label(variable)
  )]

  total_abs_contribution <- out[, sum(abs(regression_contribution), na.rm = TRUE)]
  out[, regression_pct_contribution := if (
    is.finite(total_abs_contribution) && total_abs_contribution > .Machine$double.eps
  ) {
    100 * abs(regression_contribution) / total_abs_contribution
  } else {
    NA_real_
  }]

  data.table::setcolorder(
    out,
    c(
      "variable",
      "variable_label",
      "criterion",
      "linear_coefficient",
      "mean_of_variable",
      "elasticity",
      "concentration_index_of_variable",
      "regression_contribution",
      "regression_pct_contribution"
    )
  )
  out[order(-abs(regression_pct_contribution))]
}

rineq_ci_decomp <- rineq::decomposition(
  outcome = raw_congo_model_dt$deadu5_num,
  betas = unname(linear_beta[linear_terms]),
  mm = linear_x[, linear_terms, drop = FALSE],
  ranker = raw_congo_model_dt$wealth,
  wt = raw_congo_model_dt$sample_weight,
  correction = TRUE,
  citype = "CI"
)
rineq_ci_summary <- decomposition_to_table(rineq_ci_decomp, digits = 4L)
classical_ci_table <- summarise_rineq_decomposition(rineq_ci_decomp, "CI")



rineq_cig_decomp <- rineq::decomposition(
  outcome = raw_congo_model_dt$deadu5_num,
  betas = unname(linear_beta[linear_terms]),
  mm = linear_x[, linear_terms, drop = FALSE],
  ranker = raw_congo_model_dt$wealth,
  wt = raw_congo_model_dt$sample_weight,
  correction = TRUE,
  citype = "CIg"
)
rineq_cig_summary <- decomposition_to_table(rineq_cig_decomp, digits = 4L)
classical_cig_table <- summarise_rineq_decomposition(rineq_cig_decomp, "CIg")



rineq_cic_decomp <- rineq::decomposition(
  outcome = raw_congo_model_dt$deadu5_num,
  betas = unname(linear_beta[linear_terms]),
  mm = linear_x[, linear_terms, drop = FALSE],
  ranker = raw_congo_model_dt$wealth,
  wt = raw_congo_model_dt$sample_weight,
  correction = TRUE,
  citype = "CIc"
)
rineq_cic_summary <- decomposition_to_table(rineq_cic_decomp, digits = 4L)
classical_cic_table <- summarise_rineq_decomposition(rineq_cic_decomp, "CIc")



rineq_ciw_decomp <- rineq::decomposition(
  outcome = raw_congo_model_dt$deadu5_num,
  betas = unname(linear_beta[linear_terms]),
  mm = linear_x[, linear_terms, drop = FALSE],
  ranker = raw_congo_model_dt$wealth,
  wt = raw_congo_model_dt$sample_weight,
  correction = TRUE,
  citype = "CIw"
)
rineq_ciw_summary <- decomposition_to_table(rineq_ciw_decomp, digits = 4L)
classical_ciw_table <- summarise_rineq_decomposition(rineq_ciw_decomp, "CIw")



rineq_decomp_by_type <- list(
  CI = rineq_ci_decomp,
  CIg = rineq_cig_decomp,
  CIc = rineq_cic_decomp,
  CIw = rineq_ciw_decomp
)

rineq_summary_by_type <- list(
  CI = rineq_ci_summary,
  CIg = rineq_cig_summary,
  CIc = rineq_cic_summary,
  CIw = rineq_ciw_summary
)

rineq_summary <- rbindlist(lapply(names(rineq_summary_by_type), function(type) {
  out <- copy(rineq_summary_by_type[[type]])
  out[, criterion := type]
  data.table::setcolorder(out, c("criterion", setdiff(names(out), "criterion")))
  out
}), fill = TRUE)

classical_decomposition_table <- rbindlist(
  list(classical_ci_table, classical_cig_table, classical_cic_table, classical_ciw_table),
  fill = TRUE
)

if (isTRUE(runtime_params$run_cv)) {
  tree_grid <- as.data.table(ineqTrees::ci_tree_control_grid(
    minsplit = 500L,
    minbucket = 250L,
    minprob = 0.01,
    maxdepth = 10L,
    min_gain = 0.00001,
    min_relative_gain = c(0.05, 0.10, 0.20)
  ))

  tree_tuning <-tune_ci_tree(
    formula = congo_ci_formula,
    data = congo_model_dt,
    rank_name = "wealth",
    outcome_name = "deadu5_num",
    weights = congo_model_dt$sample_weight,
    type = criterion_types,
    control_grid = tree_grid,
    v = 10L,
    strata = "deadu5_num",
    seed = 20260508,
    metric = c("validation_gain", "relative_validation_gain"),
    refit = FALSE,
    control =control_ci_tune(save_pred = TRUE)
  )

  tree_best_by_type <-ci_select_best(tree_tuning, metric = "relative_validation_gain")
  tree_selection_table <-ci_fit_summary_table(
    tree_tuning,
    selected = tree_best_by_type,
    metrics = c(
      "train_gain",
      "validation_gain",
      "train_relative_gain",
      "relative_validation_gain"
    ),
    include_percent = FALSE
  )
  tree_selection_table[
    tree_best_by_type,
    terminal_nodes := i.mean_terminal_nodes,
    on = c("type", "grid_id")
  ]
  tree_selection_table[, selection_basis := "10-fold cross-validation"]

  tree_models_by_type <- setNames(
    lapply(seq_len(nrow(tree_best_by_type)), function(i) {
      setting <- tree_best_by_type[i]
     ci_tree(
        formula = congo_ci_formula,
        data = congo_model_dt,
        rank_name = "wealth",
        outcome_name = "deadu5_num",
        weights = congo_model_dt$sample_weight,
        type = setting$type[1L],
        control =ci_control_from_row(setting)
      )
    }),
    tree_best_by_type$type
  )
} else {
  tree_tuning <- NULL
  tree_best_by_type <- NULL
  tree_models_by_type <- setNames(
    lapply(criterion_types, function(type) {
     ci_tree(
        formula = congo_ci_formula,
        data = congo_model_dt,
        rank_name = "wealth",
        outcome_name = "deadu5_num",
        weights = congo_model_dt$sample_weight,
        type = type,
        control = main_tree_control
      )
    }),
    criterion_types
  )

  tree_selection_table <- rbindlist(lapply(names(tree_models_by_type), function(type) {
    root_impurity <-ci_root_impurity(
      data = congo_model_dt,
      rank_name = "wealth",
      outcome_name = "deadu5_num",
      weights = congo_model_dt$sample_weight,
      type = type
    )
    training_gain <-ci_tree_validation_gain(
      fit = tree_models_by_type[[type]],
      new_data = congo_model_dt,
      rank_name = "wealth",
      outcome_name = "deadu5_num",
      weights = congo_model_dt$sample_weight,
      type = type,
      root_impurity = root_impurity
    )
    data.table(
      type = type,
      mean_root_objective = root_impurity,
      mean_train_gain = training_gain,
      mean_train_relative_gain = training_gain / abs(root_impurity),
      mean_validation_gain = NA_real_,
      mean_validation_relative_gain = NA_real_,
      terminal_nodes = length(partykit::nodeids(tree_models_by_type[[type]], terminal = TRUE)),
      selection_basis = "fixed-control training screen"
    )
  }))
}

selected_tree_type <- if (isTRUE(runtime_params$run_cv)) {
  tree_selection_table[order(-mean_validation_relative_gain)][1L, type]
} else {
  tree_selection_table[order(-mean_train_relative_gain)][1L, type]
}
selected_tree_fit <- tree_models_by_type[[selected_tree_type]]
selected_tree_summary <- as.data.table(ineqTrees::ci_tree_terminal_summary(selected_tree_fit))

all_tree_terminal_summaries <- rbindlist(lapply(names(tree_models_by_type), function(type) {
  out <- as.data.table(ineqTrees::ci_tree_terminal_summary(tree_models_by_type[[type]]))
  out[, criterion := type]
  data.table::setcolorder(out, "criterion")
  out
}), fill = TRUE)

tree_variable_importance <- collect_variable_importance(tree_models_by_type)

tree_best_controls <- if (!is.null(tree_best_by_type) && nrow(tree_best_by_type)) {
  tree_best_by_type[
    ,
    .(type, minsplit, minbucket, minprob, maxdepth, min_gain, min_relative_gain)
  ]
} else {
  data.table(
    type = criterion_types,
    minsplit = main_tree_control$minsplit,
    minbucket = main_tree_control$minbucket,
    minprob = main_tree_control$minprob,
    maxdepth = main_tree_control$maxdepth,
    min_gain = main_tree_control$min_gain,
    min_relative_gain = main_tree_control$min_relative_gain
  )
}

tree_complexity_metrics <- if (!is.null(tree_tuning)) {
  collect_complexity_metrics(tree_tuning)
} else {
  data.table()
}

tree_min_relative_gain_path <- if (nrow(tree_complexity_metrics)) {
  min_relative_gain_path(tree_complexity_metrics, selected = tree_best_by_type)
} else {
  data.table()
}

tree_report_metrics <- c(
  "type",
  "mean_root_objective",
  "mean_train_gain",
  "std_err_train_gain",
  "mean_train_relative_gain",
  "std_err_train_relative_gain",
  "mean_validation_gain",
  "std_err_validation_gain",
  "mean_validation_relative_gain",
  "std_err_validation_relative_gain"
)

tree_fit_wide <- if (!is.null(tree_tuning)) {
 ci_fit_summary_table(
    tree_tuning,
    selected = tree_best_by_type,
    metrics = c(
      "train_gain",
      "validation_gain",
      "train_relative_gain",
      "relative_validation_gain"
    ),
    include_percent = FALSE
  )[, ..tree_report_metrics]
} else {
  out <- copy(tree_selection_table)
  for (col in tree_report_metrics) {
    if (!col %in% names(out)) {
      out[, (col) := NA_real_]
    }
  }
  out[, ..tree_report_metrics]
}

tree_fit_table <- fit_summary_long(tree_fit_wide)





if (requireNamespace("rpart", quietly = TRUE)) {
  lecturer_rpart_models <- fit_ci_rpart_comparison(
    data = congo_model_dt,
    predictors = congo_predictors,
    outcome = "deadu5_num",
    rank = "wealth",
    weights = "sample_weight"
  )
  lecturer_rpart_summary <- summarise_ci_rpart_comparison(lecturer_rpart_models)
} else {
  lecturer_rpart_models <- list()
  lecturer_rpart_summary <- data.table(
    note = "Install rpart to fit lecturer rpart comparison trees."
  )
}





if (isTRUE(runtime_params$run_forest)) {
  if (isTRUE(runtime_params$run_forest_cv)) {
    forest_grid <- as.data.table(ineqTrees::ci_tree_control_grid(
      minsplit = 500L,
      minbucket = 250L,
      minprob = 0.01,
      maxdepth = 10L,
      min_gain = 0.00001,
      min_relative_gain = c(0.05, 0.10, 0.20),
      mtry = max(1L, 3, floor(sqrt(length(congo_predictors)))),
      ntree = 200L
    ))

    forest_tuning <- run_forest_tuning_batches(
      formula = congo_ci_formula,
      data = congo_model_dt,
      forest_grid = forest_grid,
      criterion_types = criterion_types,
      tuning_metrics = c("validation_gain", "relative_validation_gain"),
      tuning_selection_metric = "relative_validation_gain",
      rank_name = "wealth",
      outcome_name = "deadu5_num",
      weights = congo_model_dt$sample_weight,
      folds = 10L,
      workers = 4L,
      progress_steps = min(20L, nrow(forest_grid)),
      log_file = here("logs", "forest_tuning.log"),
      seed = 20260508,
      perturb = list(replace = FALSE, fraction = 0.7)
    )
    forest_best_by_type <-ci_select_best(
      forest_tuning,
      metric = "relative_validation_gain"
    )
    forest_fit_controls <- forest_best_by_type
  } else {
    forest_tuning <- NULL
    forest_best_by_type <- NULL
    forest_fit_controls <- data.table(
      type = criterion_types,
      minsplit = main_tree_control$minsplit,
      minbucket = main_tree_control$minbucket,
      minprob = main_tree_control$minprob,
      maxdepth = main_tree_control$maxdepth,
      min_gain = main_tree_control$min_gain,
      min_relative_gain = main_tree_control$min_relative_gain,
      mtry = max(1L, floor(sqrt(length(congo_predictors)))),
      ntree = 100L
    )
  }

  forest_models_by_type <- setNames(
    lapply(seq_len(nrow(forest_fit_controls)), function(i) {
      setting <- forest_fit_controls[i]
     ci_forest(
        formula = congo_ci_formula,
        data = congo_model_dt,
        rank_name = "wealth",
        outcome_name = "deadu5_num",
        weights = congo_model_dt$sample_weight,
        type = setting$type[1L],
        control =ci_control_from_row(setting),
        ntree = setting$ntree[1L],
        mtry = setting$mtry[1L],
        perturb = list(replace = FALSE, fraction = 0.7)
      )
    }),
    forest_fit_controls$type
  )
  selected_forest_fit <- forest_models_by_type[[selected_tree_type]]
  forest_variable_importance <- collect_variable_importance(forest_models_by_type)

  forest_complexity_metrics <- if (!is.null(forest_tuning)) {
    collect_complexity_metrics(forest_tuning)
  } else {
    data.table()
  }

  forest_min_relative_gain_path <- if (nrow(forest_complexity_metrics)) {
    min_relative_gain_path(forest_complexity_metrics, selected = forest_best_by_type)
  } else {
    data.table()
  }

  forest_report_metrics <- c(
    "type",
    "mean_root_objective",
    "mean_train_gain",
    "std_err_train_gain",
    "mean_train_relative_gain",
    "std_err_train_relative_gain",
    "mean_validation_gain",
    "std_err_validation_gain",
    "mean_validation_relative_gain",
    "std_err_validation_relative_gain"
  )

  if (!is.null(forest_tuning)) {
    forest_selection_table <-ci_fit_summary_table(
      forest_tuning,
      selected = forest_best_by_type,
      metrics = c(
        "train_gain",
        "validation_gain",
        "train_relative_gain",
        "relative_validation_gain"
      ),
      include_percent = FALSE
    )
    if (!all(c("ntree", "mtry") %in% names(forest_selection_table))) {
      forest_selection_table <- merge(
        forest_selection_table,
        forest_fit_controls[, .(type, ntree, mtry)],
        by = "type",
        all.x = TRUE,
        sort = FALSE
      )
    }
    forest_selection_table[
      forest_fit_controls,
      `:=`(
        mean_terminal_nodes = i.mean_terminal_nodes,
        mean_max_depth = i.maxdepth
      ),
      on = "type"
    ]
    forest_selection_table[, selection_basis := "10-fold cross-validation"]
  } else {
    forest_selection_table <- rbindlist(lapply(names(forest_models_by_type), function(type) {
      fit <- forest_models_by_type[[type]]
      root_impurity <-ci_root_impurity(
        data = congo_model_dt,
        rank_name = "wealth",
        outcome_name = "deadu5_num",
        weights = congo_model_dt$sample_weight,
        type = type
      )
      training_gain <-ci_forest_validation_gain(
        fit = fit,
        new_data = congo_model_dt,
        rank_name = "wealth",
        outcome_name = "deadu5_num",
        weights = congo_model_dt$sample_weight,
        type = type,
        root_impurity = root_impurity
      )
      forest_summary <- as.data.table(ineqTrees::ci_forest_summary(fit))
      forest_summary[
        ,
        .(
          type,
          ntree,
          mtry,
          mean_root_objective = root_impurity,
          mean_train_gain = training_gain,
          mean_train_relative_gain = training_gain / abs(root_impurity),
          mean_validation_gain = NA_real_,
          mean_validation_relative_gain = NA_real_,
          mean_terminal_nodes,
          mean_max_depth,
          selection_basis = "fixed-control training screen"
        )
      ]
    }), fill = TRUE)
  }

  forest_fit_wide <- copy(forest_selection_table)
  for (col in forest_report_metrics) {
    if (!col %in% names(forest_fit_wide)) {
      forest_fit_wide[, (col) := NA_real_]
    }
  }
  forest_fit_wide <- forest_fit_wide[, ..forest_report_metrics]
  forest_fit_table <- fit_summary_long(forest_fit_wide)
  forest_best_controls <- forest_fit_controls[
    ,
    .(type, ntree, mtry, minsplit, minbucket, minprob, maxdepth, min_gain, min_relative_gain)
  ]

  forest_surrogate_control <- main_tree_control
  forest_surrogate_control$minsplit <- 1000L
  forest_surrogate_control$minbucket <- 500L
  forest_surrogate_control$maxdepth <- 10L
  forest_surrogate_control$min_relative_gain <- 0.20
  forest_surrogates_by_type <- lapply(forest_models_by_type, function(fit) {
   ci_forest_surrogate(
      forest_fit = fit,
      type = fit$type,
      control = forest_surrogate_control
    )
  })
  selected_forest_surrogate <-ci_forest_surrogate(
    forest_fit = selected_forest_fit,
    type = selected_tree_type,
    control = forest_surrogate_control
  )
  selected_forest_surrogate_summary <- as.data.table(
   ci_tree_terminal_summary(selected_forest_surrogate$fit)
  )
} else {
  forest_models_by_type <- list()
  forest_tuning <- NULL
  forest_best_by_type <- NULL
  forest_best_controls <- data.table()
  forest_complexity_metrics <- data.table()
  forest_min_relative_gain_path <- data.table()
  forest_fit_wide <- data.table()
  forest_fit_table <- data.table()
  selected_forest_fit <- NULL
  forest_surrogates_by_type <- list()
  selected_forest_surrogate <- NULL
  selected_forest_surrogate_summary <- data.table()
  forest_selection_table <- data.table(
    type = character(),
    ntree = integer(),
    mtry = integer(),
    mean_root_objective = numeric(),
    mean_train_gain = numeric(),
    mean_train_relative_gain = numeric(),
    mean_validation_gain = numeric(),
    mean_validation_relative_gain = numeric(),
    mean_terminal_nodes = numeric(),
    mean_max_depth = numeric(),
    selection_basis = character()
  )
  forest_variable_importance <- data.table(
    criterion = character(),
    variable = character(),
    importance = numeric(),
    rank = integer()
  )
}

if (isTRUE(runtime_params$run_shap)) {
  required_predictive_packages <- c(
    "dials",
    "fastshap",
    "hardhat",
    "parsnip",
    "ranger",
    "recipes",
    "rsample",
    "tune",
    "yardstick",
    "workflows"
  )
  missing_predictive_packages <- required_predictive_packages[
    !vapply(required_predictive_packages, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing_predictive_packages)) {
    stop(
      "Install packages for the predictive workflow: ",
      paste(missing_predictive_packages, collapse = ", "),
      call. = FALSE
    )
  }
}

if (isTRUE(runtime_params$run_shap)) {
  set.seed(20260528)
  ranger_source_dt <- copy(congo_model_dt)
  ranger_source_dt[, source_row_id := .I]

  ranger_source_summary <- ranger_source_dt[
    ,
    .(weighted_n = sum(sample_weight), rows = .N),
    by = .(deadu5_num)
  ][
    ,
    `:=`(
      sample = "Before undersampling",
      weighted_percent = 100 * weighted_n / sum(weighted_n)
    )
  ]

  death_rows <- ranger_source_dt[deadu5_num == 1, source_row_id]
  nondeath_rows <- ranger_source_dt[deadu5_num == 0, source_row_id]
  death_weight <- ranger_source_dt[deadu5_num == 1, sum(sample_weight)]
  nondeath_weight <- ranger_source_dt[deadu5_num == 0, sum(sample_weight)]
  target_death_share <- runtime_params$ranger_death_target_share

  if (is.finite(target_death_share) &&
      target_death_share > 0 &&
      target_death_share < 1 &&
      length(death_rows) &&
      length(nondeath_rows) &&
      death_weight / (death_weight + nondeath_weight) < target_death_share) {
    target_nondeath_weight <- death_weight * (1 - target_death_share) / target_death_share
    nondeath_order <- sample(nondeath_rows)
    nondeath_cum_weight <- cumsum(ranger_source_dt[nondeath_order, sample_weight])
    keep_nondeath <- nondeath_order[nondeath_cum_weight <= target_nondeath_weight]
    next_idx <- length(keep_nondeath) + 1L
    if (next_idx <= length(nondeath_order)) {
      previous_weight <- if (length(keep_nondeath)) {
        nondeath_cum_weight[length(keep_nondeath)]
      } else {
        0
      }
      next_weight <- nondeath_cum_weight[next_idx]
      if (abs(next_weight - target_nondeath_weight) < abs(previous_weight - target_nondeath_weight)) {
        keep_nondeath <- c(keep_nondeath, nondeath_order[next_idx])
      }
    }
    ranger_dt <- ranger_source_dt[source_row_id %in% c(death_rows, keep_nondeath)]
  } else {
    ranger_dt <- copy(ranger_source_dt)
  }

  ranger_sampling_summary <- rbindlist(
    list(
      ranger_source_summary,
      ranger_dt[
        ,
        .(weighted_n = sum(sample_weight), rows = .N),
        by = .(deadu5_num)
      ][
        ,
        `:=`(
          sample = "After undersampling",
          weighted_percent = 100 * weighted_n / sum(weighted_n)
        )
      ]
    ),
    fill = TRUE
  )
  ranger_sampling_summary[
    ,
    outcome := fifelse(deadu5_num == 1, "Died before age 5", "Alive at age 5")
  ]
  setcolorder(ranger_sampling_summary, c("sample", "outcome"))
  ranger_dt[, source_row_id := NULL]
  ranger_dt[, deadu5_factor := factor(
    ifelse(deadu5_num == 1, "Died", "Survived"),
    levels = c("Died", "Survived")
  )]
  ranger_dt[, sample_weight_case := hardhat::importance_weights(sample_weight)]
  ranger_predictor_terms <- congo_predictors
  ranger_feature_df <- as.data.frame(ranger_dt[, ranger_predictor_terms, with = FALSE])
  names(ranger_feature_df) <- ranger_predictor_terms
  ranger_feature_lookup <- data.table(
    variable = ranger_predictor_terms,
    syntactic_variable = make.names(ranger_predictor_terms, unique = FALSE)
  )
  ranger_fit_df <- data.frame(
    deadu5_factor = ranger_dt$deadu5_factor,
    ranger_feature_df,
    sample_weight_case = ranger_dt$sample_weight_case,
    check.names = FALSE
  )
} else {
  ranger_dt <- data.table()
  ranger_sampling_summary <- data.table()
  ranger_predictor_terms <- character()
  ranger_feature_lookup <- data.table(
    variable = character(),
    syntactic_variable = character()
  )
  ranger_fit_df <- data.frame()
}

if (isTRUE(runtime_params$run_shap)) {
  ranger_recipe <- recipes::recipe(deadu5_factor ~ ., data = ranger_fit_df)
} else {
  ranger_recipe <- NULL
}

if (isTRUE(runtime_params$run_shap)) {
  ranger_model_spec <- parsnip::rand_forest(
    mtry = tune::tune(),
    min_n = tune::tune(),
    trees = 500L
  ) |>
    parsnip::set_mode("classification") |>
    parsnip::set_engine("ranger", probability = TRUE, importance = "impurity")
} else {
  ranger_model_spec <- NULL
}

if (isTRUE(runtime_params$run_shap)) {
  ranger_resamples <- rsample::vfold_cv(
    ranger_fit_df,
    v = 5,
    strata = deadu5_factor
  )
} else {
  ranger_resamples <- NULL
}

if (isTRUE(runtime_params$run_shap)) {
  ranger_workflow <- workflows::workflow() |>
    workflows::add_model(ranger_model_spec) |>
    workflows::add_recipe(ranger_recipe) |>
    workflows::add_case_weights(sample_weight_case)
} else {
  ranger_workflow <- NULL
}

if (isTRUE(runtime_params$run_shap)) {
  ranger_parameter_set <- hardhat::extract_parameter_set_dials(ranger_workflow) |>
    update(
      mtry = dials::mtry(range = c(1L, length(ranger_predictor_terms))),
      min_n = dials::min_n(range = c(10L, 150L))
    )
  ranger_grid <- dials::grid_regular(
    ranger_parameter_set,
    levels = c(mtry = 4, min_n = 3)
  )
} else {
  ranger_parameter_set <- NULL
  ranger_grid <- data.frame()
}

if (isTRUE(runtime_params$run_shap)) {
  ranger_metric_set <- yardstick::metric_set(
    yardstick::roc_auc,
    yardstick::accuracy,
    yardstick::sens
  )
  ranger_tune_results <- tune::tune_grid(
    ranger_workflow,
    resamples = ranger_resamples,
    grid = ranger_grid,
    metrics = ranger_metric_set,
    control = tune::control_grid(save_pred = TRUE)
  )
  ranger_tuning <- as.data.table(tune::collect_metrics(ranger_tune_results))
} else {
  ranger_tune_results <- NULL
  ranger_tuning <- data.table()
}



if (isTRUE(runtime_params$run_shap)) {
  selected_ranger_setting <- tune::select_best(ranger_tune_results, metric = "roc_auc")
} else {
  selected_ranger_setting <- data.table()
}



if (isTRUE(runtime_params$run_shap)) {
  final_ranger_workflow <- tune::finalize_workflow(ranger_workflow, selected_ranger_setting)
  predictive_forest_fit <- workflows::fit(final_ranger_workflow, data = ranger_fit_df)
} else {
  final_ranger_workflow <- NULL
  predictive_forest_fit <- NULL
}

predict_ranger_risk <- function(object, newdata) {
  predict(object, new_data = as.data.frame(newdata), type = "prob")$.pred_Died
}

if (isTRUE(runtime_params$run_shap)) {
  shap_rows <- sort(sample.int(nrow(ranger_fit_df), min(runtime_params$shap_n, nrow(ranger_fit_df))))
  shap_x_background <- ranger_fit_df[, ranger_predictor_terms, drop = FALSE]
  shap_x_eval <- ranger_fit_df[shap_rows, ranger_predictor_terms, drop = FALSE]
  shap_prediction <- as.numeric(predict_ranger_risk(predictive_forest_fit, shap_x_eval))
} else {
  shap_rows <- integer()
  shap_x_background <- data.frame()
  shap_x_eval <- data.frame()
  shap_prediction <- numeric()
}

if (isTRUE(runtime_params$run_shap)) {
  selected_forest_shap <- fastshap::explain(
    object = predictive_forest_fit,
    X = shap_x_background,
    pred_wrapper = predict_ranger_risk,
    newdata = shap_x_eval,
    nsim = runtime_params$shap_nsim,
    adjust = runtime_params$shap_nsim > 1L
  )
  colnames(selected_forest_shap) <- canonicalize_model_terms(
    colnames(selected_forest_shap),
    congo_predictors
  )

  shap_feature_values <- copy(as.data.table(shap_x_eval))
  shap_feature_values[, observation_id := seq_len(.N)]
  selected_forest_shap_dt <- as.data.table(selected_forest_shap)
  selected_forest_shap_dt[, observation_id := seq_len(.N)]

  shap_summary_long <- melt(
    selected_forest_shap_dt,
    id.vars = "observation_id",
    variable.name = "feature",
    value.name = "shap_value",
    variable.factor = FALSE
  )
  shap_value_long <- melt(
    shap_feature_values,
    id.vars = "observation_id",
    variable.name = "feature",
    value.name = "feature_value",
    variable.factor = FALSE
  )
  shap_summary_long <- merge(
    shap_summary_long,
    shap_value_long,
    by = c("observation_id", "feature"),
    all.x = TRUE,
    sort = FALSE
  )
  shap_summary_long[
    ,
    feature_value_numeric := suppressWarnings(as.numeric(as.character(feature_value)))
  ]
  shap_summary_long[
    is.na(feature_value_numeric),
    feature_value_numeric := as.numeric(factor(feature_value)) - 1
  ]
  shap_summary_long[
    ,
    feature_value_scaled := {
      range_value <- range(feature_value_numeric, na.rm = TRUE)
      if (all(is.finite(range_value)) && diff(range_value) > .Machine$double.eps) {
        (feature_value_numeric - range_value[1L]) / diff(range_value)
      } else {
        rep(0.5, .N)
      }
    },
    by = feature
  ]
  shap_summary_long[, feature_label := format_variable_label(feature)]
  shap_summary_long[
    ,
    mean_abs_shap := mean(abs(shap_value), na.rm = TRUE),
    by = feature
  ]
  data.table::setcolorder(
    shap_summary_long,
    c(
      "feature",
      "feature_label",
      "observation_id",
      "shap_value",
      "feature_value",
      "feature_value_numeric",
      "feature_value_scaled",
      "mean_abs_shap"
    )
  )
} else {
  selected_forest_shap <- NULL
  shap_summary_long <- data.table()
}

if (isTRUE(runtime_params$run_shap)) {
  shap_ci_decomp <-shap_conc_decomp(
    shap = selected_forest_shap,
    rank = ranger_dt$wealth[shap_rows],
    type = "CI",
    prediction = shap_prediction,
    weights = ranger_dt$sample_weight[shap_rows]
  )
} else {
  shap_ci_decomp <- NULL
}



if (isTRUE(runtime_params$run_shap)) {
  shap_cig_decomp <-shap_conc_decomp(
    shap = selected_forest_shap,
    rank = ranger_dt$wealth[shap_rows],
    type = "CIg",
    prediction = shap_prediction,
    weights = ranger_dt$sample_weight[shap_rows]
  )
} else {
  shap_cig_decomp <- NULL
}



if (isTRUE(runtime_params$run_shap)) {
  shap_cic_decomp <-shap_conc_decomp(
    shap = selected_forest_shap,
    rank = ranger_dt$wealth[shap_rows],
    type = "CIc",
    prediction = shap_prediction,
    weights = ranger_dt$sample_weight[shap_rows]
  )
} else {
  shap_cic_decomp <- NULL
}



if (isTRUE(runtime_params$run_shap)) {
  shap_l_decomp <-shap_conc_decomp(
    shap = selected_forest_shap,
    rank = ranger_dt$wealth[shap_rows],
    type = "L",
    prediction = shap_prediction,
    weights = ranger_dt$sample_weight[shap_rows]
  )
} else {
  shap_l_decomp <- NULL
}



if (isTRUE(runtime_params$run_shap)) {
  selected_shap_contributions <- rbindlist(list(
    as.data.table(shap_ci_decomp$contributions)[, criterion := "CI"],
    as.data.table(shap_cig_decomp$contributions)[, criterion := "CIg"],
    as.data.table(shap_cic_decomp$contributions)[, criterion := "CIc"],
    as.data.table(shap_l_decomp$contributions)[, criterion := "L"]
  ), fill = TRUE)
  selected_shap_contributions[, feature := canonicalize_model_terms(feature, congo_predictors)]
  data.table::setcolorder(selected_shap_contributions, c("criterion", setdiff(names(selected_shap_contributions), "criterion")))

  selected_shap_diagnostics <- rbindlist(list(
    as.data.table(shap_ci_decomp$diagnostics)[, criterion := "CI"],
    as.data.table(shap_cig_decomp$diagnostics)[, criterion := "CIg"],
    as.data.table(shap_cic_decomp$diagnostics)[, criterion := "CIc"],
    as.data.table(shap_l_decomp$diagnostics)[, criterion := "L"]
  ), fill = TRUE)
  data.table::setcolorder(selected_shap_diagnostics, c("criterion", setdiff(names(selected_shap_diagnostics), "criterion")))
} else {
  selected_shap_contributions <- data.table(
    criterion = character(),
    feature = character(),
    D_k_SHAP = numeric(),
    pct_contribution = numeric(),
    abs_contribution = numeric()
  )
  selected_shap_diagnostics <- data.table()
}



linear_rank <- classical_decomposition_table[
  criterion == "CI"
][
  order(-abs(regression_pct_contribution)),
  .(variable, method = "Linear decomposition", rank = seq_len(.N))
]
tree_rank <- rank_from_named_vector(selected_tree_fit$variable.importance, "Selected CI-tree")
forest_rank <- if (!is.null(selected_forest_fit)) {
  rank_from_named_vector(selected_forest_fit$variable.importance, "Selected CI-forest")
} else {
  data.table(variable = character(), method = "Selected CI-forest", rank = integer())
}
shap_rank <- if (nrow(selected_shap_contributions)) {
  selected_shap_contributions[
    criterion == "CI"
  ][
    order(-abs_contribution),
    .(variable = feature, method = "Selected forest SHAP", rank = seq_len(.N))
  ]
} else {
  data.table(variable = character(), method = "Selected forest SHAP", rank = integer())
}
importance_long <- rbindlist(list(linear_rank, tree_rank, forest_rank, shap_rank), fill = TRUE)

selected_variables <- classical_decomposition_table[
  criterion == "CI"
][
  order(-abs(regression_pct_contribution))
]$variable[seq_len(min(10L, classical_decomposition_table[criterion == "CI", .N]))]
selected_variables <- selected_variables[!is.na(selected_variables) & nzchar(selected_variables)]

method_comparison_table <- data.table(variable = selected_variables)
method_comparison_table[, variable_label := format_variable_label(variable)]
regression_wide <- dcast(
  classical_decomposition_table[
    criterion %in% c("CI", "CIg", "CIc"),
    .(variable, criterion, regression_pct_contribution)
  ],
  variable ~ criterion,
  value.var = "regression_pct_contribution"
)
for (criterion in c("CI", "CIg", "CIc")) {
  if (!criterion %in% names(regression_wide)) {
    regression_wide[, (criterion) := NA_real_]
  }
}
method_comparison_table <- merge(
  method_comparison_table,
  regression_wide[
    ,
    .(
      variable,
      regression_ci_pct_contribution = CI,
      regression_cig_pct_contribution = CIg,
      regression_cic_pct_contribution = CIc
    )
  ],
  by = "variable",
  all.x = TRUE
)
if (nrow(selected_shap_contributions)) {
  shap_wide <- dcast(
    selected_shap_contributions[
      criterion %in% criterion_types,
      .(variable = feature, criterion, pct_contribution)
    ],
    variable ~ criterion,
    value.var = "pct_contribution"
  )
  for (criterion in criterion_types) {
    if (!criterion %in% names(shap_wide)) {
      shap_wide[, (criterion) := NA_real_]
    }
  }
  shap_wide <- shap_wide[
    ,
    .(
      variable,
      forest_shap_ci_pct_contribution = CI,
      forest_shap_cig_pct_contribution = CIg,
      forest_shap_cic_pct_contribution = CIc,
      forest_shap_l_pct_contribution = L
    )
  ]
} else {
  shap_wide <- data.table(
    variable = selected_variables,
    forest_shap_ci_pct_contribution = NA_real_,
    forest_shap_cig_pct_contribution = NA_real_,
    forest_shap_cic_pct_contribution = NA_real_,
    forest_shap_l_pct_contribution = NA_real_
  )
}
method_comparison_table <- merge(
  method_comparison_table,
  shap_wide,
  by = "variable",
  all.x = TRUE
)
method_comparison_table <- method_comparison_table[
  order(-abs(regression_ci_pct_contribution), -abs(forest_shap_ci_pct_contribution))
]

importance_wide <- dcast(
  importance_long,
  rank ~ method,
  value.var = "variable",
  fun.aggregate = function(x) paste(unique(x), collapse = ", ")
)
for (col in setdiff(names(importance_wide), "rank")) {
  importance_wide[, (col) := format_variable_label(get(col))]
}
importance_wide <- importance_wide[rank <= 10L]





sim_file <- here("data", "known_subgroup_sim_data.rda")
if (file.exists(sim_file)) {
  load(sim_file)
  known_subgroup_fit <- fit_known_subgroup_ci_tree(known_subgroup_sim_dt, type = selected_tree_type)
  known_subgroup_validation <- validate_known_subgroup_recovery(known_subgroup_fit, known_subgroup_sim_dt)
  known_subgroup_summary <- summarise_known_subgroup_tree(known_subgroup_fit, known_subgroup_sim_dt)
} else {
  known_subgroup_sim_dt <- make_known_subgroup_data(n = 5000L, seed = 20260528)
  known_subgroup_fit <- fit_known_subgroup_ci_tree(known_subgroup_sim_dt, type = selected_tree_type)
  known_subgroup_validation <- validate_known_subgroup_recovery(known_subgroup_fit, known_subgroup_sim_dt)
  known_subgroup_summary <- summarise_known_subgroup_tree(known_subgroup_fit, known_subgroup_sim_dt)
}





drc_results_objects <- list(
  generation_params = runtime_params,
  outcome_labels = outcome_labels,
  congo_model_dt = congo_model_dt,
  raw_congo_model_dt = raw_congo_model_dt,
  congo_var_labels = congo_var_labels,
  raw_congo_var_labels = raw_congo_var_labels,
  congo_predictors = congo_predictors,
  raw_congo_predictors = raw_congo_predictors,
  indicator_level_map = indicator_level_map,
  congo_ci_formula = congo_ci_formula,
  data_summary = drc_data_summary_table,
  outcome_distribution = outcome_plot_dt,
  candidate_predictors = candidate_predictor_table,
  root_ci = drc_root_ci_table,
  regression_model_summary = regression_model_summary,
  classical_decomposition = classical_decomposition_table,
  rineq_summary = rineq_summary,
  tree_tuning = tree_tuning,
  tree_best_by_type = tree_best_by_type,
  tree_best_controls = tree_best_controls,
  tree_selection = tree_selection_table,
  tree_fit_table = tree_fit_table,
  tree_complexity_metrics = tree_complexity_metrics,
  tree_min_relative_gain_path = tree_min_relative_gain_path,
  tree_models_by_type = tree_models_by_type,
  selected_tree_type = selected_tree_type,
  selected_tree_fit = selected_tree_fit,
  selected_tree_summary = selected_tree_summary,
  all_tree_terminal_summaries = all_tree_terminal_summaries,
  tree_variable_importance = tree_variable_importance,
  lecturer_rpart_models = lecturer_rpart_models,
  lecturer_rpart_summary = lecturer_rpart_summary,
  forest_tuning = forest_tuning,
  forest_best_by_type = forest_best_by_type,
  forest_best_controls = forest_best_controls,
  forest_fit_table = forest_fit_table,
  forest_complexity_metrics = forest_complexity_metrics,
  forest_min_relative_gain_path = forest_min_relative_gain_path,
  selected_forest_fit = selected_forest_fit,
  forest_surrogates_by_type = forest_surrogates_by_type,
  selected_forest_surrogate = selected_forest_surrogate,
  selected_forest_surrogate_summary = selected_forest_surrogate_summary,
  forest_selection = forest_selection_table,
  forest_variable_importance = forest_variable_importance,
  predictive_forest_tuning = ranger_tuning,
  predictive_forest_sampling = ranger_sampling_summary,
  predictive_forest_feature_lookup = ranger_feature_lookup,
  predictive_forest_fit = predictive_forest_fit,
  shap_contributions = selected_shap_contributions,
  shap_diagnostics = selected_shap_diagnostics,
  shap_summary_long = shap_summary_long,
  method_comparison = method_comparison_table,
  importance_wide = importance_wide,
  known_subgroup_sim_dt = known_subgroup_sim_dt,
  known_subgroup_fit = known_subgroup_fit,
  simulation_validation = known_subgroup_validation,
  simulation_summary = known_subgroup_summary
)

save(drc_results_objects, file = results_object_file)

data.table(
  saved_to = results_object_file,
  object_count = length(drc_results_objects),
  object_names = paste(names(drc_results_objects), collapse = ", ")
)
