params <-
list(run_cv = FALSE, run_forest = TRUE, run_shap = TRUE, render_appendix_tree_plots = TRUE, 
    shap_n = 300L, shap_nsim = 8L)

knitr::opts_chunk$set(
  fig.width = 10,
  fig.height = 6,
  dpi = 500,
  dev = "png",
  warning = FALSE,
  message = FALSE
)

library(tidyverse)
library(data.table)
library(grid)
library(knitr)
library(ineqTrees)

source(file.path("R", "report_helpers.R"))
source(file.path("R", "ci_known_subgroup_demo.R"))

if (!exists("params", inherits = FALSE)) {
  params <- list(
    run_cv = FALSE,
    run_forest = TRUE,
    run_shap = TRUE,
    render_appendix_tree_plots = TRUE,
    shap_n = 300,
    shap_nsim = 8
  )
}

if (requireNamespace("kableExtra", quietly = TRUE)) {
  kableExtra::use_latex_packages()
}

criterion_types <- c("CI", "CIg", "CIc", "L")
results_object_file <- file.path("data", "drc_report_results_objects.rda")

load(file.path("data", "congo_model_data.rda"))
congo_model_dt <- as.data.table(congo_model_dt)
congo_predictors <- names(congo_var_labels)
congo_ci_formula <- stats::as.formula(
  paste("cbind(wealth, deadu5_num) ~", paste(congo_predictors, collapse = " + "))
)

main_tree_control <- ineqTrees::ci_tree_control(
  minsplit = 500L,
  minbucket = 250L,
  minprob = 0.01,
  maxdepth = 10L,
  min_gain = 0.00001,
  min_relative_gain = 0.05
)

weighted_term_ci <- function(x, rank, weights, type = "CI") {
  ci_fun <- ineqTrees::ci_factory(type)
  ci_fun(cbind(rank, x), weights)
}

map_model_terms_to_predictors <- function(term_names, predictors) {
  vapply(term_names, function(term) {
    hit <- predictors[term == predictors | startsWith(term, predictors)]
    if (length(hit)) hit[which.max(nchar(hit))] else term
  }, character(1L))
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

format_variable_label <- function(x) {
  fifelse(x %in% names(congo_var_labels), unname(congo_var_labels[x]), x)
}

# Data-summary objects.
outcome_plot_dt <- congo_model_dt[
  ,
  .(weighted_n = sum(sample_weight)),
  by = .(deadu5_num)
][
  ,
  outcome_label := fifelse(deadu5_num == 1, "Died before age 5", "Alive or died after age 5")
][
  ,
  weighted_percent := 100 * weighted_n / sum(weighted_n)
]

candidate_predictor_table <- rbindlist(lapply(congo_predictors, function(var) {
  x <- congo_model_dt[[var]]
  level_dt <- congo_model_dt[
    ,
    .(weighted_n = sum(sample_weight)),
    by = var
  ][order(-weighted_n)]
  level_dt[, weighted_percent := 100 * weighted_n / sum(weighted_n)]

  data.table(
    variable = var,
    label = unname(congo_var_labels[var]),
    class = paste(class(x), collapse = "/"),
    n_levels = uniqueN(x),
    most_common_level = as.character(level_dt[[var]][1L]),
    most_common_percent = level_dt$weighted_percent[1L]
  )
}), fill = TRUE)

drc_root_ci_table <- rbindlist(lapply(criterion_types, function(type) {
  ci_fun <- ineqTrees::ci_factory(type)
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
    sprintf(
      "%.2f%%",
      100 * weighted_mean_safe(congo_model_dt$deadu5_num, congo_model_dt$sample_weight)
    ),
    "DHS wealth index",
    paste(unname(congo_var_labels), collapse = "; ")
  )
)

# Classical linear decomposition-style objects.
linear_decomp_formula <- stats::as.formula(
  paste("deadu5_num ~", paste(congo_predictors, collapse = " + "))
)
linear_decomp_fit <- stats::lm(
  linear_decomp_formula,
  data = congo_model_dt,
  weights = sample_weight
)

linear_x <- stats::model.matrix(linear_decomp_fit)
linear_beta <- stats::coef(linear_decomp_fit)
linear_terms <- setdiff(colnames(linear_x), "(Intercept)")
linear_mu_y <- weighted_mean_safe(congo_model_dt$deadu5_num, congo_model_dt$sample_weight)

linear_term_decomp <- rbindlist(lapply(linear_terms, function(term) {
  x <- linear_x[, term]
  beta <- linear_beta[[term]]
  mean_x <- weighted_mean_safe(x, congo_model_dt$sample_weight)
  ci_x <- weighted_term_ci(
    x = x,
    rank = congo_model_dt$wealth,
    weights = congo_model_dt$sample_weight,
    type = "CI"
  )
  wagstaff <- beta * mean_x * ci_x / linear_mu_y
  data.table(
    term = term,
    variable = map_model_terms_to_predictors(term, congo_predictors),
    linear_coefficient = beta,
    mean_of_variable = mean_x,
    concentration_index_of_variable = ci_x,
    wagstaff_contribution = wagstaff,
    erreygers_contribution = 4 * linear_mu_y * wagstaff
  )
}), fill = TRUE)

classical_decomposition_table <- linear_term_decomp[
  ,
  .(
    linear_coefficient = sum(linear_coefficient * mean_of_variable, na.rm = TRUE) /
      sum(mean_of_variable, na.rm = TRUE),
    mean_of_variable = sum(mean_of_variable, na.rm = TRUE),
    concentration_index_of_variable = stats::weighted.mean(
      concentration_index_of_variable,
      w = pmax(abs(wagstaff_contribution), .Machine$double.eps),
      na.rm = TRUE
    ),
    wagstaff_contribution = sum(wagstaff_contribution, na.rm = TRUE),
    erreygers_contribution = sum(erreygers_contribution, na.rm = TRUE)
  ),
  by = variable
][
  ,
  variable_label := format_variable_label(variable)
][
  order(-abs(wagstaff_contribution))
]

if (requireNamespace("rineq", quietly = TRUE)) {
  rineq_decomp <- rineq::contribution(linear_decomp_fit, congo_model_dt$wealth)
  rineq_summary <- as.data.table(as.data.frame(summary(rineq_decomp)), keep.rownames = "term")
} else {
  rineq_summary <- data.table(
    term = character(),
    note = "Install rineq to reproduce the package-based Wagstaff decomposition."
  )
}

# Tree model objects. If CV is enabled, CV selects the best control row by
# criterion; otherwise the fixed reporting control is used for all criteria.
if (isTRUE(params$run_cv)) {
  tree_grid <- as.data.table(ineqTrees::ci_tree_control_grid(
    minsplit = 500L,
    minbucket = 250L,
    minprob = 0.01,
    maxdepth = 10L,
    min_gain = 0.00001,
    min_relative_gain = c(0.05, 0.10, 0.20)
  ))

  tree_tuning <- ineqTrees::tune_ci_tree(
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
    control = ineqTrees::control_ci_tune(save_pred = TRUE)
  )

  tree_best_by_type <- ineqTrees::ci_select_best(
    tree_tuning,
    metric = "relative_validation_gain"
  )

  tree_selection_table <- ineqTrees::ci_fit_summary_table(
    tree_tuning,
    selected = tree_best_by_type,
    metrics = c("validation_gain", "relative_validation_gain"),
    include_percent = TRUE
  )
  tree_selection_table[, selection_basis := "10-fold cross-validation"]

  tree_models_by_type <- setNames(
    lapply(seq_len(nrow(tree_best_by_type)), function(i) {
      setting <- tree_best_by_type[i]
      ineqTrees::ci_tree(
        formula = congo_ci_formula,
        data = congo_model_dt,
        rank_name = "wealth",
        outcome_name = "deadu5_num",
        weights = congo_model_dt$sample_weight,
        type = setting$type[1L],
        control = ineqTrees::ci_control_from_row(setting)
      )
    }),
    tree_best_by_type$type
  )
} else {
  tree_tuning <- NULL
  tree_best_by_type <- NULL
  tree_models_by_type <- setNames(
    lapply(criterion_types, function(type) {
      ineqTrees::ci_tree(
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
    root_impurity <- ineqTrees::ci_root_impurity(
      data = congo_model_dt,
      rank_name = "wealth",
      outcome_name = "deadu5_num",
      weights = congo_model_dt$sample_weight,
      type = type
    )
    training_gain <- ineqTrees::ci_tree_validation_gain(
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

selected_tree_type <- if (isTRUE(params$run_cv)) {
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

# Forest and SHAP objects are generated, but only selected summaries are emitted
# in the main results section.
if (isTRUE(params$run_forest)) {
  forest_models_by_type <- setNames(
    lapply(criterion_types, function(type) {
      ineqTrees::ci_forest(
        formula = congo_ci_formula,
        data = congo_model_dt,
        rank_name = "wealth",
        outcome_name = "deadu5_num",
        weights = congo_model_dt$sample_weight,
        type = type,
        control = main_tree_control,
        ntree = 100L,
        mtry = max(1L, floor(sqrt(length(congo_predictors)))),
        perturb = list(replace = FALSE, fraction = 0.7)
      )
    }),
    criterion_types
  )
  selected_forest_fit <- forest_models_by_type[[selected_tree_type]]
  forest_variable_importance <- collect_variable_importance(forest_models_by_type)
} else {
  forest_models_by_type <- list()
  selected_forest_fit <- NULL
  forest_variable_importance <- data.table(
    criterion = character(),
    variable = character(),
    importance = numeric(),
    rank = integer()
  )
}

if (isTRUE(params$run_shap) && !is.null(selected_forest_fit)) {
  if (!requireNamespace("fastshap", quietly = TRUE)) {
    stop("Install fastshap to compute SHAP results.", call. = FALSE)
  }

  set.seed(20260528)
  shap_rows <- sort(sample.int(
    nrow(congo_model_dt),
    min(params$shap_n, nrow(congo_model_dt))
  ))
  shap_x_background <- as.data.frame(congo_model_dt[, ..congo_predictors])
  shap_x_eval <- as.data.frame(congo_model_dt[shap_rows, ..congo_predictors])
  shap_prediction <- as.numeric(predict_ci_forest_risk(selected_forest_fit, shap_x_eval))

  selected_forest_shap <- fastshap::explain(
    object = selected_forest_fit,
    X = shap_x_background,
    pred_wrapper = predict_ci_forest_risk,
    newdata = shap_x_eval,
    nsim = params$shap_nsim,
    adjust = TRUE
  )

  selected_shap_decomp <- ineqTrees::shap_conc_decomp(
    shap = selected_forest_shap,
    rank = congo_model_dt$wealth[shap_rows],
    type = selected_tree_type,
    prediction = shap_prediction,
    weights = congo_model_dt$sample_weight[shap_rows]
  )
  selected_shap_contributions <- as.data.table(selected_shap_decomp$contributions)
  selected_shap_diagnostics <- as.data.table(selected_shap_decomp$diagnostics)
} else {
  shap_rows <- integer()
  selected_forest_shap <- NULL
  selected_shap_decomp <- NULL
  selected_shap_contributions <- data.table(
    feature = character(),
    D_k_SHAP = numeric(),
    pct_contribution = numeric(),
    abs_contribution = numeric()
  )
  selected_shap_diagnostics <- data.table()
}

linear_rank <- classical_decomposition_table[
  order(-abs(wagstaff_contribution)),
  .(variable, method = "Linear decomposition", rank = seq_len(.N))
]
tree_rank <- rank_from_named_vector(
  selected_tree_fit$variable.importance,
  "Selected CI-tree"
)
forest_rank <- if (!is.null(selected_forest_fit)) {
  rank_from_named_vector(selected_forest_fit$variable.importance, "Selected CI-forest")
} else {
  data.table(variable = character(), method = "Selected CI-forest", rank = integer())
}
shap_rank <- if (nrow(selected_shap_contributions)) {
  selected_shap_contributions[
    order(-abs_contribution),
    .(variable = feature, method = "Selected forest SHAP", rank = seq_len(.N))
  ]
} else {
  data.table(variable = character(), method = "Selected forest SHAP", rank = integer())
}
importance_long <- rbindlist(list(linear_rank, tree_rank, forest_rank, shap_rank), fill = TRUE)

selected_variables <- unique(c(
  classical_decomposition_table[order(-abs(wagstaff_contribution)), head(variable, 6L)],
  tree_rank[rank <= 6L, variable],
  forest_rank[rank <= 6L, variable],
  shap_rank[rank <= 6L, variable]
))

method_comparison_table <- data.table(variable = selected_variables)
method_comparison_table[, variable_label := format_variable_label(variable)]
method_comparison_table <- merge(
  method_comparison_table,
  classical_decomposition_table[
    ,
    .(
      variable,
      classical_contribution = wagstaff_contribution,
      classical_rank = frank(-abs(wagstaff_contribution), ties.method = "first")
    )
  ],
  by = "variable",
  all.x = TRUE
)
method_comparison_table <- merge(
  method_comparison_table,
  tree_rank[, .(variable, ci_tree_rank = rank)],
  by = "variable",
  all.x = TRUE
)
method_comparison_table <- merge(
  method_comparison_table,
  forest_rank[, .(variable, ci_forest_rank = rank)],
  by = "variable",
  all.x = TRUE
)
method_comparison_table <- merge(
  method_comparison_table,
  selected_shap_contributions[
    ,
    .(
      variable = feature,
      shap_pct_contribution = pct_contribution,
      shap_rank = frank(-abs_contribution, ties.method = "first")
    )
  ],
  by = "variable",
  all.x = TRUE
)
method_comparison_table <- method_comparison_table[order(classical_rank, ci_tree_rank, shap_rank)]

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

# Simulation positive-control objects.
sim_file <- file.path("data", "known_subgroup_sim_data.rda")
if (file.exists(sim_file)) {
  load(sim_file)
  known_subgroup_fit <- fit_known_subgroup_ci_tree(
    data = known_subgroup_sim_dt,
    type = selected_tree_type
  )
  known_subgroup_validation <- validate_known_subgroup_recovery(
    fit = known_subgroup_fit,
    data = known_subgroup_sim_dt
  )
  known_subgroup_summary <- summarise_known_subgroup_tree(
    fit = known_subgroup_fit,
    data = known_subgroup_sim_dt
  )
} else {
  known_subgroup_sim_dt <- make_known_subgroup_data(n = 5000L, seed = 20260528)
  known_subgroup_fit <- fit_known_subgroup_ci_tree(
    data = known_subgroup_sim_dt,
    type = selected_tree_type
  )
  known_subgroup_validation <- validate_known_subgroup_recovery(
    fit = known_subgroup_fit,
    data = known_subgroup_sim_dt
  )
  known_subgroup_summary <- summarise_known_subgroup_tree(
    fit = known_subgroup_fit,
    data = known_subgroup_sim_dt
  )
}

drc_results_objects <- list(
  data_summary = drc_data_summary_table,
  outcome_distribution = outcome_plot_dt,
  candidate_predictors = candidate_predictor_table,
  root_ci = drc_root_ci_table,
  classical_decomposition = classical_decomposition_table,
  rineq_summary = rineq_summary,
  tree_selection = tree_selection_table,
  selected_tree_type = selected_tree_type,
  selected_tree_summary = selected_tree_summary,
  all_tree_terminal_summaries = all_tree_terminal_summaries,
  tree_variable_importance = tree_variable_importance,
  forest_variable_importance = forest_variable_importance,
  shap_contributions = selected_shap_contributions,
  shap_diagnostics = selected_shap_diagnostics,
  method_comparison = method_comparison_table,
  importance_wide = importance_wide,
  simulation_validation = known_subgroup_validation,
  simulation_summary = known_subgroup_summary
)

save(
  drc_results_objects,
  file = results_object_file
)

report_compact_table(
  drc_results_objects$data_summary,
  caption = "Analysis sample and target inequality quantity",
  column_widths = c("4cm", "10cm")
)

report_compact_table(
  drc_results_objects$root_ci,
  digits = 4,
  caption = "Root-node inequality in under-five mortality by concentration-index criterion"
)

report_compact_table(
  drc_results_objects$method_comparison[
    ,
    .(
      variable = variable_label,
      classical_contribution,
      classical_rank,
      ci_tree_rank,
      ci_forest_rank,
      shap_pct_contribution,
      shap_rank
    )
  ],
  digits = 3,
  caption = "Selected comparison of classical decomposition and tree-based importance",
  column_widths = c("4.4cm", rep(NA_character_, 6L))
)

report_compact_table(
  drc_results_objects$importance_wide,
  caption = "Ranked variables across the decomposition and tree-based methods",
  column_widths = c("1cm", rep("3.5cm", ncol(drc_results_objects$importance_wide) - 1L))
)

selected_tree_selection_row <- drc_results_objects$tree_selection[type == selected_tree_type]
keep_selection_cols <- intersect(
  c(
    "type",
    "selection_basis",
    "mean_root_objective",
    "mean_train_gain",
    "mean_train_relative_gain",
    "mean_validation_gain",
    "mean_validation_relative_gain",
    "terminal_nodes"
  ),
  names(selected_tree_selection_row)
)

report_compact_table(
  selected_tree_selection_row[, ..keep_selection_cols],
  digits = 4,
  caption = "Selection table for the reported concentration-index tree"
)

ci_report_tree_plot(
  fit = selected_tree_fit,
  data = congo_model_dt,
  outcome_name = "deadu5_num",
  outcome_label = "U5 death",
  ci_type = selected_tree_type,
  var_labels = congo_var_labels
)

## flowchart TD
##   A[Start with root node] --> B[Compute socioeconomic rank]
##   B --> C[Evaluate candidate predictor splits]
##   C --> D[Compute parent concentration-index impurity]
##   D --> E[Compute weighted child-node impurity]
##   E --> F[Calculate CI gain]
##   F --> G{Split admissible?}
##   G -->|Yes| H[Choose largest admissible gain]
##   H --> I[Create child nodes]
##   I --> B
##   G -->|No| J[Declare terminal node]
##   J --> K[Report subgroup rule, outcome level, and within-node CI]

report_compact_table(
  drc_results_objects$simulation_validation,
  caption = "Positive-control validation of recovery of the planted subgroup"
)

plot_known_subgroup_tree(
  fit = known_subgroup_fit,
  data = known_subgroup_sim_dt,
  type = selected_tree_type
)

report_compact_table(
  drc_results_objects$candidate_predictors,
  digits = 2,
  caption = "Candidate predictors and dominant categories in the DRC analysis data",
  column_widths = c("2.4cm", "3.5cm", rep(NA_character_, 4L))
)

report_compact_table(
  drc_results_objects$tree_selection,
  digits = 4,
  caption = "All tree criteria considered before selecting the reported tree"
)

ggplot(drc_results_objects$outcome_distribution, aes(x = outcome_label, y = weighted_percent)) +
  geom_col(fill = "#3268A8", width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", weighted_percent)), vjust = -0.35, size = 3.5) +
  labs(x = NULL, y = "Weighted percent") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggplot(congo_model_dt, aes(x = wealth, weight = sample_weight)) +
  geom_histogram(bins = 35, fill = "#00BFC4", color = "white", linewidth = 0.2) +
  labs(x = "Wealth index", y = "Weighted count") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

report_compact_table(
  drc_results_objects$classical_decomposition[
    ,
    .(
      variable = variable_label,
      linear_coefficient,
      mean_of_variable,
      concentration_index_of_variable,
      wagstaff_contribution,
      erreygers_contribution
    )
  ],
  digits = 4,
  caption = "Full classical linear decomposition-style table",
  column_widths = c("4.2cm", rep(NA_character_, 5L))
)

report_compact_table(
  drc_results_objects$rineq_summary,
  digits = 4,
  caption = "Package-based rineq decomposition output, when rineq is installed"
)

report_compact_table(
  drc_results_objects$all_tree_terminal_summaries[
    ,
    .(criterion, node, n, weight, depth, ci, outcome_percent, rule)
  ],
  digits = 3,
  caption = "Terminal-node summaries for all fitted concentration-index trees",
  column_widths = c(rep(NA_character_, 7L), "7cm")
)

for (criterion in names(tree_models_by_type)) {
  if (isTRUE(params$render_appendix_tree_plots)) {
    cat("\n\n#### Criterion: ", criterion, "\n\n", sep = "")
    ci_report_tree_plot(
      fit = tree_models_by_type[[criterion]],
      data = congo_model_dt,
      outcome_name = "deadu5_num",
      outcome_label = "U5 death",
      ci_type = criterion,
      var_labels = congo_var_labels
    )
  }
}

tree_importance_out <- copy(drc_results_objects$tree_variable_importance)
tree_importance_out[, variable_label := format_variable_label(variable)]

report_compact_table(
  tree_importance_out[, .(criterion, rank, variable = variable_label, importance)],
  digits = 4,
  caption = "Tree variable importance by concentration-index criterion"
)

forest_importance_out <- copy(drc_results_objects$forest_variable_importance)
forest_importance_out[, variable_label := format_variable_label(variable)]

report_compact_table(
  forest_importance_out[, .(criterion, rank, variable = variable_label, importance)],
  digits = 4,
  caption = "Forest variable importance by concentration-index criterion"
)

shap_out <- copy(drc_results_objects$shap_contributions)
shap_out[, feature_label := format_variable_label(feature)]

report_compact_table(
  shap_out[
    ,
    .(feature = feature_label, D_k_SHAP, pct_contribution, abs_contribution)
  ],
  digits = 4,
  caption = "SHAP concentration-index decomposition for the selected forest"
)

report_compact_table(
  drc_results_objects$shap_diagnostics,
  digits = 4,
  caption = "Diagnostics for the selected forest SHAP concentration-index decomposition"
)

report_compact_table(
  drc_results_objects$simulation_summary[
    ,
    .(node, n, rule, outcome_percent, planted_subgroup_share)
  ],
  digits = 3,
  caption = "Terminal nodes from the positive-control simulation tree",
  column_widths = c(rep(NA_character_, 2L), "7cm", rep(NA_character_, 2L))
)
