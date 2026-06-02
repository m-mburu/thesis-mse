# Reusable report helpers for concentration-index tree analyses.
library(tidyverse)
library(kableExtra)
library(data.table)
library(ineqTrees)

#' Return a fallback value when the input is NULL.
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Draw a reproducible analysis sample unless the full data should be kept.
sample_analysis_rows <- function(
  data,
  use_full = FALSE, n = 10000L,
  seed = 20260516
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  data <- as.data.table(data)
  if (isTRUE(use_full) || nrow(data) <= n) {
    return(data.table::copy(data))
  }

  set.seed(seed)
  rows <- sort(sample.int(nrow(data), n))
  data[rows]
}


#' Select a bounded tuning grid while preserving requested strata.
select_tuning_grid <- function(
  grid,
  max_rows,
  seed,
  stratify_cols = character()
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  grid <- as.data.table(grid)
  if (!is.finite(max_rows) || nrow(grid) <= max_rows) {
    return(data.table::copy(grid))
  }

  set.seed(seed)
  stratify_cols <- intersect(stratify_cols, names(grid))
  selected <- integer()

  for (col in stratify_cols) {
    for (value in unique(grid[[col]])) {
      idx <- which(grid[[col]] == value)
      available <- setdiff(idx, selected)
      if (!length(available)) {
        available <- idx
      }
      selected <- unique(c(selected, sample(available, 1L)))
      if (length(selected) >= max_rows) {
        return(grid[selected[seq_len(max_rows)]])
      }
    }
  }

  remaining <- setdiff(seq_len(nrow(grid)), selected)
  if (length(selected) < max_rows && length(remaining)) {
    selected <- c(selected, sample(remaining, max_rows - length(selected)))
  }

  grid[selected]
}

#' Suggest minsplit and minbucket grids scaled to the analysis sample size.
suggest_node_size_grid <- function(
  n,
  minsplit_prop_range = c(0.02, 0.30),
  grid_size = 3L,
  minsplit_bounds = c(100L, Inf),
  minbucket_bounds = c(50L, Inf),
  minbucket_fraction = 0.5,
  round_to = 1L
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(length(minsplit_prop_range) == 2L)
  stopifnot(length(minsplit_bounds) == 2L)
  stopifnot(length(minbucket_bounds) == 2L)

  grid_size <- as.integer(grid_size[1L])
  if (is.na(grid_size) || grid_size < 1L) {
    stop("`grid_size` must be a positive integer.", call. = FALSE)
  }

  prop_range <- sort(as.numeric(minsplit_prop_range))
  if (!all(is.finite(prop_range)) || any(prop_range <= 0) ||
    prop_range[2L] >= 1) {
    stop(
      "`minsplit_prop_range` must contain two proportions between 0 and 1.",
      call. = FALSE
    )
  }

  minsplit_prop <- if (grid_size == 1L) {
    mean(prop_range)
  } else {
    seq(prop_range[1L], prop_range[2L], length.out = grid_size)
  }

  round_count <- function(x) {
    as.integer(round(x / round_to) * round_to)
  }

  minsplit <- round_count(n * minsplit_prop)
  minsplit <- pmax(minsplit_bounds[1L], pmin(minsplit_bounds[2L], minsplit))

  minbucket_prop <- minsplit_prop * minbucket_fraction
  minbucket <- round_count(n * minbucket_prop)
  minbucket <- pmax(minbucket_bounds[1L], pmin(minbucket_bounds[2L], minbucket))

  pairs <- unique(data.table::data.table(
    minsplit_prop = minsplit_prop,
    minbucket_prop = minbucket_prop,
    minsplit = as.integer(minsplit),
    minbucket = as.integer(minbucket)
  ))

  list(
    minsplit = sort(unique(as.integer(minsplit))),
    minbucket = sort(unique(as.integer(minbucket))),
    pairs = pairs
  )
}

#' Compute a weighted mean after removing missing values and nonpositive weights.
weighted_mean_safe <- function(x, w) {
  keep <- stats::complete.cases(x, w) & w > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  stats::weighted.mean(x[keep], w[keep])
}

#' Create compact tree edge labels that show split values only.
tree_edge_panel_values_only <- function(
  obj,
  var_labels = NULL,
  digits = 3,
  fill = "white",
  justmin = 4,
  just = c("alternate", "increasing", "decreasing", "equal"),
  max_items = 3L,
  max_chars = 24L,
  plural_label = "levels",
  ...
) {
  meta <- obj$data

  compact_split_value <- function(label) {
    label <- gsub("\\s+", " ", trimws(label))
    if (!grepl(",", label, fixed = TRUE)) {
      return(paste(strwrap(label, width = max_chars), collapse = "\n"))
    }

    parts <- trimws(strsplit(label, ",", fixed = TRUE)[[1L]])
    if (length(parts) <= max_items && nchar(label) <= max_chars) {
      return(paste(strwrap(label, width = max_chars), collapse = "\n"))
    }

    shown <- paste(
      parts[seq_len(min(max_items, length(parts)))],
      collapse = ", "
    )
    paste0(
      length(parts), " ", plural_label, "\n",
      paste(strwrap(paste0(shown, ", ..."), width = max_chars), collapse = "\n")
    )
  }

  justfun <- function(i, split_labels) {
    myjust <- if (mean(nchar(split_labels)) > justmin) {
      match.arg(just, c("alternate", "increasing", "decreasing", "equal"))
    } else {
      "equal"
    }

    k <- length(split_labels)
    rval <- switch(myjust,
      equal = rep.int(0, k),
      alternate = rep(c(0.5, -0.5), length.out = k),
      increasing = seq(from = -k / 2, to = k / 2, by = 1),
      decreasing = seq(from = k / 2, to = -k / 2, by = -1)
    )

    grid::unit(0.5, "npc") + grid::unit(rval[i], "lines")
  }

  function(node, i) {
    split_info <- partykit::character_split(
      partykit::split_node(node),
      meta,
      digits = digits
    )
    split_labels <- vapply(
      split_info$levels,
      compact_split_value,
      character(1L)
    )

    y <- justfun(i, split_labels)
    label <- split_labels[[i]]
    label_lines <- strsplit(label, "\n", fixed = TRUE)[[1L]]
    widest_line <- label_lines[which.max(nchar(label_lines))]

    grid::grid.rect(
      y = y,
      gp = grid::gpar(fill = fill, col = NA),
      width = grid::unit(1.05, "strwidth", widest_line),
      height = grid::unit(max(1L, length(label_lines)), "lines")
    )
    grid::grid.text(label, y = y, just = "center")
  }
}
class(tree_edge_panel_values_only) <- "grapcon_generator"

#' Create compact tree inner-node labels sized to each node.
tree_inner_panel_compact_labeled <- function(
  obj,
  var_labels = NULL,
  id = FALSE,
  show_p = TRUE,
  gp = grid::gpar(),
  fill = "white",
  max_chars = 34L,
  min_chars = 8L,
  ...
) {
  meta <- obj$data

  pretty_var_name <- function(var_name) {
    if (!is.null(var_labels) && var_name %in% names(var_labels)) {
      return(unname(var_labels[[var_name]]))
    }
    gsub("_", " ", var_name)
  }

  extract_label <- function(node) {
    if (partykit::is.terminal(node)) {
      return(list(label = "", p = ""))
    }

    split_info <- partykit::character_split(partykit::split_node(node), meta)
    varlab <- pretty_var_name(split_info$name)

    plab <- ""
    if (show_p) {
      pvalue <- tryCatch({
        value <- partykit::info_node(node)$p.value
        if (length(value) != 1L) NA_real_ else as.numeric(value)
      }, error = function(e) NA_real_)
      plab <- if (is.na(pvalue)) {
        ""
      } else if (pvalue < 0.001) {
        "p < 0.001"
      } else {
        paste0("p = ", format(round(pvalue, 3), nsmall = 3))
      }
    }

    list(label = varlab, p = plab)
  }

  function(node) {
    lab <- extract_label(node)
    label_lines <- strwrap(lab$label, width = max_chars)
    if (!length(label_lines)) {
      label_lines <- ""
    }

    all_lines <- c(label_lines, if (nzchar(lab$p)) lab$p)
    widest_line <- all_lines[which.max(nchar(all_lines))]
    if (nchar(widest_line) < min_chars) {
      widest_line <- paste(rep("a", min_chars), collapse = "")
    }

    label_text <- paste(label_lines, collapse = "\n")
    line_count <- length(label_lines) + if (nzchar(lab$p)) 1L else 0L

    grid::pushViewport(grid::viewport(gp = gp))
    grid::pushViewport(
      grid::viewport(
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.5, "npc"),
        width = grid::unit(1.12, "strwidth", widest_line),
        height = grid::unit(line_count + 0.8, "lines")
      )
    )

    grid::grid.roundrect(
      r = grid::unit(0.12, "snpc"),
      gp = grid::gpar(fill = fill)
    )
    grid::grid.text(
      label_text,
      y = grid::unit(if (nzchar(lab$p)) 0.62 else 0.5, "npc")
    )
    if (nzchar(lab$p)) {
      grid::grid.text(
        lab$p,
        y = grid::unit(0.22, "npc"),
        gp = grid::gpar(cex = 0.85)
      )
    }

    grid::upViewport(2)
  }
}
class(tree_inner_panel_compact_labeled) <- "grapcon_generator"

#' Format a numeric gain value with fixed decimal places.
format_gain <- function(x, digits = 4L) {
  formatC(x, format = "f", digits = digits)
}

#' Format a proportion as a percentage with fixed decimal places.
format_percent <- function(x, digits = 2L) {
  formatC(100 * x, format = "f", digits = digits)
}

#' Load generated DRC report objects and expose the fields used by the report.
load_drc_report_results <- function(results_object_file) {
  if (!file.exists(results_object_file)) {
    stop(
      "Missing ", results_object_file,
      ". Render Generate_DRC_Results_Objects.qmd first.",
      call. = FALSE
    )
  }

  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(results_object_file, envir = load_env)
  if (!"drc_results_objects" %in% loaded_names) {
    stop(
      "Expected `drc_results_objects` in ", results_object_file, ".",
      call. = FALSE
    )
  }

  drc_results_objects <- load_env$drc_results_objects
  list(
    drc_results_objects = drc_results_objects,
    congo_model_dt = drc_results_objects$congo_model_dt,
    congo_var_labels = drc_results_objects$congo_var_labels,
    congo_predictors = drc_results_objects$congo_predictors,
    congo_ci_formula = drc_results_objects$congo_ci_formula,
    outcome_labels = drc_results_objects$outcome_labels %||% c(
      alive = "Alive at age 5",
      died = "Died before age 5"
    ),
    selected_tree_type = drc_results_objects$selected_tree_type,
    selected_tree_fit = drc_results_objects$selected_tree_fit,
    tree_models_by_type = drc_results_objects$tree_models_by_type,
    lecturer_rpart_models = drc_results_objects$lecturer_rpart_models,
    lecturer_rpart_summary =
      drc_results_objects$lecturer_rpart_summary %||% data.table(),
    known_subgroup_sim_dt = drc_results_objects$known_subgroup_sim_dt,
    known_subgroup_fit = drc_results_objects$known_subgroup_fit
  )
}

#' Build ranger analysis data with optional target-share undersampling.
build_ranger_analysis_data <- function(
  data,
  predictor_terms,
  undersample = NULL,
  outcome_name = "deadu5_num",
  weight_name = "sample_weight",
  seed = 20260528,
  outcome_labels = c(
    alive = "Alive at age 5",
    died = "Died before age 5"
  )
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  if (!requireNamespace("hardhat", quietly = TRUE)) {
    stop("Install hardhat to build ranger case weights.", call. = FALSE)
  }

  if (!is.null(undersample)) {
    undersample <- as.numeric(undersample[1L])
    if (!is.finite(undersample) || undersample <= 0 || undersample >= 1) {
      stop(
        "`undersample` must be NULL or a proportion between 0 and 1.",
        call. = FALSE
      )
    }
  }

  ranger_source_dt <- copy(as.data.table(data))
  ranger_source_dt[, source_row_id := .I]

  sampling_summary <- ranger_source_dt[
    ,
    .(weighted_n = sum(get(weight_name), na.rm = TRUE), rows = .N),
    by = .(deadu5_num = get(outcome_name))
  ][
    ,
    `:=`(
      sample = if (is.null(undersample)) {
        "Analysis data"
      } else {
        "Before undersampling"
      },
      weighted_percent = 100 * weighted_n / sum(weighted_n),
      outcome = fifelse(
        deadu5_num == 1,
        outcome_labels[["died"]],
        outcome_labels[["alive"]]
      )
    )
  ]

  ranger_dt <- copy(ranger_source_dt)
  if (!is.null(undersample)) {
    death_rows <- ranger_source_dt[get(outcome_name) == 1, source_row_id]
    nondeath_rows <- ranger_source_dt[get(outcome_name) == 0, source_row_id]
    death_weight <- ranger_source_dt[
      get(outcome_name) == 1,
      sum(get(weight_name), na.rm = TRUE)
    ]
    nondeath_weight <- ranger_source_dt[
      get(outcome_name) == 0,
      sum(get(weight_name), na.rm = TRUE)
    ]
    current_death_share <- death_weight / (death_weight + nondeath_weight)

    if (is.finite(current_death_share) &&
      current_death_share < undersample &&
      length(death_rows) &&
      length(nondeath_rows)) {
      target_nondeath_weight <- death_weight * (1 - undersample) / undersample
      set.seed(seed)
      nondeath_order <- sample(nondeath_rows)
      nondeath_cum_weight <- cumsum(
        ranger_source_dt[nondeath_order, get(weight_name)]
      )
      keep_nondeath <- nondeath_order[nondeath_cum_weight <= target_nondeath_weight]
      next_idx <- length(keep_nondeath) + 1L
      if (next_idx <= length(nondeath_order)) {
        previous_weight <- if (length(keep_nondeath)) {
          nondeath_cum_weight[length(keep_nondeath)]
        } else {
          0
        }
        next_weight <- nondeath_cum_weight[next_idx]
        if (abs(next_weight - target_nondeath_weight) <
          abs(previous_weight - target_nondeath_weight)) {
          keep_nondeath <- c(keep_nondeath, nondeath_order[next_idx])
        }
      }
      ranger_dt <- ranger_source_dt[
        source_row_id %in% c(death_rows, keep_nondeath)
      ]
    }

    sampling_summary <- rbindlist(
      list(
        sampling_summary,
        ranger_dt[
          ,
          .(weighted_n = sum(get(weight_name), na.rm = TRUE), rows = .N),
          by = .(deadu5_num = get(outcome_name))
        ][
          ,
          `:=`(
            sample = "After undersampling",
            weighted_percent = 100 * weighted_n / sum(weighted_n),
            outcome = fifelse(
              deadu5_num == 1,
              outcome_labels[["died"]],
              outcome_labels[["alive"]]
            )
          )
        ]
      ),
      fill = TRUE
    )
  }

  setcolorder(sampling_summary, c("sample", "outcome"))
  ranger_dt[, source_row_id := NULL]
  ranger_dt[, deadu5_factor := factor(
    fifelse(get(outcome_name) == 1, "Died", "Survived"),
    levels = c("Died", "Survived")
  )]
  ranger_dt[, sample_weight_case := hardhat::importance_weights(get(weight_name))]

  ranger_feature_df <- as.data.frame(ranger_dt[, predictor_terms, with = FALSE])
  names(ranger_feature_df) <- predictor_terms
  ranger_feature_lookup <- data.table(
    variable = predictor_terms,
    syntactic_variable = make.names(predictor_terms, unique = FALSE)
  )
  ranger_fit_df <- data.frame(
    deadu5_factor = ranger_dt$deadu5_factor,
    ranger_feature_df,
    sample_weight_case = ranger_dt$sample_weight_case,
    check.names = FALSE
  )

  list(
    ranger_dt = ranger_dt,
    ranger_sampling_summary = sampling_summary,
    ranger_predictor_terms = predictor_terms,
    ranger_feature_lookup = ranger_feature_lookup,
    ranger_fit_df = ranger_fit_df
  )
}

#' Format model variable names using the DRC report label lookup.
format_variable_label <- function(x, var_labels = congo_var_labels) {
  predictor_names <- names(var_labels)
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(value)) {
      return(NA_character_)
    }
    if (value %in% predictor_names) {
      return(unname(var_labels[value]))
    }

    hits <- predictor_names[
      value == predictor_names | startsWith(value, predictor_names)
    ]
    if (length(hits)) {
      parent <- hits[which.max(nchar(hits))]
      suffix <- trimws(sub(paste0("^", parent), "", value))
      suffix <- gsub("^[:._]+", "", suffix)
      if (nzchar(suffix)) {
        return(paste0(unname(var_labels[parent]), ": ", suffix))
      }
      return(unname(var_labels[parent]))
    }

    value
  }, character(1L))
}

#' Plot a top-n variable-importance panel for a single model.
make_importance_plot <- function(
  importance_dt,
  title,
  fill,
  top_n = 10L,
  var_labels = congo_var_labels
) {
  importance_dt <- as.data.table(importance_dt)
  if (!nrow(importance_dt) || !"importance" %in% names(importance_dt)) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "Not available") +
        labs(title = title) +
        theme_void(base_size = 11)
    )
  }

  importance_dt <- importance_dt[is.finite(importance) & importance != 0]
  if (!nrow(importance_dt)) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No non-zero importance") +
        labs(title = title) +
        theme_void(base_size = 11)
    )
  }

  importance_dt[, variable_label := format_variable_label(variable, var_labels)]
  importance_dt[
    ,
    variable_label := vapply(
      variable_label,
      function(value) paste(strwrap(value, width = 34), collapse = "\n"),
      character(1L)
    )
  ]
  setorder(importance_dt, -importance)
  importance_dt <- head(importance_dt, top_n)
  importance_dt[, variable_label := make.unique(variable_label)]
  importance_dt[, variable_label := factor(variable_label, levels = rev(variable_label))]

  ggplot(importance_dt, aes(x = importance, y = variable_label)) +
    geom_col(fill = fill, width = 0.68) +
    labs(title = title, x = "Importance", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

#' Extract variable importance from a fitted surrogate tree object.
extract_surrogate_importance <- function(surrogate_object) {
  if (is.null(surrogate_object) ||
    !length(surrogate_object) ||
    is.null(surrogate_object$fit$variable.importance)) {
    return(data.table(variable = character(), importance = numeric()))
  }

  importance <- surrogate_object$fit$variable.importance
  data.table(
    variable = names(importance),
    importance = as.numeric(importance)
  )[order(-importance)]
}

#' Plot a criterion-specific forest surrogate tree from the DRC report objects.
plot_forest_surrogate <- function(
  type,
  fontsize = 6.2,
  results_objects = drc_results_objects,
  var_labels = congo_var_labels
) {
  surrogate <- results_objects$forest_surrogates_by_type[[type]]
  if (is.null(surrogate) || !length(surrogate)) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0, y = 0,
          label = paste("No", type, "forest surrogate available")
        ) +
        theme_void(base_size = 11)
    )
  }

  ci_report_tree_plot(
    fit = surrogate$fit,
    data = surrogate$data,
    outcome_name = surrogate$prediction_name,
    outcome_label = "Forest predicted U5 death",
    ci_type = type,
    var_labels = var_labels,
    fontsize = fontsize
  )
}

#' Plot a lecturer-style rpart comparison tree with report labels.
plot_lecturer_rpart <- function(
  method_name,
  fontsize = 4.4,
  models = lecturer_rpart_models,
  data = congo_model_dt,
  var_labels = congo_var_labels
) {
  fit <- models[[method_name]]
  if (is.null(fit) || !length(fit)) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0, y = 0,
          label = paste("No", method_name, "rpart tree available")
        ) +
        theme_void(base_size = 11)
    )
  }

  ci_report_tree_plot(
    fit = fit,
    data = data,
    outcome_name = "deadu5_num",
    outcome_label = "U5 death",
    ci_type = "CI",
    var_labels = var_labels,
    fontsize = fontsize
  )
}

#' Plot variable-importance rankings across concentration-index criteria.
plot_importance_by_criterion <- function(
  importance_dt,
  title,
  fill = "#3268A8",
  var_labels = congo_var_labels
) {
  plot_dt <- as.data.table(importance_dt)
  if (!nrow(plot_dt)) {
    return(
      ggplot() +
        annotate(
          "text",
          x = 0, y = 0,
          label = "No variable importance available"
        ) +
        labs(title = title, x = NULL, y = NULL) +
        theme_void(base_size = 11)
    )
  }

  plot_dt[, variable_label := format_variable_label(variable, var_labels)]
  plot_dt[
    ,
    variable_label := vapply(
      variable_label,
      function(value) paste(strwrap(value, width = 32), collapse = "\n"),
      character(1L)
    )
  ]
  plot_dt[
    ,
    facet_label := factor(
      paste(criterion, variable_label, sep = " | "),
      levels = rev(paste(criterion, variable_label, sep = " | "))
    )
  ]

  ggplot(plot_dt, aes(x = facet_label, y = importance)) +
    geom_col(width = 0.7, fill = fill) +
    coord_flip() +
    facet_wrap(vars(criterion), scales = "free_y", ncol = 2) +
    scale_x_discrete(labels = function(x) sub("^.* \\| ", "", x)) +
    labs(x = NULL, y = "Total split gain", title = title) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text.y = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

#' Build the weighted baseline table for the DRC appendix.
weighted_baseline_table <- function(
  results_objects,
  outcome_labels,
  caption = "Weighted distribution of analysis variables by under-five mortality status"
) {
  if (!requireNamespace("table1", quietly = TRUE)) {
    stop(
      "Install table1 to render the weighted appendix baseline table.",
      call. = FALSE
    )
  }

  baseline_dt <- copy(as.data.table(results_objects$raw_congo_model_dt))
  baseline_labels <- results_objects$raw_congo_var_labels
  baseline_vars <- setdiff(results_objects$raw_congo_predictors, "reg")
  baseline_dt[, u5_mortality := factor(
    fifelse(
      deadu5_num == 1,
      outcome_labels[["died"]],
      outcome_labels[["alive"]]
    ),
    levels = c(outcome_labels[["alive"]], outcome_labels[["died"]])
  )]

  baseline_table_df <- as.data.frame(
    baseline_dt[, c(baseline_vars, "u5_mortality"), with = FALSE]
  )
  for (variable in baseline_vars) {
    table1::label(baseline_table_df[[variable]]) <-
      unname(baseline_labels[variable])
  }

  baseline_weights <- baseline_dt$sample_weight
  render_weighted_categorical <- function(x, ...) {
    row_index <- table1:::indices.indexed(x)
    if (is.null(row_index)) {
      row_index <- seq_along(x)
    }

    keep <- !is.na(x)
    x_keep <- x[keep]
    w_keep <- baseline_weights[row_index][keep]
    levels_keep <- if (is.factor(x)) {
      levels(x)
    } else {
      sort(unique(as.character(x_keep)))
    }

    out <- vapply(levels_keep, function(level) {
      weighted_n <- sum(w_keep[as.character(x_keep) == level], na.rm = TRUE)
      weighted_denominator <- sum(w_keep, na.rm = TRUE)
      weighted_percent <- if (weighted_denominator > 0) {
        100 * weighted_n / weighted_denominator
      } else {
        NA_real_
      }
      sprintf(
        "%s (%.1f%%)",
        format(round(weighted_n, 1), big.mark = ",", scientific = FALSE),
        weighted_percent
      )
    }, character(1L))
    names(out) <- levels_keep
    c(" " = "", out)
  }

  render_weighted_strata <- function(strata, ...) {
    out <- vapply(seq_along(strata), function(i) {
      row_index <- table1:::indices.indexed(strata[[i]])
      weighted_n <- sum(baseline_weights[row_index], na.rm = TRUE)
      sprintf(
        "%s\n(weighted N=%s)",
        names(strata)[i],
        format(round(weighted_n), big.mark = ",", scientific = FALSE)
      )
    }, character(1L))
    names(out) <- out
    out
  }

  baseline_formula <- as.formula(
    paste("~", paste(baseline_vars, collapse = " + "), "| u5_mortality")
  )
  baseline_table1 <- table1::table1(
    baseline_formula,
    data = table1::indexed(baseline_table_df),
    overall = "Overall",
    rowlabelhead = "Variable",
    render.categorical = render_weighted_categorical,
    render.strat = render_weighted_strata
  )
  baseline_table_out <- as.data.frame(baseline_table1)
  names(baseline_table_out) <- c(
    "Variable",
    outcome_labels[["alive"]],
    outcome_labels[["died"]],
    "Overall"
  )

  baseline_table_kable <- report_compact_table(
    baseline_table_out,
    caption = caption,
    align = c("l", "r", "r", "r"),
    column_widths = c("4.5cm", rep("3.1cm", 3L)),
    escape = TRUE,
    bold_header = FALSE
  )

  if (requireNamespace("kableExtra", quietly = TRUE)) {
    baseline_table_kable <- kableExtra::add_header_above(
      baseline_table_kable,
      c(" " = 1, "Under-five mortality" = 2, "Overall" = 1),
      bold = TRUE
    )
    baseline_table_kable <- kableExtra::row_spec(
      baseline_table_kable,
      0,
      bold = TRUE
    )
  }

  baseline_table_kable
}

#' Create a report-ready kable table with shared styling defaults.
report_kable <- function(
  x,
  caption = NULL,
  digits = 4,
  format = NULL,
  align = NULL,
  booktabs = TRUE,
  longtable = TRUE,
  font_size = 11,
  position = "center",
  scale_down = TRUE,
  hold_position = TRUE,
  column_widths = NULL,
  bold_header = TRUE,
  bootstrap_options = c("striped", "condensed"),
  full_width = FALSE,
  ...
) {
  out <- knitr::kable(
    x,
    caption = caption,
    digits = digits,
    format = format,
    align = align,
    booktabs = booktabs,
    longtable = longtable,
    ...
  )
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    latex_options <- c(
      if (isTRUE(hold_position)) "hold_position",
      if (isTRUE(scale_down)) "scale_down"
    )
    out <- kableExtra::kable_styling(
      out,
      bootstrap_options = bootstrap_options,
      latex_options = latex_options,
      full_width = full_width,
      font_size = font_size,
      position = position
    )
    if (isTRUE(bold_header)) {
      out <- kableExtra::row_spec(out, 0, bold = TRUE)
    }
    if (!is.null(column_widths)) {
      for (i in seq_along(column_widths)) {
        if (!is.na(column_widths[[i]]) && nzchar(column_widths[[i]])) {
          out <- kableExtra::column_spec(out, i, width = column_widths[[i]])
        }
      }
    }
  }
  out
}

#' Create a LaTeX report table with compact defaults.
report_latex_table <- function(
  x,
  caption = NULL,
  digits = 4,
  align = NULL,
  column_widths = NULL,
  longtable = TRUE,
  font_size = 9,
  ...
) {
  report_kable(
    x,
    caption = caption,
    digits = digits,
    format = "latex",
    align = align,
    longtable = longtable,
    column_widths = column_widths,
    font_size = font_size,
    ...
  )
}

#' Create a compact report table for HTML or document output.
report_compact_table <- function(
  x,
  caption = NULL,
  digits = 4,
  align = NULL,
  column_widths = NULL,
  ...
) {
  report_kable(
    x,
    caption = caption,
    digits = digits,
    align = align,
    column_widths = column_widths,
    bootstrap_options = c("striped", "condensed", "hover"),
    full_width = FALSE,
    ...
  )
}

#' Render a fitted-model summary table for reports.
fit_summary_kable <- function(x, caption, digits = 4) {
  report_kable(x, caption = caption, digits = digits)
}

#' Convert fit summary metrics into a long table for reporting.
fit_summary_long <- function(x) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  x <- as.data.table(x)
  metric_levels <- c(
    "Root impurity",
    "Training gain",
    "Training gain (% root)",
    "Validation gain",
    "Validation gain (% root)"
  )

  out <- rbindlist(list(
    x[, .(
      criterion = type, metric = "Root impurity",
      estimate = mean_root_objective, std_error = NA_real_
    )],
    x[, .(
      criterion = type, metric = "Training gain",
      estimate = mean_train_gain, std_error = std_err_train_gain
    )],
    x[, .(
      criterion = type, metric = "Training gain (% root)",
      estimate = 100 * mean_train_relative_gain,
      std_error = 100 * std_err_train_relative_gain
    )],
    x[, .(
      criterion = type, metric = "Validation gain",
      estimate = mean_validation_gain, std_error = std_err_validation_gain
    )],
    x[, .(
      criterion = type, metric = "Validation gain (% root)",
      estimate = 100 * mean_validation_relative_gain,
      std_error = 100 * std_err_validation_relative_gain
    )]
  ))

  out[, criterion := factor(criterion, levels = unique(x$type))]
  out[, metric := factor(metric, levels = metric_levels)]
  setorder(out, criterion, metric)
  out[, criterion := as.character(criterion)]
  out[, metric := as.character(metric)]
  out[]
}

#' Plot a CI tree with weighted terminal-node summary statistics.
ci_report_tree_plot <- function(
  fit,
  data,
  outcome_name,
  outcome_label = "Outcome",
  ci_type = "CI",
  var_labels = NULL,
  rank_name = "wealth",
  weight_name = "sample_weight",
  fontsize = 11
) {
  stopifnot(requireNamespace("grid", quietly = TRUE))
  stopifnot(requireNamespace("ineqTrees", quietly = TRUE))
  stopifnot(requireNamespace("partykit", quietly = TRUE))

  if (inherits(fit, "rpart")) {
    fit <- partykit::as.party(fit)
  }

  ci_fun <- ineqTrees::ci_factory(ci_type)
  plot_fun <- if (inherits(fit, "ci_tree")) {
    plot
  } else if (inherits(fit, "party")) {
    getFromNamespace("plot.ci_tree", "ineqTrees")
  } else {
    plot
  }

  plot_fun(
    fit,
    gp = grid::gpar(fontsize = fontsize),
    data = as.data.frame(data),
    var_labels = var_labels,
    terminal_stats = list(
      weighted_n = function(df) sum(df[[weight_name]], na.rm = TRUE),
      outcome = function(df) {
        weighted_mean_safe(
          df[[outcome_name]],
          df[[weight_name]]
        )
      },
      mean_rank = function(df) {
        weighted_mean_safe(
          df[[rank_name]],
          df[[weight_name]]
        )
      },
      ci = function(df) {
        ci_fun(cbind(
          df[[rank_name]],
          df[[outcome_name]]
        ), df[[weight_name]])
      }
    ),
    stat_labels = list(
      weighted_n = "weighted n",
      outcome = outcome_label,
      mean_rank = paste("mean", rank_name),
      ci = ci_type
    ),
    stat_formatters = list(
      weighted_n = function(x) {
        format(round(x),
          big.mark = ",", scientific = FALSE
        )
      },
      outcome = function(x) sprintf("%.2f%%", 100 * x),
      mean_rank = function(x) sprintf("%.2f", x),
      ci = function(x) sprintf("%.3f", x)
    ),
    terminal_fill = "#d9d9d9",
    edge_panel = tree_edge_panel_values_only,
    inner_panel = tree_inner_panel_compact_labeled,
    tp_args = list(width_lines = 11, height_lines = 5.2),
    tnex = 0.85
  )
}

demo_tree_plot <- ci_report_tree_plot

#' Collect tuning complexity and performance metrics into one table.
collect_complexity_metrics <- function(tuning) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("ineqTrees", quietly = TRUE))
  metric_cols <- c(
    "train_gain",
    "validation_gain",
    "train_relative_gain",
    "relative_validation_gain"
  )

  out <- ineqTrees::ci_fit_summary_table(
    tuning,
    selected = NULL,
    metrics = metric_cols,
    include_percent = FALSE
  )

  terminal_nodes <- unique(data.table::as.data.table(tuning$summary)[
    ,
    .(grid_id, type, mean_terminal_nodes)
  ])

  out <- merge(
    as.data.table(out),
    terminal_nodes,
    by = c("grid_id", "type"),
    all.x = TRUE,
    sort = FALSE
  )

  keep_cols <- intersect(
    c(
      "grid_id", "type", "ntree", "mtry", "minsplit", "minbucket",
      "minprob", "maxdepth", "min_gain", "min_relative_gain",
      "mean_terminal_nodes", "mean_train_gain", "std_err_train_gain",
      "mean_validation_gain", "std_err_validation_gain",
      "mean_train_relative_gain", "std_err_train_relative_gain",
      "mean_validation_relative_gain", "std_err_validation_relative_gain",
      "mean_root_objective", "std_err_root_objective"
    ),
    names(out)
  )

  out <- out[, ..keep_cols]
  setorderv(out, intersect(c("type", "grid_id"), names(out)))
  out[]
}

#' Summarize tuning paths across minimum relative split gain values.
min_relative_gain_path <- function(complexity_dt, selected = NULL) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  required_cols <- c(
    "type",
    "grid_id",
    "min_relative_gain",
    "mean_validation_relative_gain",
    "mean_terminal_nodes"
  )
  missing_cols <- setdiff(required_cols, names(complexity_dt))
  if (length(missing_cols)) {
    stop(
      "Missing columns in complexity metrics: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  out <- copy(data.table::as.data.table(complexity_dt))
  out[
    ,
    `:=`(
      validation_percent_gain = 100 * mean_validation_relative_gain,
      validation_percent_se = if (
        "std_err_validation_relative_gain" %in% names(out)
      ) {
        100 * std_err_validation_relative_gain
      } else {
        NA_real_
      }
    )
  ]
  out <- out[
    is.finite(min_relative_gain) &
      is.finite(validation_percent_gain) &
      is.finite(mean_terminal_nodes)
  ]
  if (!nrow(out)) {
    return(out)
  }

  out[, selected := FALSE]
  if (!is.null(selected)) {
    selected_keys <- unique(
      as.data.table(selected)[, .(type, grid_id)]
    )
    out[selected_keys, selected := TRUE, on = c("type", "grid_id")]
  } else {
    out[
      ,
      `:=`(
        selected = get("validation_percent_gain") ==
          max(get("validation_percent_gain"), na.rm = TRUE)
      ),
      by = type
    ]
  }

  out <- out[
    ,
    {
      n_settings <- sum(is.finite(validation_percent_gain))
      se_values <- validation_percent_se[is.finite(validation_percent_se)]
      list(
        n_settings = n_settings,
        mean_validation_gain = mean(mean_validation_gain, na.rm = TRUE),
        mean_root_objective = mean(mean_root_objective, na.rm = TRUE),
        validation_percent_gain = mean(validation_percent_gain, na.rm = TRUE),
        validation_percent_se = if (length(se_values) && n_settings > 0L) {
          sqrt(sum(se_values^2)) / n_settings
        } else {
          NA_real_
        },
        validation_percent_from_mean_root = {
          mean_gain <- mean(mean_validation_gain, na.rm = TRUE)
          mean_root <- mean(mean_root_objective, na.rm = TRUE)
          if (is.finite(mean_gain) &&
            is.finite(mean_root) &&
            abs(mean_root) > .Machine$double.eps) {
            100 * mean_gain / abs(mean_root)
          } else {
            NA_real_
          }
        },
        min_validation_percent_gain = min(
          validation_percent_gain,
          na.rm = TRUE
        ),
        max_validation_percent_gain = max(
          validation_percent_gain,
          na.rm = TRUE
        ),
        mean_terminal_nodes = mean(mean_terminal_nodes, na.rm = TRUE),
        terminal_nodes_se = if (sum(is.finite(mean_terminal_nodes)) > 1L) {
          stats::sd(mean_terminal_nodes, na.rm = TRUE) /
            sqrt(sum(is.finite(mean_terminal_nodes)))
        } else {
          NA_real_
        },
        selected = any(selected)
      )
    },
    by = .(type, min_relative_gain)
  ]

  setorder(out, type, min_relative_gain)
  out[]
}

#' Plot validation percent gain across minimum relative split gain values.
plot_validation_percent_gain_path <- function(path_dt, title) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  if (!nrow(path_dt)) {
    return(
      ggplot() +
        annotate("text",
          x = 0, y = 0,
          label = "No finite tuning metrics"
        ) +
        labs(title = title, x = NULL, y = NULL) +
        theme_void(base_size = 11)
    )
  }

  plot_dt <- copy(data.table::as.data.table(path_dt))
  plot_dt[!is.finite(validation_percent_se), validation_percent_se := NA_real_]
  selected_dt <- plot_dt[selected == TRUE]

  ggplot(
    plot_dt,
    aes(
      x = min_relative_gain,
      y = validation_percent_gain, group = type
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2, color = "gray55", linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.6,
      color = "#00BFC4", alpha = 0.65, na.rm = TRUE
    ) +
    geom_errorbar(
      aes(
        ymin = validation_percent_gain - validation_percent_se,
        ymax = validation_percent_gain + validation_percent_se
      ),
      width = 0,
      alpha = 0.35,
      color = "#00BFC4",
      na.rm = TRUE
    ) +
    geom_point(size = 2, color = "#00BFC4", na.rm = TRUE) +
    geom_point(
      data = selected_dt,
      aes(x = min_relative_gain, y = validation_percent_gain),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = 0.9,
      fill = "white",
      color = "black",
      na.rm = TRUE
    ) +
    facet_wrap(ggplot2::vars(type), scales = "free_y") +
    labs(
      x = "Minimum relative split gain",
      y = "Validation gain (% root)",
      title = title
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

#' Plot mean terminal-node counts across minimum relative split gain values.
plot_terminal_nodes_path <- function(path_dt, title) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  if (!nrow(path_dt)) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No finite tuning metrics") +
        labs(title = title, x = NULL, y = NULL) +
        theme_void(base_size = 11)
    )
  }

  path_dt <- as.data.table(path_dt)
  selected_dt <- path_dt[selected == TRUE]

  ggplot(
    path_dt,
    aes(x = min_relative_gain, y = mean_terminal_nodes, group = type)
  ) +
    geom_line(linewidth = 0.6, color = "#3268A8", alpha = 0.75, na.rm = TRUE) +
    geom_errorbar(
      aes(
        ymin = mean_terminal_nodes - terminal_nodes_se,
        ymax = mean_terminal_nodes + terminal_nodes_se
      ),
      width = 0,
      alpha = 0.35,
      color = "#3268A8",
      na.rm = TRUE
    ) +
    geom_point(size = 2, color = "#3268A8", na.rm = TRUE) +
    geom_point(
      data = selected_dt,
      aes(x = min_relative_gain, y = mean_terminal_nodes),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = 0.9,
      fill = "white",
      color = "black",
      na.rm = TRUE
    ) +
    facet_wrap(ggplot2::vars(type), scales = "free_y") +
    labs(
      x = "Minimum relative split gain",
      y = "Mean terminal nodes",
      title = title
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

#' Extract top positive variable-importance scores from fitted models.
collect_variable_importance <- function(models, top_n = 12L) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  empty_importance <- data.table(
    criterion = character(),
    variable = character(),
    importance = numeric(),
    rank = integer()
  )

  out <- rbindlist(lapply(names(models), function(criterion_name) {
    importance <- models[[criterion_name]]$variable.importance
    if (is.null(importance) || !length(importance)) {
      return(empty_importance)
    }

    out <- data.table(
      criterion = criterion_name,
      variable = names(importance),
      importance = as.numeric(importance)
    )
    out <- out[is.finite(importance) & importance > 0]
    out <- out[order(-importance)]
    out[, rank := seq_len(.N)]
    out[rank <= top_n]
  }), fill = TRUE)

  if (!nrow(out)) {
    return(empty_importance)
  }

  out[]
}

#' Plot variable-importance rankings by criterion.
plot_variable_importance <- function(importance_dt, title, var_labels = NULL) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  plot_dt <- copy(data.table::as.data.table(importance_dt))
  if (!nrow(plot_dt)) {
    return(
      ggplot() +
        annotate("text",
          x = 0, y = 0,
          label = "No variable importance available"
        ) +
        labs(title = title, x = NULL, y = NULL) +
        theme_void(base_size = 11)
    )
  }

  if (is.null(var_labels)) {
    plot_dt[, variable_pretty := variable]
  } else {
    plot_dt[
      ,
      variable_pretty := ifelse(
        variable %in% names(var_labels),
        unname(var_labels[variable]),
        variable
      )
    ]
  }
  plot_dt[
    ,
    variable_label := factor(
      paste(criterion, variable_pretty, sep = " | "),
      levels = rev(paste(criterion, variable_pretty, sep = " | "))
    )
  ]

  ggplot(plot_dt, aes(x = variable_label, y = importance)) +
    geom_col(width = 0.7, fill = "#3268A8") +
    coord_flip() +
    facet_wrap(ggplot2::vars(criterion), scales = "free_y", ncol = 2) +
    scale_x_discrete(labels = function(x) sub("^.* \\| ", "", x)) +
    labs(x = NULL, y = "Total split gain", title = title) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

#' Predict CI forest risk using weighted terminal-node outcomes.
predict_ci_forest_risk <- function(
  object, newdata,
  outcome_name = object$outcome_name %||% "deadu5_num"
) {
  stats::predict(
    object,
    newdata = as.data.frame(newdata),
    type = "response",
    OOB = FALSE,
    FUN = function(y, w) {
      stats::weighted.mean(y[, outcome_name], w, na.rm = TRUE)
    }
  )
}

#' Run CI forest tuning in logged batches and combine the results.
run_forest_tuning_batches <- function(
  formula,
  data,
  forest_grid,
  criterion_types = c("CI", "CIg", "CIc", "L"),
  tuning_metrics = c("validation_gain", "relative_validation_gain"),
  tuning_selection_metric = tuning_metrics[2L],
  rank_name = "wealth",
  outcome_name = "deadu5_num",
  weights = data$sample_weight,
  folds = 10L,
  workers = 4L,
  progress_steps = 20L,
  log_file = "logs/forest_tuning.log",
  seed = 20260508,
  perturb = list(replace = FALSE, fraction = 0.7)
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("future", quietly = TRUE))
  stopifnot(requireNamespace("ineqTrees", quietly = TRUE))

  forest_grid <- as.data.table(forest_grid)
  log_dir <- dirname(log_file)
  if (!identical(log_dir, ".") && !dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(log_file)) {
    file.remove(log_file)
  }

  log_msg <- function(fmt, ...) {
    line <- sprintf(
      "[%s] %s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      sprintf(fmt, ...)
    )
    message(line)
    cat(line, "\n", file = log_file, append = TRUE)
    flush.console()
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)

  total_grid <- nrow(forest_grid)
  progress_steps <- max(1L, min(as.integer(progress_steps), total_grid))
  batch_size <- max(1L, ceiling(total_grid / progress_steps))
  batches <- split(
    seq_len(total_grid),
    ceiling(seq_len(total_grid) / batch_size)
  )
  total_fold_fits <- total_grid * length(criterion_types) * folds

  log_msg(
    "Forest tuning: %d grid rows, %d criteria, %d folds (~%s fold-level fits).",
    total_grid,
    length(criterion_types),
    folds,
    format(total_fold_fits, big.mark = ",", scientific = FALSE)
  )

  shift_grid_id <- function(dt, global_rows) {
    dt <- as.data.table(dt)
    if (nrow(dt) && "grid_id" %in% names(dt)) {
      keep <- !is.na(dt$grid_id)
      batch_n <- length(global_rows)
      local_grid_id <- as.integer(dt$grid_id[keep])
      dt[keep, grid_id := global_rows[((local_grid_id - 1L) %% batch_n) + 1L]]
    }
    dt
  }

  batch_results <- vector("list", length(batches))

  for (batch_idx in seq_along(batches)) {
    batch_rows <- batches[[batch_idx]]
    batch_grid <- forest_grid[batch_rows]
    batch_start <- Sys.time()

    log_msg(
      "Forest tuning: starting batch %d/%d, grid rows %d-%d (%d rows).",
      batch_idx,
      length(batches),
      min(batch_rows),
      max(batch_rows),
      length(batch_rows)
    )

    batch_fit <- ineqTrees::tune_ci_forest(
      formula = formula,
      data = data,
      rank_name = rank_name,
      outcome_name = outcome_name,
      weights = weights,
      type = criterion_types,
      control_grid = batch_grid,
      v = folds,
      strata = outcome_name,
      seed = seed,
      metric = tuning_metrics,
      refit = FALSE,
      perturb = perturb,
      parallel_over = "tuning",
      future.seed = TRUE,
      control = ineqTrees::control_ci_tune(save_pred = TRUE)
    )

    for (name in c(
      "fold_results", "summary",
      "predictions", "extracts", "fits", "notes"
    )) {
      batch_fit[[name]] <- shift_grid_id(batch_fit[[name]], batch_rows)
    }
    batch_results[[batch_idx]] <- batch_fit

    batch_minutes <- as.numeric(difftime(Sys.time(),
      batch_start,
      units = "mins"
    ))
    done <- sum(lengths(batches[seq_len(batch_idx)]))
    log_msg(
      paste0(
        "Forest tuning: completed batch %d/%d in %.1f minutes. ",
        "%.1f%% complete (%d/%d grid rows)."
      ),
      batch_idx,
      length(batches),
      batch_minutes,
      100 * done / total_grid,
      done,
      total_grid
    )
  }

  combined_summary <- rbindlist(lapply(
    batch_results,
    `[[`, "summary"
  ), fill = TRUE)
  setorderv(combined_summary, c("metric", "type", "grid_id"), c(1L, 1L, 1L))

  out <- list(
    fold_results = rbindlist(lapply(
      batch_results,
      `[[`, "fold_results"
    ), fill = TRUE),
    summary = combined_summary,
    best_params = data.table(),
    metric = tuning_metrics,
    selection_metric = tuning_selection_metric,
    model = "forest",
    fold_id = batch_results[[1L]]$fold_id,
    resamples = batch_results[[1L]]$resamples,
    control_grid = as.data.table(forest_grid),
    predictions = rbindlist(lapply(
      batch_results,
      `[[`, "predictions"
    ), fill = TRUE),
    extracts = rbindlist(lapply(
      batch_results,
      `[[`, "extracts"
    ), fill = TRUE),
    fits = rbindlist(lapply(
      batch_results,
      `[[`, "fits"
    ), fill = TRUE),
    notes = rbindlist(lapply(
      batch_results,
      `[[`, "notes"
    ), fill = TRUE),
    validation_roots = batch_results[[1L]]$validation_roots,
    control = batch_results[[1L]]$control
  )
  class(out) <- c("ci_forest_tuning", "ci_tree_tuning", class(out))
  out$best_params <- ineqTrees::ci_select_best(out)
  log_msg("Forest tuning: finished all batches.")
  out
}
