# Reusable report helpers for concentration-index tree analyses.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

sample_analysis_rows <- function(data, use_full = FALSE, n = 10000L, seed = 20260516) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  data <- data.table::as.data.table(data)
  if (isTRUE(use_full) || nrow(data) <= n) {
    return(data.table::copy(data))
  }

  set.seed(seed)
  rows <- sort(sample.int(nrow(data), n))
  data[rows]
}

select_tuning_grid <- function(grid, max_rows, seed, stratify_cols = character()) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  grid <- data.table::as.data.table(grid)
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

suggest_node_size_grid <- function(
    n,
    minsplit_prop_range = c(0.02, 0.30),
    grid_size = 3L,
    minsplit_bounds = c(100L, Inf),
    minbucket_bounds = c(50L, Inf),
    minbucket_fraction = 0.5,
    round_to = 1L) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(length(minsplit_prop_range) == 2L)
  stopifnot(length(minsplit_bounds) == 2L)
  stopifnot(length(minbucket_bounds) == 2L)

  grid_size <- as.integer(grid_size[1L])
  if (is.na(grid_size) || grid_size < 1L) {
    stop("`grid_size` must be a positive integer.", call. = FALSE)
  }

  prop_range <- sort(as.numeric(minsplit_prop_range))
  if (any(!is.finite(prop_range)) || any(prop_range <= 0) || prop_range[2L] >= 1) {
    stop("`minsplit_prop_range` must contain two proportions between 0 and 1.", call. = FALSE)
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

weighted_mean_safe <- function(x, w) {
  keep <- stats::complete.cases(x, w) & w > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  stats::weighted.mean(x[keep], w[keep])
}

format_gain <- function(x, digits = 4L) {
  formatC(x, format = "f", digits = digits)
}

format_percent <- function(x, digits = 2L) {
  formatC(100 * x, format = "f", digits = digits)
}

report_kable <- function(
    x,
    caption = NULL,
    digits = 4,
    format = NULL,
    align = NULL,
    booktabs = TRUE,
    longtable = FALSE,
    font_size = 9,
    position = "center",
    scale_down = TRUE,
    hold_position = TRUE,
    column_widths = NULL,
    bold_header = TRUE,
    bootstrap_options = c("striped", "condensed"),
    full_width = FALSE,
    ...) {
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

report_latex_table <- function(
    x,
    caption = NULL,
    digits = 4,
    align = NULL,
    column_widths = NULL,
    longtable = TRUE,
    font_size = 9,
    ...) {
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

report_compact_table <- function(
    x,
    caption = NULL,
    digits = 4,
    align = NULL,
    column_widths = NULL,
    ...) {
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

fit_summary_kable <- function(x, caption, digits = 4) {
  report_kable(x, caption = caption, digits = digits)
}

fit_summary_long <- function(x) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  x <- data.table::as.data.table(x)
  metric_levels <- c(
    "Root impurity",
    "Training gain",
    "Training gain (% root)",
    "Validation gain",
    "Validation gain (% root)"
  )

  out <- data.table::rbindlist(list(
    x[, .(criterion = type, metric = "Root impurity",
          estimate = mean_root_objective, std_error = NA_real_)],
    x[, .(criterion = type, metric = "Training gain",
          estimate = mean_train_gain, std_error = std_err_train_gain)],
    x[, .(criterion = type, metric = "Training gain (% root)",
          estimate = 100 * mean_train_relative_gain,
          std_error = 100 * std_err_train_relative_gain)],
    x[, .(criterion = type, metric = "Validation gain",
          estimate = mean_validation_gain, std_error = std_err_validation_gain)],
    x[, .(criterion = type, metric = "Validation gain (% root)",
          estimate = 100 * mean_validation_relative_gain,
          std_error = 100 * std_err_validation_relative_gain)]
  ))

  out[, criterion := factor(criterion, levels = unique(x$type))]
  out[, metric := factor(metric, levels = metric_levels)]
  data.table::setorder(out, criterion, metric)
  out[, criterion := as.character(criterion)]
  out[, metric := as.character(metric)]
  out[]
}

ci_report_tree_plot <- function(
    fit,
    data,
    outcome_name,
    outcome_label = "Outcome",
    ci_type = "CI",
    var_labels = NULL,
    rank_name = "wealth",
    weight_name = "sample_weight",
    fontsize = 6.5) {
  stopifnot(requireNamespace("grid", quietly = TRUE))
  stopifnot(requireNamespace("ineqTrees", quietly = TRUE))

  ci_fun <- ineqTrees::ci_factory(ci_type)
  plot(
    fit,
    gp = grid::gpar(fontsize = fontsize),
    data = as.data.frame(data),
    var_labels = var_labels,
    terminal_stats = list(
      weighted_n = function(df) sum(df[[weight_name]], na.rm = TRUE),
      outcome = function(df) weighted_mean_safe(df[[outcome_name]], df[[weight_name]]),
      mean_rank = function(df) weighted_mean_safe(df[[rank_name]], df[[weight_name]]),
      ci = function(df) ci_fun(cbind(df[[rank_name]], df[[outcome_name]]), df[[weight_name]])
    ),
    stat_labels = list(
      weighted_n = "weighted n",
      outcome = outcome_label,
      mean_rank = paste("mean", rank_name),
      ci = ci_type
    ),
    stat_formatters = list(
      weighted_n = function(x) format(round(x), big.mark = ",", scientific = FALSE),
      outcome = function(x) sprintf("%.2f%%", 100 * x),
      mean_rank = function(x) sprintf("%.2f", x),
      ci = function(x) sprintf("%.3f", x)
    ),
    terminal_fill = "#d9d9d9",
    tp_args = list(width_lines = 11, height_lines = 5.2),
    tnex = 0.85
  )
}

demo_tree_plot <- ci_report_tree_plot

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
    data.table::as.data.table(out),
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
  data.table::setorderv(out, intersect(c("type", "grid_id"), names(out)))
  out[]
}

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
    stop("Missing columns in complexity metrics: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  out <- data.table::copy(data.table::as.data.table(complexity_dt))
  out[
    ,
    `:=`(
      validation_percent_gain = 100 * mean_validation_relative_gain,
      validation_percent_se = if ("std_err_validation_relative_gain" %in% names(out)) {
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
    selected_keys <- unique(data.table::as.data.table(selected)[, .(type, grid_id)])
    out[selected_keys, selected := TRUE, on = c("type", "grid_id")]
  } else {
    out[, selected := validation_percent_gain == max(validation_percent_gain, na.rm = TRUE),
        by = type]
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
        min_validation_percent_gain = min(validation_percent_gain, na.rm = TRUE),
        max_validation_percent_gain = max(validation_percent_gain, na.rm = TRUE),
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

  data.table::setorder(out, type, min_relative_gain)
  out[]
}

plot_validation_percent_gain_path <- function(path_dt, title) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  if (!nrow(path_dt)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No finite tuning metrics") +
        ggplot2::labs(title = title, x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = 12)
    )
  }

  plot_dt <- data.table::copy(data.table::as.data.table(path_dt))
  plot_dt[!is.finite(validation_percent_se), validation_percent_se := NA_real_]
  selected_dt <- plot_dt[selected == TRUE]

  ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = min_relative_gain, y = validation_percent_gain, group = type)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "gray55", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.6, color = "#00BFC4", alpha = 0.65, na.rm = TRUE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = validation_percent_gain - validation_percent_se,
        ymax = validation_percent_gain + validation_percent_se
      ),
      width = 0,
      alpha = 0.35,
      color = "#00BFC4",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2, color = "#00BFC4", na.rm = TRUE) +
    ggplot2::geom_point(
      data = selected_dt,
      ggplot2::aes(x = min_relative_gain, y = validation_percent_gain),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = 0.9,
      fill = "white",
      color = "black",
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(ggplot2::vars(type), scales = "free_y") +
    ggplot2::labs(
      x = "Minimum relative split gain",
      y = "Validation gain (% root)",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

plot_terminal_nodes_path <- function(path_dt, title) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  if (!nrow(path_dt)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No finite tuning metrics") +
        ggplot2::labs(title = title, x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = 12)
    )
  }

  path_dt <- data.table::as.data.table(path_dt)
  selected_dt <- path_dt[selected == TRUE]

  ggplot2::ggplot(
    path_dt,
    ggplot2::aes(x = min_relative_gain, y = mean_terminal_nodes, group = type)
  ) +
    ggplot2::geom_line(linewidth = 0.6, color = "#3268A8", alpha = 0.75, na.rm = TRUE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_terminal_nodes - terminal_nodes_se,
        ymax = mean_terminal_nodes + terminal_nodes_se
      ),
      width = 0,
      alpha = 0.35,
      color = "#3268A8",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2, color = "#3268A8", na.rm = TRUE) +
    ggplot2::geom_point(
      data = selected_dt,
      ggplot2::aes(x = min_relative_gain, y = mean_terminal_nodes),
      inherit.aes = FALSE,
      shape = 21,
      size = 3,
      stroke = 0.9,
      fill = "white",
      color = "black",
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(ggplot2::vars(type), scales = "free_y") +
    ggplot2::labs(
      x = "Minimum relative split gain",
      y = "Mean terminal nodes",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

collect_variable_importance <- function(models, top_n = 12L) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  empty_importance <- data.table::data.table(
    criterion = character(),
    variable = character(),
    importance = numeric(),
    rank = integer()
  )

  out <- data.table::rbindlist(lapply(names(models), function(criterion_name) {
    importance <- models[[criterion_name]]$variable.importance
    if (is.null(importance) || !length(importance)) {
      return(empty_importance)
    }

    out <- data.table::data.table(
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

plot_variable_importance <- function(importance_dt, title, var_labels = NULL) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  plot_dt <- data.table::copy(data.table::as.data.table(importance_dt))
  if (!nrow(plot_dt)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No variable importance available") +
        ggplot2::labs(title = title, x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = 12)
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

  ggplot2::ggplot(plot_dt, ggplot2::aes(x = variable_label, y = importance)) +
    ggplot2::geom_col(width = 0.7, fill = "#3268A8") +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(ggplot2::vars(criterion), scales = "free_y", ncol = 2) +
    ggplot2::scale_x_discrete(labels = function(x) sub("^.* \\| ", "", x)) +
    ggplot2::labs(x = NULL, y = "Total split gain", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

predict_ci_forest_risk <- function(object, newdata, outcome_name = object$outcome_name %||% "deadu5_num") {
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
    perturb = list(replace = FALSE, fraction = 0.7)) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("future", quietly = TRUE))
  stopifnot(requireNamespace("ineqTrees", quietly = TRUE))

  forest_grid <- data.table::as.data.table(forest_grid)
  log_dir <- dirname(log_file)
  if (!identical(log_dir, ".") && !dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(log_file)) {
    file.remove(log_file)
  }

  log_msg <- function(fmt, ...) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(fmt, ...))
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
  batches <- split(seq_len(total_grid), ceiling(seq_len(total_grid) / batch_size))
  total_fold_fits <- total_grid * length(criterion_types) * folds

  log_msg(
    "Forest tuning: %d grid rows, %d criteria, %d folds (~%s fold-level fits).",
    total_grid,
    length(criterion_types),
    folds,
    format(total_fold_fits, big.mark = ",", scientific = FALSE)
  )

  shift_grid_id <- function(dt, global_rows) {
    dt <- data.table::as.data.table(dt)
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

    for (name in c("fold_results", "summary", "predictions", "extracts", "fits", "notes")) {
      batch_fit[[name]] <- shift_grid_id(batch_fit[[name]], batch_rows)
    }
    batch_results[[batch_idx]] <- batch_fit

    batch_minutes <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))
    done <- sum(lengths(batches[seq_len(batch_idx)]))
    log_msg(
      "Forest tuning: completed batch %d/%d in %.1f minutes. %.1f%% complete (%d/%d grid rows).",
      batch_idx,
      length(batches),
      batch_minutes,
      100 * done / total_grid,
      done,
      total_grid
    )
  }

  combined_summary <- data.table::rbindlist(lapply(batch_results, `[[`, "summary"), fill = TRUE)
  data.table::setorderv(combined_summary, c("metric", "type", "grid_id"), c(1L, 1L, 1L))

  out <- list(
    fold_results = data.table::rbindlist(lapply(batch_results, `[[`, "fold_results"), fill = TRUE),
    summary = combined_summary,
    best_params = data.table::data.table(),
    metric = tuning_metrics,
    selection_metric = tuning_selection_metric,
    model = "forest",
    fold_id = batch_results[[1L]]$fold_id,
    resamples = batch_results[[1L]]$resamples,
    control_grid = data.table::as.data.table(forest_grid),
    predictions = data.table::rbindlist(lapply(batch_results, `[[`, "predictions"), fill = TRUE),
    extracts = data.table::rbindlist(lapply(batch_results, `[[`, "extracts"), fill = TRUE),
    fits = data.table::rbindlist(lapply(batch_results, `[[`, "fits"), fill = TRUE),
    notes = data.table::rbindlist(lapply(batch_results, `[[`, "notes"), fill = TRUE),
    validation_roots = batch_results[[1L]]$validation_roots,
    control = batch_results[[1L]]$control
  )
  class(out) <- c("ci_forest_tuning", "ci_tree_tuning", class(out))
  out$best_params <- ineqTrees::ci_select_best(out)
  log_msg("Forest tuning: finished all batches.")
  out
}
