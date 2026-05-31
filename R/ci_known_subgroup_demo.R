# Small sourceable demo for checking whether a CI tree recovers a planted subgroup.

known_subgroup_labels <- function() {
  c(
    rural = "Residence",
    ed = "Mother education",
    reg = "Province",
    birth = "Birth group",
    male = "Child sex"
  )
}

make_known_subgroup_data <- function(n = 3000L, seed = 20260528) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  set.seed(seed)
  wealth <- stats::rbeta(n, shape1 = 1.3, shape2 = 1.3)
  poor <- 1 - wealth
  rural <- stats::runif(n) < stats::plogis(-0.8 + 2.5 * poor)
  no_education <- stats::runif(n) < stats::plogis(-1.2 + 2.1 * poor + 0.7 * rural)

  p_kasai <- 0.10 + 0.22 * poor
  p_kasai_central <- 0.08 + 0.20 * poor
  p_kinshasa <- 0.38 - 0.28 * poor
  p_kongo_central <- pmax(0.05, 1 - p_kasai - p_kasai_central - p_kinshasa)
  reg_probability <- cbind(p_kinshasa, p_kasai, p_kasai_central, p_kongo_central)
  reg_probability <- reg_probability / rowSums(reg_probability)
  reg_levels <- c("Kinshasa", "Kasai", "Kasai Central", "Kongo Central")
  reg <- vapply(seq_len(n), function(i) {
    sample(reg_levels, size = 1L, prob = reg_probability[i, ])
  }, character(1L))

  dt <- data.table::data.table(
    wealth = wealth,
    rural = ifelse(rural, "Rural", "Urban"),
    ed = ifelse(no_education, "b no education", "a education"),
    reg = reg,
    birth = sample(
      c("a first", "b 2-4 short", "c 2-4 long", "d 5+ short"),
      n,
      replace = TRUE,
      prob = c(0.25, 0.25, 0.35, 0.15)
    ),
    male = sample(c("Female", "Male"), n, replace = TRUE)
  )

  dt[, sample_weight := stats::runif(.N, min = 0.5, max = 2.0)]
  dt[, planted_subgroup := rural == "Rural" &
    ed == "b no education" &
    reg %in% c("Kasai", "Kasai Central")]

  poor_rank <- 1 - ((data.table::frankv(dt$wealth, ties.method = "average") - 0.5) / n)
  eta <- stats::qlogis(0.035) +
    0.25 * (dt$rural == "Rural") +
    0.25 * (dt$ed == "b no education") +
    3.20 * dt$planted_subgroup +
    2.60 * dt$planted_subgroup * poor_rank

  dt[, true_risk := pmin(pmax(stats::plogis(eta), 0.005), 0.65)]
  dt[, deadu5_num := stats::rbinom(.N, size = 1L, prob = true_risk)]

  factor_cols <- c("rural", "ed", "reg", "birth", "male")
  dt[, (factor_cols) := lapply(.SD, factor), .SDcols = factor_cols]
  dt[]
}

fit_known_subgroup_ci_tree <- function(
    data,
    type = "L",
    minsplit = 150L,
    minbucket = 60L,
    maxdepth = 4L,
    min_relative_gain = 0.001) {
  stopifnot(requireNamespace("ineqTrees", quietly = TRUE))

  predictors <- names(known_subgroup_labels())
  formula <- stats::as.formula(
    paste("cbind(wealth, deadu5_num) ~", paste(predictors, collapse = " + "))
  )

  ineqTrees::ci_tree(
    formula = formula,
    data = data,
    rank_name = "wealth",
    outcome_name = "deadu5_num",
    weights = data$sample_weight,
    type = type,
    control = ineqTrees::ci_tree_control(
      minsplit = minsplit,
      minbucket = minbucket,
      maxdepth = maxdepth,
      min_relative_gain = min_relative_gain
    )
  )
}

summarise_known_subgroup_tree <- function(fit, data) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stopifnot(requireNamespace("ineqTrees", quietly = TRUE))

  out <- data.table::as.data.table(ineqTrees::ci_tree_terminal_summary(fit))
  node_id <- stats::predict(fit, newdata = as.data.frame(data), type = "node")
  check <- data.table::as.data.table(data)
  check[, node := as.integer(node_id)]
  check[
    ,
    planted_subgroup_share := stats::weighted.mean(planted_subgroup, sample_weight),
    by = node
  ]
  planted <- unique(check[, .(node, planted_subgroup_share)])
  merge(out, planted, by = "node", all.x = TRUE)[order(node)]
}

validate_known_subgroup_recovery <- function(
    fit,
    data,
    planted_variables = c("rural", "ed", "reg"),
    planted_share_threshold = 0.80) {
  summary <- summarise_known_subgroup_tree(fit, data)
  used_variables <- names(fit$variable.importance)

  data.table::data.table(
    max_terminal_planted_share = max(summary$planted_subgroup_share, na.rm = TRUE),
    recovered_pure_planted_node = max(summary$planted_subgroup_share, na.rm = TRUE) >=
      planted_share_threshold,
    planted_variables_used = paste(intersect(planted_variables, used_variables), collapse = ", "),
    all_planted_variables_used = all(planted_variables %in% used_variables)
  )
}

plot_known_subgroup_tree <- function(fit, data, type = "L") {
  stopifnot(requireNamespace("ineqTrees", quietly = TRUE))
  stopifnot(requireNamespace("grid", quietly = TRUE))

  ci_fun <- ineqTrees::ci_factory(type)
  plot(
    fit,
    gp = grid::gpar(fontsize = 7),
    data = as.data.frame(data),
    var_labels = known_subgroup_labels(),
    terminal_stats = list(
      n = function(df) sum(df$sample_weight, na.rm = TRUE),
      u5_death = function(df) stats::weighted.mean(df$deadu5_num, df$sample_weight),
      planted = function(df) stats::weighted.mean(df$planted_subgroup, df$sample_weight),
      ci = function(df) ci_fun(cbind(df$wealth, df$deadu5_num), df$sample_weight)
    ),
    stat_labels = list(
      n = "weighted n",
      u5_death = "U5 death",
      planted = "planted share",
      ci = type
    ),
    stat_formatters = list(
      n = function(x) format(round(x), big.mark = ",", scientific = FALSE),
      u5_death = function(x) sprintf("%.1f%%", 100 * x),
      planted = function(x) sprintf("%.1f%%", 100 * x),
      ci = function(x) sprintf("%.3f", x)
    ),
    terminal_fill = "#d9d9d9",
    tp_args = list(width_lines = 11, height_lines = 5),
    tnex = 0.85
  )
}

run_known_subgroup_demo <- function(n = 3000L, seed = 20260528, type = "L") {
  data <- make_known_subgroup_data(n = n, seed = seed)
  fit <- fit_known_subgroup_ci_tree(data = data, type = type)
  list(
    data = data,
    fit = fit,
    summary = summarise_known_subgroup_tree(fit, data),
    validation = validate_known_subgroup_recovery(fit, data)
  )
}
