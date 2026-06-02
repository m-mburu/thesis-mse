# Modular helpers for concentration-index tree and forest simulation studies.

ci_sim_predictors <- function() {
  c(
    "unskilled", "male", "birth", "agemoth", "rural",
    "ed", "ped", "mocc", "pocc", "reg"
  )
}

ci_sim_feature_labels <- function() {
  c(
    unskilled = "Low-skill household",
    male = "Child sex",
    birth = "Birth order and interval",
    agemoth = "Maternal age",
    rural = "Residence",
    ed = "Maternal education",
    ped = "Partner education",
    mocc = "Maternal occupation",
    pocc = "Partner occupation",
    reg = "Province"
  )
}

ci_sim_simple_predictors <- function() {
  c("mother_education", "residence", "child_sex", "province")
}

ci_sim_simple_feature_labels <- function() {
  c(
    mother_education = "Mother's education",
    residence = "Residence",
    child_sex = "Child sex",
    province = "Province"
  )
}

ci_sim_simple_province_labels <- function() {
  c("Kinshasa", "Kongo Central", "Kasai", "Nord-Kivu")
}

ci_sim_make_simple_data <- function(
    n = 1000L,
    scenario = c("null_random", "known_subgroup"),
    seed = 20260602,
    outcome_name = "health_outcome") {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  scenario <- match.arg(scenario)
  set.seed(seed)

  dt <- data.table::data.table(
    mother_education = sample(
      c("No education", "Primary", "Secondary or higher"),
      size = n,
      replace = TRUE,
      prob = c(0.30, 0.40, 0.30)
    ),
    residence = sample(
      c("Urban", "Rural"),
      size = n,
      replace = TRUE,
      prob = c(0.35, 0.65)
    ),
    child_sex = sample(c("Female", "Male"), size = n, replace = TRUE),
    province = sample(
      ci_sim_simple_province_labels(),
      size = n,
      replace = TRUE,
      prob = c(0.25, 0.25, 0.25, 0.25)
    ),
    wealth = stats::runif(n),
    sample_weight = 1
  )

  wealth_rank <- ci_sim_fractional_rank(dt$wealth)
  poor <- wealth_rank <= 0.40

  # health_outcome = 1 is an adverse health outcome in this toy example.
  probability <- switch(
    scenario,
    null_random = rep(0.10, n),
    known_subgroup = 0.05 +
      0.20 * (dt$residence == "Rural" & dt$province == "Kasai" & poor) +
      0.12 * (dt$mother_education == "No education" & poor) +
      0.08 * (dt$residence == "Rural" & dt$province == "Nord-Kivu")
  )
  probability <- pmin(pmax(probability, 0.01), 0.60)

  dt[, (outcome_name) := stats::rbinom(.N, size = 1L, prob = probability)]
  dt[, `:=`(
    true_probability = probability,
    scenario = scenario,
    wealth_group = factor(
      ifelse(poor, "Poorer 40%", "Richer 60%"),
      levels = c("Poorer 40%", "Richer 60%")
    ),
    mother_education = factor(
      mother_education,
      levels = c("Secondary or higher", "Primary", "No education")
    ),
    residence = factor(residence, levels = c("Urban", "Rural")),
    child_sex = factor(child_sex, levels = c("Female", "Male")),
    province = factor(province, levels = ci_sim_simple_province_labels())
  )]

  data.table::setcolorder(
    dt,
    c(
      ci_sim_simple_predictors(),
      "wealth",
      "wealth_group",
      outcome_name,
      "true_probability",
      "scenario",
      "sample_weight"
    )
  )
  dt
}

ci_sim_make_simple_examples <- function(
    n = 1000L,
    seed = 20260602,
    outcome_name = "health_outcome") {
  list(
    random_no_relationship = ci_sim_make_simple_data(
      n = n,
      scenario = "null_random",
      seed = seed,
      outcome_name = outcome_name
    ),
    subgroup_relationship = ci_sim_make_simple_data(
      n = n,
      scenario = "known_subgroup",
      seed = seed + 1L,
      outcome_name = outcome_name
    )
  )
}

ci_sim_province_labels <- function() {
  c(
    "1" = "Kinshasa", "2" = "Kwango", "3" = "Kwilu",
    "4" = "Mai-Ndombe", "5" = "Kongo Central", "6" = "Equateur",
    "7" = "Mongala", "8" = "Nord-Ubangi", "9" = "Sud-Ubangi",
    "10" = "Tshuapa", "11" = "Kasai", "12" = "Kasai Central",
    "13" = "Kasai Oriental", "14" = "Lomami", "15" = "Sankuru",
    "16" = "Haut-Katanga", "17" = "Haut-Lomami", "18" = "Lualaba",
    "19" = "Tanganyika", "20" = "Maniema", "21" = "Nord-Kivu",
    "22" = "Bas-Uele", "23" = "Haut-Uele", "24" = "Ituri",
    "25" = "Tshopo", "26" = "Sud-Kivu"
  )
}

ci_sim_prepare_congo_base <- function(
    drc_data,
    n = 5000L,
    use_full = FALSE,
    seed = 20260524) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  dt <- data.table::as.data.table(drc_data)
  province_labels <- ci_sim_province_labels()

  dt[, `:=`(
    b5_num = as.numeric(b5),
    b7_num = as.numeric(b7),
    v191_num = as.numeric(v191),
    v025_num = as.numeric(v025),
    v024_num = as.integer(v024),
    v133_num = as.numeric(v133),
    v012_num = as.numeric(v012),
    v701_num = as.numeric(v701),
    v717_num = as.numeric(v717),
    v705_num = as.numeric(v705),
    bord_num = as.numeric(bord),
    b11_num = as.numeric(b11),
    b4_num = as.numeric(b4),
    v005_num = as.numeric(v005)
  )]

  dt[, `:=`(
    sample_weight = v005_num / 1000000,
    wealth = v191_num / 100000,
    reg = province_labels[as.character(v024_num)],
    rural = data.table::fcase(v025_num == 1, FALSE, v025_num == 2, TRUE, default = NA),
    ed = data.table::fcase(
      v133_num == 0, "b no education",
      !is.na(v133_num), "a education",
      default = NA_character_
    ),
    ped = data.table::fcase(
      v701_num == 0, "b no education",
      v701_num %in% c(1, 2, 3), "a education",
      default = NA_character_
    ),
    mocc = data.table::fcase(
      v717_num %in% c(0, 6, 9), "c Household, unskilled manual, not working",
      v717_num %in% c(4, 5, 10), "d Agriculture",
      v717_num %in% c(1, 2, 3, 7, 8, 96), "a other",
      v717_num == 97, "c Household, unskilled manual, not working",
      v717_num >= 98, NA_character_,
      default = NA_character_
    ),
    pocc = data.table::fcase(
      v705_num %in% c(0, 6, 9), "c Household, unskilled manual, not working",
      v705_num %in% c(4, 5, 10), "d Agriculture",
      v705_num %in% c(1, 2, 3, 7, 8, 96), "a other",
      v705_num == 97, "c Household, unskilled manual, not working",
      v705_num >= 98, NA_character_,
      default = NA_character_
    ),
    agemoth = data.table::fcase(
      v012_num < 20, "less than 20",
      v012_num >= 20, "a20 or more",
      default = NA_character_
    ),
    male = data.table::fcase(b4_num == 1, TRUE, b4_num == 2, FALSE, default = NA),
    birth = data.table::fcase(
      bord_num == 1, "a first",
      bord_num %in% 2:4 & !is.na(b11_num) & b11_num < 24, "b 2-4 short",
      bord_num %in% 2:4 & !is.na(b11_num) & b11_num >= 24, "c 2-4 long",
      bord_num > 4 & !is.na(b11_num) & b11_num < 24, "d 5+ short",
      bord_num > 4 & !is.na(b11_num) & b11_num >= 24, "e 5+ long",
      default = NA_character_
    ),
    unskilled = data.table::fcase(
      v717_num %in% c(0, 6, 9) | v705_num %in% c(0, 6, 9), TRUE,
      !is.na(v717_num) | !is.na(v705_num), FALSE,
      default = NA
    )
  )]

  out <- stats::na.omit(dt[, c("wealth", ci_sim_predictors(), "sample_weight"), with = FALSE])
  out <- out[sample_weight > 0]

  out[, `:=`(
    unskilled = factor(unskilled, levels = c(FALSE, TRUE), labels = c("No", "Yes")),
    male = factor(male, levels = c(FALSE, TRUE), labels = c("Female", "Male")),
    birth = factor(
      birth,
      levels = c("a first", "b 2-4 short", "c 2-4 long", "d 5+ short", "e 5+ long")
    ),
    agemoth = factor(agemoth, levels = c("a20 or more", "less than 20")),
    rural = factor(rural, levels = c(FALSE, TRUE), labels = c("Urban", "Rural")),
    ed = factor(ed, levels = c("a education", "b no education")),
    ped = factor(ped, levels = c("a education", "b no education")),
    mocc = factor(
      mocc,
      levels = c("a other", "c Household, unskilled manual, not working", "d Agriculture")
    ),
    pocc = factor(
      pocc,
      levels = c("a other", "c Household, unskilled manual, not working", "d Agriculture")
    ),
    reg = factor(reg, levels = unname(province_labels))
  )]

  if (!isTRUE(use_full) && nrow(out) > n) {
    set.seed(seed)
    out <- out[sort(sample.int(nrow(out), n))]
  }

  data.table::copy(out)
}

ci_sim_fractional_rank <- function(x) {
  (data.table::frankv(x, ties.method = "average") - 0.5) / length(x)
}

ci_sim_add_outcome <- function(
    data,
    scenario = c(
      "null",
      "simple_subgroup",
      "rank_mechanism",
      "level_mechanism",
      "complex_interaction"
    ),
    seed = 20260524,
    outcome_name = "deadu5_num") {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  scenario <- match.arg(scenario)
  dt <- data.table::copy(data.table::as.data.table(data))
  wealth_rank <- ci_sim_fractional_rank(dt$wealth)
  wealth_z <- as.numeric(scale(dt$wealth))
  poor_rank <- 0.5 - wealth_rank
  poor_level <- -wealth_z

  high_risk_province <- dt$reg %in% c(
    "Kasai", "Kasai Central", "Kasai Oriental", "Lomami", "Sankuru"
  )
  rural <- dt$rural == "Rural"
  no_maternal_education <- dt$ed == "b no education"
  short_or_high_order_birth <- dt$birth %in% c("b 2-4 short", "d 5+ short", "e 5+ long")
  agriculture_or_unskilled_partner <- dt$pocc %in% c(
    "c Household, unskilled manual, not working", "d Agriculture"
  )

  eta <- switch(
    scenario,
    null = stats::qlogis(0.075),
    simple_subgroup = stats::qlogis(0.060) +
      0.25 * rural +
      0.15 * high_risk_province +
      1.80 * (rural & high_risk_province) * poor_rank,
    rank_mechanism = stats::qlogis(0.065) +
      1.45 * rural * poor_rank +
      1.20 * no_maternal_education * poor_rank,
    level_mechanism = stats::qlogis(0.065) +
      0.90 * rural * poor_level +
      0.75 * no_maternal_education * poor_level,
    complex_interaction = stats::qlogis(0.055) +
      0.25 * rural +
      0.20 * no_maternal_education +
      0.15 * short_or_high_order_birth +
      1.65 * (rural & high_risk_province) * poor_rank +
      1.20 * (no_maternal_education & short_or_high_order_birth) * poor_rank +
      0.85 * agriculture_or_unskilled_partner * poor_rank
  )

  probability <- pmin(pmax(stats::plogis(eta), 0.005), 0.500)
  set.seed(seed)
  dt[, (outcome_name) := stats::rbinom(.N, size = 1L, prob = probability)]
  dt[, `:=`(
    true_probability = probability,
    scenario = scenario
  )]

  dt
}

ci_sim_formula <- function(
    predictors = ci_sim_predictors(),
    rank_name = "wealth",
    outcome_name = "deadu5_num") {
  stats::as.formula(
    paste(
      sprintf("cbind(%s, %s)", rank_name, outcome_name),
      "~",
      paste(predictors, collapse = " + ")
    )
  )
}

ci_sim_predict_tree <- function(fit, data, outcome_name = "deadu5_num", weights = NULL) {
  ineqTrees::predict_ci_tree_terminal_mean(
    fit = fit,
    train_data = data,
    new_data = data,
    outcome_name = outcome_name,
    weights = weights
  )
}

ci_sim_predict_forest <- function(fit, data) {
  stats::predict(
    fit,
    newdata = as.data.frame(data),
    type = "response",
    OOB = FALSE,
    FUN = function(y, w) stats::weighted.mean(y[, fit$outcome_name], w, na.rm = TRUE)
  )
}

ci_sim_weighted_brier <- function(truth, estimate, weights) {
  stats::weighted.mean((truth - estimate)^2, weights, na.rm = TRUE)
}

ci_sim_model_summary <- function(
    fit,
    model,
    data,
    type,
    rank_name = "wealth",
    outcome_name = "deadu5_num",
    weights = data$sample_weight) {
  root_impurity <- ineqTrees::ci_root_impurity(
    data = data,
    rank_name = rank_name,
    outcome_name = outcome_name,
    weights = weights,
    type = type
  )

  if (identical(model, "tree")) {
    pred <- ci_sim_predict_tree(fit, data, outcome_name, weights)
    gain <- ineqTrees::ci_tree_validation_gain(
      fit = fit,
      new_data = data,
      rank_name = rank_name,
      outcome_name = outcome_name,
      weights = weights,
      type = type,
      root_impurity = root_impurity
    )
    terminal_nodes <- length(partykit::nodeids(fit, terminal = TRUE))
    mean_terminal_nodes <- terminal_nodes
  } else {
    pred <- ci_sim_predict_forest(fit, data)
    gain <- ineqTrees::ci_forest_validation_gain(
      fit = fit,
      new_data = data,
      rank_name = rank_name,
      outcome_name = outcome_name,
      weights = weights,
      type = type,
      root_impurity = root_impurity
    )
    forest_summary <- ineqTrees::ci_forest_summary(fit)
    terminal_nodes <- NA_integer_
    mean_terminal_nodes <- forest_summary$mean_terminal_nodes
  }

  prediction_ci <- ineqTrees::ci_factory(type)(
    cbind(rank = data[[rank_name]], outcome = pred),
    weights
  )

  data.table::data.table(
    scenario = unique(data$scenario)[1],
    model = model,
    type = type,
    n = nrow(data),
    weighted_outcome_mean = stats::weighted.mean(data[[outcome_name]], weights),
    root_impurity = as.numeric(root_impurity),
    gain = as.numeric(gain),
    relative_gain = if (is.finite(root_impurity) && abs(root_impurity) > 0) {
      as.numeric(gain / abs(root_impurity))
    } else {
      NA_real_
    },
    prediction_ci = as.numeric(prediction_ci),
    brier = ci_sim_weighted_brier(data[[outcome_name]], pred, weights),
    terminal_nodes = terminal_nodes,
    mean_terminal_nodes = mean_terminal_nodes
  )
}

ci_sim_fit_models <- function(
    data,
    types = c("CI", "CIg", "CIc", "L"),
    predictors = ci_sim_predictors(),
    rank_name = "wealth",
    outcome_name = "deadu5_num",
    weights = data$sample_weight,
    tree_control = ineqTrees::ci_tree_control(
      minsplit = 250L,
      minbucket = 125L,
      maxdepth = 4L,
      min_relative_gain = 0.05
    ),
    forest_control = tree_control,
    ntree = 50L,
    mtry = max(1L, floor(sqrt(length(predictors)))),
    fit_forest = TRUE,
    seed = 20260524) {
  formula <- ci_sim_formula(predictors, rank_name, outcome_name)
  fits <- list()
  summaries <- list()

  for (type in types) {
    set.seed(seed)
    tree_fit <- ineqTrees::ci_tree(
      formula = formula,
      data = data,
      rank_name = rank_name,
      outcome_name = outcome_name,
      weights = weights,
      type = type,
      control = tree_control
    )
    fits[[paste(type, "tree", sep = "_")]] <- tree_fit
    summaries[[paste(type, "tree", sep = "_")]] <- ci_sim_model_summary(
      fit = tree_fit,
      model = "tree",
      data = data,
      type = type,
      rank_name = rank_name,
      outcome_name = outcome_name,
      weights = weights
    )

    if (isTRUE(fit_forest)) {
      set.seed(seed)
      forest_fit <- ineqTrees::ci_forest(
        formula = formula,
        data = data,
        rank_name = rank_name,
        outcome_name = outcome_name,
        weights = weights,
        type = type,
        control = forest_control,
        ntree = ntree,
        mtry = mtry,
        perturb = list(replace = FALSE, fraction = 0.632)
      )
      fits[[paste(type, "forest", sep = "_")]] <- forest_fit
      summaries[[paste(type, "forest", sep = "_")]] <- ci_sim_model_summary(
        fit = forest_fit,
        model = "forest",
        data = data,
        type = type,
        rank_name = rank_name,
        outcome_name = outcome_name,
        weights = weights
      )
    }
  }

  list(
    fits = fits,
    summary = data.table::rbindlist(summaries, fill = TRUE)
  )
}

ci_sim_permutation_importance <- function(
    fit,
    model = c("tree", "forest"),
    data,
    predictors = ci_sim_predictors(),
    outcome_name = "deadu5_num",
    weights = data$sample_weight,
    repeats = 3L,
    seed = 20260524) {
  model <- match.arg(model)
  dt <- data.table::as.data.table(data)
  base_pred <- if (identical(model, "tree")) {
    ci_sim_predict_tree(fit, dt, outcome_name, weights)
  } else {
    ci_sim_predict_forest(fit, dt)
  }
  base_brier <- ci_sim_weighted_brier(dt[[outcome_name]], base_pred, weights)

  out <- vector("list", length(predictors) * repeats)
  idx <- 1L
  for (feature in predictors) {
    for (rep_id in seq_len(repeats)) {
      permuted_dt <- data.table::copy(dt)
      set.seed(seed + rep_id + match(feature, predictors) * 1000L)
      permuted_dt[[feature]] <- sample(permuted_dt[[feature]])
      pred <- if (identical(model, "tree")) {
        ci_sim_predict_tree(fit, permuted_dt, outcome_name, weights)
      } else {
        ci_sim_predict_forest(fit, permuted_dt)
      }
      out[[idx]] <- data.table::data.table(
        feature = feature,
        repeat_id = rep_id,
        base_brier = base_brier,
        permuted_brier = ci_sim_weighted_brier(dt[[outcome_name]], pred, weights),
        importance = ci_sim_weighted_brier(dt[[outcome_name]], pred, weights) - base_brier
      )
      idx <- idx + 1L
    }
  }

  data.table::rbindlist(out)[
    ,
    .(
      base_brier = mean(base_brier),
      permuted_brier = mean(permuted_brier),
      importance = mean(importance)
    ),
    by = feature
  ][order(-importance)]
}

ci_sim_run_scenarios <- function(
    base_data,
    scenarios = c(
      "null",
      "simple_subgroup",
      "rank_mechanism",
      "level_mechanism",
      "complex_interaction"
    ),
    types = c("CI", "CIg", "CIc", "L"),
    ntree = 50L,
    seed = 20260524,
    ...) {
  results <- vector("list", length(scenarios))
  names(results) <- scenarios

  for (i in seq_along(scenarios)) {
    scenario_data <- ci_sim_add_outcome(
      base_data,
      scenario = scenarios[[i]],
      seed = seed + i
    )
    fit <- ci_sim_fit_models(
      data = scenario_data,
      types = types,
      ntree = ntree,
      seed = seed + i,
      ...
    )
    results[[i]] <- list(
      data = scenario_data,
      fits = fit$fits,
      summary = fit$summary
    )
  }

  results
}

ci_sim_collect_summary <- function(results) {
  data.table::rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
}

ci_sim_collect_importance <- function(
    results,
    scenarios = names(results),
    model = "forest",
    type = "L",
    repeats = 3L,
    seed = 20260524) {
  data.table::rbindlist(lapply(scenarios, function(scenario) {
    fit_name <- paste(type, model, sep = "_")
    importance <- ci_sim_permutation_importance(
      fit = results[[scenario]]$fits[[fit_name]],
      model = model,
      data = results[[scenario]]$data,
      repeats = repeats,
      seed = seed
    )
    importance[, `:=`(
      scenario = scenario,
      model = model,
      type = type
    )]
    data.table::setcolorder(importance, c("scenario", "model", "type", "feature"))
    importance
  }), fill = TRUE)
}
