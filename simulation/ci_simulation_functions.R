# The simulation has two scenarios:
# 1. null: wealth, outcome, and subgroup variables are generated independently.
#    This gives a negative-control setting where the tree should not find a
#    meaningful subgroup structure.
# 2. poor_subgroup_concentration: explicit subgroups are made poorer and given
#    a higher binary-outcome probability. This gives a positive-control setting
#    where the expected subgroup pattern is known before the tree is fitted.

CI_SIM_PREDICTORS <- c("mother_education", "residence", "child_sex", "province")
CI_SIM_PROVINCES <- c("Kinshasa", "Kongo Central", "Kasai", "Nord-Kivu")

#' Generate the observed covariates used in the simulation.
#'
#' This function creates the variables that the concentration-index tree is
#' allowed to split on. Wealth and the health outcome are not added here because
#' their relationship to the covariates is scenario-specific.
#'
#' @param n Number of observations to generate.
#' @param seed Random seed used for the covariate draws.
#'
#' @return A `data.table` with province, residence, maternal education, child
#'   sex, and unit sample weights.
ci_sim_make_base_data <- function(
  n = 1000L,
  seed = 20260602
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  set.seed(seed)

  # Start with observed covariates only. Wealth and the health outcome are added
  # later because their relationship to the covariates differs by simulation
  # scenario.
  dt <- data.table::data.table(
    province = sample(
      CI_SIM_PROVINCES,
      size = n,
      replace = TRUE,
      prob = c(0.25, 0.25, 0.25, 0.25)
    ),
    residence = sample(
      c("Urban", "Rural"),
      size = n,
      replace = TRUE,
      prob = c(0.35, 0.65)
    ),
    mother_education = sample(
      c("No education", "Primary", "Secondary or higher"),
      size = n,
      replace = TRUE,
      prob = c(0.30, 0.40, 0.30)
    ),
    child_sex = sample(c("Female", "Male"), size = n, replace = TRUE),
    sample_weight = 1
  )

  dt[, `:=`(
    mother_education = factor(
      mother_education,
      levels = c("Secondary or higher", "Primary", "No education")
    ),
    residence = factor(residence, levels = c("Urban", "Rural")),
    child_sex = factor(child_sex, levels = c("Female", "Male")),
    province = factor(province, levels = CI_SIM_PROVINCES)
  )]

  data.table::setcolorder(
    dt,
    c(CI_SIM_PREDICTORS, "sample_weight")
  )
  dt
}

#' Add wealth, subgroup labels, and a binary health outcome.
#'
#' This function turns the base covariate table into either a negative-control
#' null dataset or a positive-control poor-subgroup dataset. In the null
#' scenario, wealth and the outcome are independent of the subgroup variables.
#' In the poor-subgroup scenario, selected covariate combinations are placed
#' lower in the socioeconomic ranking and given higher outcome probabilities.
#'
#' @param data Base dataset created by `ci_sim_make_base_data()`.
#' @param scenario Simulation design. Use `"null"` for no planted subgroup
#'   signal, or `"poor_subgroup_concentration"` for the positive-control
#'   dataset.
#' @param seed Random seed used for wealth and outcome generation.
#' @param outcome_name Name of the binary outcome column to create.
#' @param baseline_risk Baseline probability of the adverse health outcome
#'   before any planted subgroup increase is added.
#' @param poor_share Proportion of observations treated as poor after ranking
#'   wealth from lowest to highest.
#' @param concentration_strength Multiplier for the planted outcome-risk
#'   increases. Larger values make the outcome more strongly concentrated among
#'   poor planted subgroups.
#' @param high_risk_provinces Province levels used to define the high-risk
#'   subgroup.
#' @param high_risk_residence Residence level used to define the high-risk
#'   subgroup.
#' @param high_risk_education Maternal-education level used to define the
#'   high-risk subgroup.
#' @param max_probability Upper bound applied to the simulated outcome
#'   probability.
#'
#' @return A `data.table` containing the original covariates, simulated wealth,
#'   socioeconomic rank indicators, planted subgroup labels, the binary outcome,
#'   and the true outcome probability.
ci_sim_add_binary_outcome <- function(
  data,
  scenario = c("null", "poor_subgroup_concentration"),
  seed = 20260602,
  outcome_name = "health_outcome",
  baseline_risk = 0.06,
  poor_share = 0.40,
  concentration_strength = 1.00,
  high_risk_provinces = c("Kasai", "Nord-Kivu"),
  high_risk_residence = "Rural",
  high_risk_education = "No education",
  max_probability = 0.75
) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  scenario <- match.arg(scenario)
  dt <- data.table::copy(data.table::as.data.table(data))

  concentration_strength <- max(0, as.numeric(concentration_strength))
  poor_share <- min(max(as.numeric(poor_share), 0.01), 0.99)

  # In the null scenario, socioeconomic position and mortality are independent
  # of the subgroup variables. This is the setting where a fitted tree should
  # have little reason to split the data into interpretable risk groups.
  if (identical(scenario, "null")) {
    set.seed(seed)
    dt[, wealth := pmax(stats::rnorm(.N, mean = 1.00, sd = 0.25), 0.01)]
    high_risk_subgroup <- rep(FALSE, nrow(dt))
    planted_subgroup <- rep("No planted subgroup", nrow(dt))
    probability <- rep(baseline_risk, nrow(dt))
  } else {
    # In the positive-control scenario, the subgroup structure is planted before
    # the outcome is generated. The tree can only recover this pattern through
    # the observed covariates: province, residence, and maternal education.
    planted_subgroup <- data.table::fcase(
      dt$province %in% high_risk_provinces &
        dt$residence == high_risk_residence &
        dt$mother_education == high_risk_education,
      "Rural, no education, high-risk province",
      dt$province %in% high_risk_provinces &
        dt$residence == high_risk_residence,
      "Rural high-risk province",
      dt$residence == high_risk_residence &
        dt$mother_education == high_risk_education,
      "Rural, no education",
      default = "Reference"
    )

    # The same planted groups are placed lower in the socioeconomic ranking.
    # This makes the health outcome concentrated among poorer individuals within
    # recognisable subgroups, which is the structure the concentration-index
    # tree is designed to separate.
    wealth_mean <- data.table::fcase(
      planted_subgroup == "Rural, no education, high-risk province", 0.25,
      planted_subgroup == "Rural high-risk province", 0.65,
      planted_subgroup == "Rural, no education", 0.80,
      default = 1.60
    )

    set.seed(seed)
    dt[, wealth := pmax(stats::rnorm(.N, mean = wealth_mean, sd = 0.25), 0.01)]

    high_risk_subgroup <- planted_subgroup != "Reference"
  }

  wealth_rank <- (data.table::frankv(dt$wealth, ties.method = "average") - 0.5) /
    nrow(dt)
  poor <- wealth_rank <= poor_share

  # The concentration_strength parameter controls how strongly the binary
  # outcome is concentrated in the poor planted subgroups. Increasing it should
  # make the recovered subgroup structure easier to see.
  if (identical(scenario, "poor_subgroup_concentration")) {
    probability <- baseline_risk + concentration_strength * data.table::fcase(
      planted_subgroup == "Rural, no education, high-risk province" & poor, 0.55,
      planted_subgroup == "Rural, no education, high-risk province", 0.35,
      planted_subgroup == "Rural high-risk province" & poor, 0.30,
      planted_subgroup == "Rural, no education" & poor, 0.25,
      poor, 0.08,
      default = 0.00
    )
  }

  probability <- pmin(pmax(probability, 0.001), max_probability)

  set.seed(seed + 1L)
  dt[, (outcome_name) := stats::rbinom(.N, size = 1L, prob = probability)]

  # Keep the true simulation quantities in the data. They are not used as split
  # variables, but they make it possible to check whether the fitted tree
  # recovers the intended data-generating structure.
  dt[, `:=`(
    wealth_rank = wealth_rank,
    poor = poor,
    planted_subgroup = factor(
      planted_subgroup,
      levels = c(
        "Reference",
        "Rural, no education",
        "Rural high-risk province",
        "Rural, no education, high-risk province",
        "No planted subgroup"
      )
    ),
    high_risk_subgroup = high_risk_subgroup,
    poor_high_risk_subgroup = poor & high_risk_subgroup,
    wealth_group = factor(
      ifelse(poor, "Poorer group", "Richer group"),
      levels = c("Poorer group", "Richer group")
    ),
    true_probability = probability,
    scenario = scenario,
    concentration_strength = concentration_strength
  )]

  data.table::setcolorder(
    dt,
    c(
      CI_SIM_PREDICTORS,
      "wealth",
      "wealth_rank",
      "wealth_group",
      "poor",
      "planted_subgroup",
      "high_risk_subgroup",
      "poor_high_risk_subgroup",
      outcome_name,
      "true_probability",
      "scenario",
      "concentration_strength",
      "sample_weight"
    )
  )
  dt
}

#' Generate one complete simulation dataset.
#'
#' This is a convenience wrapper that first creates the observed covariates and
#' then adds scenario-specific wealth and outcome values.
#'
#' @inheritParams ci_sim_make_base_data
#' @inheritParams ci_sim_add_binary_outcome
#'
#' @return A complete simulated `data.table` for one scenario.
ci_sim_make_simple_data <- function(
  n = 1000L,
  scenario = c("null", "poor_subgroup_concentration"),
  seed = 20260602,
  outcome_name = "health_outcome",
  baseline_risk = 0.06,
  poor_share = 0.40,
  concentration_strength = 1.00,
  high_risk_provinces = c("Kasai", "Nord-Kivu"),
  high_risk_residence = "Rural",
  high_risk_education = "No education",
  max_probability = 0.75
) {
  scenario <- match.arg(scenario)
  base_data <- ci_sim_make_base_data(n = n, seed = seed)
  ci_sim_add_binary_outcome(
    data = base_data,
    scenario = scenario,
    seed = seed + 1L,
    outcome_name = outcome_name,
    baseline_risk = baseline_risk,
    poor_share = poor_share,
    concentration_strength = concentration_strength,
    high_risk_provinces = high_risk_provinces,
    high_risk_residence = high_risk_residence,
    high_risk_education = high_risk_education,
    max_probability = max_probability
  )
}

#' Generate the null and poor-subgroup simulation datasets.
#'
#' This function returns the two datasets used in the simulation report: a null
#' negative-control dataset and a positive-control dataset where the outcome is
#' concentrated among poor planted subgroups.
#'
#' @inheritParams ci_sim_make_simple_data
#'
#' @return A named list with `null` and `poor_subgroup_concentration`
#'   `data.table` objects.
ci_sim_make_simple_examples <- function(
  n = 1000L,
  seed = 20260602,
  outcome_name = "health_outcome",
  baseline_risk = 0.06,
  poor_share = 0.40,
  concentration_strength = 1.00,
  high_risk_provinces = c("Kasai", "Nord-Kivu"),
  high_risk_residence = "Rural",
  high_risk_education = "No education",
  max_probability = 0.75
) {
  list(
    null = ci_sim_make_simple_data(
      n = n,
      scenario = "null",
      seed = seed,
      outcome_name = outcome_name,
      baseline_risk = baseline_risk,
      poor_share = poor_share,
      concentration_strength = concentration_strength,
      high_risk_provinces = high_risk_provinces,
      high_risk_residence = high_risk_residence,
      high_risk_education = high_risk_education,
      max_probability = max_probability
    ),
    poor_subgroup_concentration = ci_sim_make_simple_data(
      n = n,
      scenario = "poor_subgroup_concentration",
      seed = seed + 1000L,
      outcome_name = outcome_name,
      baseline_risk = baseline_risk,
      poor_share = poor_share,
      concentration_strength = concentration_strength,
      high_risk_provinces = high_risk_provinces,
      high_risk_residence = high_risk_residence,
      high_risk_education = high_risk_education,
      max_probability = max_probability
    )
  )
}
