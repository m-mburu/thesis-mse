# Generate the simple simulation datasets used by simulation/simulation_fits.qmd.

library(data.table)

source(here::here("simulation", "ci_simulation_functions.R"))

simulation_data <- ci_sim_make_simple_examples(
  n = 3000L,
  seed = 20260603,
  outcome_name = "health_outcome",
  baseline_risk = 0.06,
  poor_share = 0.40,
  concentration_strength = 1.00
)

simulation_summary <- rbindlist(lapply(names(simulation_data), function(scenario_name) {
  dt <- simulation_data[[scenario_name]]
  data.table(
    scenario = scenario_name,
    n = nrow(dt),
    outcome_rate = mean(dt$health_outcome),
    poor_outcome_rate = mean(dt$health_outcome[dt$poor]),
    richer_outcome_rate = mean(dt$health_outcome[!dt$poor]),
    high_risk_subgroup_outcome_rate = mean(dt$health_outcome[dt$high_risk_subgroup]),
    poor_high_risk_subgroup_outcome_rate = mean(dt$health_outcome[dt$poor_high_risk_subgroup])
  )
}), fill = TRUE)

if (!dir.exists(here::here("analysis"))) {
  dir.create(here::here("analysis"), recursive = TRUE)
}

saveRDS(
  list(
    data = simulation_data,
    summary = simulation_summary
  ),
  here::here("analysis", "ci_simulation_results.rds")
)

print(simulation_summary)
