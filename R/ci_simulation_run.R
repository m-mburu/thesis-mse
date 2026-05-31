# Run the concentration-index tree and forest simulation scenarios.
# The Quarto report can run the same workflow directly; this script is useful
# when you want a saved RDS object for later inspection.

library(data.table)
library(ineqTrees)

source(here::here("R", "ci_simulation_functions.R"))

set.seed(20260524)

load(here::here("data", "drc_data.rda"))

simulation_base_dt <- ci_sim_prepare_congo_base(
  drc_data = drc_data,
  n = 2000L,
  use_full = FALSE,
  seed = 20260524
)

simulation_tree_control <- ci_tree_control(
  minsplit = 250L,
  minbucket = 125L,
  maxdepth = 4L,
  min_gain = 0,
  min_relative_gain = 0.05
)

simulation_results <- ci_sim_run_scenarios(
  base_data = simulation_base_dt,
  scenarios = c(
    "null",
    "simple_subgroup",
    "rank_mechanism",
    "level_mechanism",
    "complex_interaction"
  ),
  types = c("CI", "CIg", "CIc", "L"),
  tree_control = simulation_tree_control,
  forest_control = simulation_tree_control,
  ntree = 30L,
  mtry = max(1L, floor(sqrt(length(ci_sim_predictors())))),
  seed = 20260524
)

simulation_summary <- ci_sim_collect_summary(simulation_results)

simulation_importance <- rbindlist(list(
  ci_sim_collect_importance(
    simulation_results,
    model = "forest",
    type = "CI",
    repeats = 2L,
    seed = 20260525
  ),
  ci_sim_collect_importance(
    simulation_results,
    model = "forest",
    type = "L",
    repeats = 2L,
    seed = 20260524
  )
), fill = TRUE)

out <- list(
  base_data = simulation_base_dt,
  results = simulation_results,
  summary = simulation_summary,
  importance = simulation_importance
)

saveRDS(out, here::here("analysis", "ci_simulation_results.rds"))

