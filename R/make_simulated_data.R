# Generate the positive-control simulated dataset used in the public report.

source(file.path("R", "ci_known_subgroup_demo.R"))

set.seed(20260528)

known_subgroup_sim_dt <- make_known_subgroup_data(
  n = 5000L,
  seed = 20260528
)

known_subgroup_fit <- fit_known_subgroup_ci_tree(
  data = known_subgroup_sim_dt,
  type = "L"
)

known_subgroup_validation <- validate_known_subgroup_recovery(
  fit = known_subgroup_fit,
  data = known_subgroup_sim_dt
)

known_subgroup_summary <- summarise_known_subgroup_tree(
  fit = known_subgroup_fit,
  data = known_subgroup_sim_dt
)

if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}

save(
  known_subgroup_sim_dt,
  known_subgroup_validation,
  known_subgroup_summary,
  file = file.path("data", "known_subgroup_sim_data.rda")
)

message("Saved data/known_subgroup_sim_data.rda")
