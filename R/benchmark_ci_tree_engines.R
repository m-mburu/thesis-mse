# Benchmark ineqTrees CI tree split engines on the saved DRC analysis object.

suppressPackageStartupMessages({
  library(data.table)
  library(ineqTrees)
  library(microbenchmark)
  library(here)
})

env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = NA_character_)))
  if (is.na(value) || value <= 0L) default else value
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

benchmark_sample <- function(data, n, seed = 20260601L) {
  data <- as.data.table(data)
  if (nrow(data) <= n) {
    return(copy(data))
  }

  set.seed(seed)
  data[sort(sample.int(nrow(data), n))]
}

fit_ci_tree_engine <- function(data, formula, type, engine) {
  control <- ineqTrees::ci_tree_control(
    minsplit = 500L,
    minbucket = 250L,
    minprob = 0.02,
    maxdepth = 5L,
    min_gain = 0.00001,
    min_relative_gain = 0.20,
    split_engine = engine
  )

  ineqTrees::ci_tree(
    formula = formula,
    data = data,
    rank_name = "wealth",
    outcome_name = "deadu5_num",
    weights = data$sample_weight,
    type = type,
    control = control
  )
}

summarise_benchmark <- function(bench) {
  out <- as.data.table(summary(bench, unit = "s"))
  out[, expr := as.character(expr)]
  out[
    ,
    c("criterion", "engine") := tstrsplit(expr, "_engine_", fixed = TRUE)
  ]

  wide <- dcast(
    out,
    criterion ~ engine,
    value.var = c("median", "mean", "min", "max")
  )

  if (all(c("median_R", "median_cpp") %in% names(wide))) {
    wide[
      ,
      `:=`(
        median_speedup_cpp_vs_R = median_R / median_cpp,
        median_time_reduction_pct = 100 * (median_R - median_cpp) / median_R
      )
    ]
  }

  list(summary = out[], comparison = wide[])
}

warm_up_engines <- function(expressions) {
  for (expr in expressions) {
    invisible(eval(expr, envir = parent.frame()))
  }
}

write_findings <- function(summary_dt, comparison_dt, log_file, csv_file, params) {
  fwrite(summary_dt, csv_file)

  lines <- c(
    "CI tree split-engine benchmark",
    paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Data rows benchmarked: ", format(params$n, big.mark = ",")),
    paste0("Repetitions per expression: ", params$times),
    paste0("Warm-up pass: ", if (isTRUE(params$warmup)) "yes" else "no"),
    paste0("Criteria: ", paste(params$types, collapse = ", ")),
    paste0("CSV summary: ", normalizePath(csv_file, winslash = "/", mustWork = FALSE)),
    "",
    "Median runtime comparison (seconds):"
  )

  if (nrow(comparison_dt)) {
    display_cols <- intersect(
      c(
        "criterion", "median_cpp", "median_R", "median_speedup_cpp_vs_R",
        "median_time_reduction_pct", "mean_cpp", "mean_R"
      ),
      names(comparison_dt)
    )
    lines <- c(
      lines,
      capture.output(print(comparison_dt[, ..display_cols], digits = 4))
    )

    if ("median_speedup_cpp_vs_R" %in% names(comparison_dt)) {
      finite_speedups <- comparison_dt[
        is.finite(median_speedup_cpp_vs_R),
        median_speedup_cpp_vs_R
      ]
      if (length(finite_speedups)) {
        lines <- c(
          lines,
          "",
          sprintf(
            "Finding: the cpp split engine was %.2fx faster on median across criteria.",
            stats::median(finite_speedups)
          )
        )
      }
    }
  } else {
    lines <- c(lines, "No benchmark rows were produced.")
  }

  writeLines(lines, log_file)
}

results_object_file <- here::here("data", "drc_report_results_objects.rda")
if (!file.exists(results_object_file)) {
  stop("Missing results object: ", results_object_file, call. = FALSE)
}

load(results_object_file)

bench_times <- env_int("CI_ENGINE_BENCH_TIMES", 3L)
bench_n <- env_int("CI_ENGINE_BENCH_N", 5000L)
bench_seed <- env_int("CI_ENGINE_BENCH_SEED", 20260601L)
bench_warmup <- tolower(env_chr("CI_ENGINE_BENCH_WARMUP", "true")) %in%
  c("1", "true", "yes", "y")
bench_types <- strsplit(env_chr("CI_ENGINE_BENCH_TYPES", "CI,CIg,CIc,L"), ",", fixed = TRUE)[[1]]
bench_types <- trimws(bench_types)
bench_types <- bench_types[nzchar(bench_types)]

bench_data <- benchmark_sample(
  drc_results_objects$congo_model_dt,
  n = bench_n,
  seed = bench_seed
)
bench_formula <- drc_results_objects$congo_ci_formula

expressions <- list()
for (type in bench_types) {
  for (engine in c("cpp", "R")) {
    label <- paste(type, "engine", engine, sep = "_")
    expressions[[label]] <- substitute(
      fit_ci_tree_engine(bench_data, bench_formula, TYPE, ENGINE),
      list(TYPE = type, ENGINE = engine)
    )
  }
}

if (isTRUE(bench_warmup)) {
  warm_up_engines(expressions)
}

bench <- microbenchmark::microbenchmark(
  list = expressions,
  times = bench_times,
  unit = "s"
)

bench_summary <- summarise_benchmark(bench)

log_dir <- here::here("logs")
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(log_dir, paste0("ci_tree_engine_benchmark_", timestamp, ".log"))
csv_file <- file.path(log_dir, paste0("ci_tree_engine_benchmark_", timestamp, ".csv"))

write_findings(
  summary_dt = bench_summary$summary,
  comparison_dt = bench_summary$comparison,
  log_file = log_file,
  csv_file = csv_file,
  params = list(
    n = nrow(bench_data),
    times = bench_times,
    warmup = bench_warmup,
    types = bench_types
  )
)

cat("Wrote benchmark log: ", normalizePath(log_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Wrote benchmark CSV: ", normalizePath(csv_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
