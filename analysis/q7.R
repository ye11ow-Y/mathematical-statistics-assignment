# Run from the assignment root: Rscript analysis/q7.R
# Conditional predictive modelling with fully disclosed assumptions.
# Base R and stats only. No figures are created.
options(digits = 12, warn = 2)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
N_students <- 146L
master_seed <- 20261108L
simulation_B <- 10000L
out_dir <- file.path("output", "q7")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, filename) write.csv(x, file.path(out_dir, filename), row.names = FALSE)

# Complete transcription of the displayed fragments. No omitted values are imputed.
fragments <- list(
  c(114.02, 94.02, 72.68, 51.18, 50.85, 46.47, 46.00, 45.73, 44.00, 43.95, 42.80, 40.95),
  c(20.78, 20.18),
  c(18.38, 17.77, 17.05, 16.78),
  c(3.32, 3.28),
  c(1.47, 1.43, 1.42, 1.20, 1.10),
  c(0.35, 0.33, 0.28, 0.25, 0.22, 0.10, 0.10, 0.07, 0.02, -0.07, -0.07, -2.98)
)
history <- do.call(rbind, lapply(seq_along(fragments), function(j) {
  data.frame(fragment = j, within_fragment = seq_along(fragments[[j]]),
             hours = fragments[[j]], ellipsis_after = FALSE)
}))
history$display_index <- seq_len(nrow(history))
history$ellipsis_after[cumsum(lengths(fragments))[-length(fragments)]] <- TRUE
save_csv(history, "historical_fragments.csv")
provenance <- data.frame(
  item = c("current_N", "deadline", "relative_time", "historical_values", "historical_omissions",
           "historical_N", "historical_window", "scenario_parameters", "scenario_weights", "conditional_iid"),
  value = c(as.character(N_students), "2026-11-12 17:25 Asia/Shanghai", "X=(deadline-submission)/1 hour; descending",
            "37 shown values; preserve repeated and negative values", "5 ellipses; counts unknown",
            "unknown", "2025-10-31 14:00 to 2025-11-14 17:25", "explicit assumptions in scenarios.csv",
            "central=.5; earlier=.25; later=.25", "assumption; one scenario shared by whole cohort"),
  source = c("user; confirmed", rep("exam", 6), rep("modelling assumption", 3))
)
save_csv(provenance, "provenance.csv")

# Within each scenario: late negative exponential; near-deadline power density
# on [0,24]; early shifted exponential above 24 hours.
scenarios <- data.frame(
  scenario = c("central", "earlier", "later"),
  weight = c(0.5, 0.25, 0.25),
  p_late = c(0.02, 0.01, 0.04),
  p_early = c(0.20, 0.30, 0.12),
  a_near = c(0.50, 0.75, 0.35),
  tau_late = c(1, 0.75, 1.5),
  tau_early = c(24, 24, 24),
  threshold_hours = 24
)
scenarios$p_near <- 1 - scenarios$p_late - scenarios$p_early
save_csv(scenarios, "scenarios.csv")

scenario_cdf <- function(x, row) {
  result <- numeric(length(x))
  late <- x < 0
  near <- x >= 0 & x <= row$threshold_hours
  early <- x > row$threshold_hours
  result[late] <- row$p_late * exp(x[late] / row$tau_late)
  result[near] <- row$p_late + row$p_near * (x[near] / row$threshold_hours)^row$a_near
  result[early] <- 1 - row$p_early * exp(-(x[early] - row$threshold_hours) / row$tau_early)
  pmin(1, pmax(0, result))
}

scenario_quantile <- function(p, row) {
  result <- numeric(length(p))
  late <- p < row$p_late
  near <- p >= row$p_late & p <= 1 - row$p_early
  early <- p > 1 - row$p_early
  result[late] <- row$tau_late * log(p[late] / row$p_late)
  result[near] <- row$threshold_hours * ((p[near] - row$p_late) / row$p_near)^(1 / row$a_near)
  result[early] <- row$threshold_hours - row$tau_early * log((1 - p[early]) / row$p_early)
  result
}

# A separate simple quantile baseline. Hours use displayed historical scales;
# all probability positions, the zero anchor and finite endpoints are assumptions.
baseline <- data.frame(
  probability = c(0, 0.02, 0.20, 0.35, 0.75, 0.92, 0.98, 0.995, 1),
  hours = c(-3.5, 0, 1.47, 3.32, 20.78, 40.95, 72.68, 114.02, 130),
  hour_source = c("assumed endpoint", "deadline", rep("displayed historical value", 6), "assumed endpoint")
)
save_csv(baseline, "baseline_quantile_knots.csv")
baseline_quantile <- function(p) approx(baseline$probability, baseline$hours, xout = p, rule = 2)$y

# For descending rank i, ascending rank is j=N+1-i and
# F_s(X_[i]) ~ Beta(j,i). A cohort-level scenario mixture must be taken AFTER
# this Beta CDF transformation, rather than mixing student CDFs first.
rank <- seq_len(N_students)
ascending_rank <- N_students + 1L - rank
beta_median <- qbeta(0.5, ascending_rank, rank)
component_medians <- sapply(seq_len(nrow(scenarios)), function(s) {
  scenario_quantile(beta_median, scenarios[s, ])
})
colnames(component_medians) <- scenarios$scenario
predictive_order_cdf <- function(x, i) {
  sum(vapply(seq_len(nrow(scenarios)), function(s) {
    scenarios$weight[s] * pbeta(scenario_cdf(x, scenarios[s, ]), N_students + 1L - i, i)
  }, numeric(1)))
}
final_prediction <- vapply(rank, function(i) {
  uniroot(function(x) predictive_order_cdf(x, i) - 0.5,
          interval = c(-100, 2000), tol = 1e-10)$root
}, numeric(1))
prediction <- data.frame(
  rank = rank, N_students = N_students,
  hours = final_prediction, hours_rounded = round(final_prediction, 2),
  beta_probability_median = beta_median,
  baseline_order_median = baseline_quantile(beta_median),
  central_order_median = component_medians[, "central"],
  earlier_order_median = component_medians[, "earlier"],
  later_order_median = component_medians[, "later"],
  predictive_cdf_at_prediction = vapply(rank, function(i) predictive_order_cdf(final_prediction[i], i), numeric(1))
)
save_csv(prediction, "predictions.csv")
save_csv(data.frame(rank = rank, hours = sprintf("%.2f", prediction$hours_rounded)), "predictions_for_submission.csv")
writeLines(paste(sprintf("%.2f", prediction$hours_rounded), collapse = ", "), file.path(out_dir, "prediction_values.txt"))

# Sampling-reference envelopes for the later plot, plus scenario medians.
# These describe hypothetical conditional models, not calibration to real data.
reference_bands <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(s) {
  data.frame(rank = rank, scenario = scenarios$scenario[s],
             lower90 = scenario_quantile(qbeta(0.05, ascending_rank, rank), scenarios[s, ]),
             median = component_medians[, s],
             upper90 = scenario_quantile(qbeta(0.95, ascending_rank, rank), scenarios[s, ]))
}))
save_csv(reference_bands, "scenario_order_reference.csv")
scenario_probabilities <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(s) {
  row <- scenarios[s, ]
  data.frame(scenario = row$scenario,
             p_late = row$p_late,
             p_within1h_before = scenario_cdf(1, row) - row$p_late,
             p_within6h_before = scenario_cdf(6, row) - row$p_late,
             p_within24h_before = row$p_near,
             p_earlier_than24h = row$p_early,
             expected_late_count = N_students * row$p_late,
             expected_earlier24_count = N_students * row$p_early)
}))
save_csv(scenario_probabilities, "scenario_probabilities.csv")

# Model-internal comparison: generate complete hypothetical cohorts under each
# named scenario, then compute W1. These are not observed forecasting errors.
candidate_predictions <- cbind(
  baseline = prediction$baseline_order_median,
  central = prediction$central_order_median,
  scenario_predictive_median = prediction$hours
)
model_risks <- list()
distance_rows <- list()
seeds <- data.frame(scenario = scenarios$scenario, seed = master_seed + 701L:703L,
                    repetitions = simulation_B, N_students = N_students)
save_csv(seeds, "seeds.csv")
for (s in seq_len(nrow(scenarios))) {
  set.seed(seeds$seed[s])
  losses <- matrix(NA_real_, nrow = simulation_B, ncol = ncol(candidate_predictions),
                   dimnames = list(NULL, colnames(candidate_predictions)))
  for (b in seq_len(simulation_B)) {
    actual_hypothetical <- sort(scenario_quantile(runif(N_students), scenarios[s, ]), decreasing = TRUE)
    losses[b, ] <- colMeans(abs(candidate_predictions - actual_hypothetical))
  }
  for (k in seq_len(ncol(losses))) {
    model_risks[[length(model_risks) + 1L]] <- data.frame(
      generating_scenario = scenarios$scenario[s], candidate = colnames(losses)[k],
      hypothetical_expected_W1 = mean(losses[, k]),
      mcse = sd(losses[, k]) / sqrt(simulation_B), repetitions = simulation_B)
  }
  gain <- losses[, "central"] - losses[, "scenario_predictive_median"]
  distance_rows[[s]] <- data.frame(
    scenario = scenarios$scenario[s],
    W1_from_final_to_component_median = mean(abs(final_prediction - component_medians[, s])),
    central_minus_final_risk = mean(gain), paired_mcse = sd(gain) / sqrt(simulation_B))
}
model_risks <- do.call(rbind, model_risks)
scenario_sensitivity <- do.call(rbind, distance_rows)
save_csv(model_risks, "hypothetical_risk_comparison.csv")
save_csv(scenario_sensitivity, "scenario_sensitivity.csv")
weighted_risks <- do.call(rbind, lapply(colnames(candidate_predictions), function(candidate) {
  rows <- model_risks[model_risks$candidate == candidate, ]
  data.frame(candidate = candidate,
             weighted_hypothetical_W1 = sum(scenarios$weight * rows$hypothetical_expected_W1),
             mcse = sqrt(sum((scenarios$weight * rows$mcse)^2)))
}))
save_csv(weighted_risks, "weighted_hypothetical_risks.csv")
summary <- data.frame(
  N_students = N_students, earliest_hours = prediction$hours_rounded[1],
  latest_hours = tail(prediction$hours_rounded, 1),
  positive_count_rounded = sum(prediction$hours_rounded > 0),
  zero_count_rounded = sum(prediction$hours_rounded == 0),
  negative_count_rounded = sum(prediction$hours_rounded < 0),
  within1h_before_rounded = sum(prediction$hours_rounded >= 0 & prediction$hours_rounded <= 1),
  within6h_before_rounded = sum(prediction$hours_rounded >= 0 & prediction$hours_rounded <= 6),
  within24h_before_rounded = sum(prediction$hours_rounded >= 0 & prediction$hours_rounded <= 24),
  earlier_than24h_rounded = sum(prediction$hours_rounded > 24),
  median_hours_rounded = median(prediction$hours_rounded),
  rounding_W1 = mean(abs(prediction$hours_rounded - prediction$hours)))
save_csv(summary, "prediction_summary.csv")

# Full Word-friendly table, four rank/value pairs per row. Read left to right,
# then proceed to the next row. No truncation or ellipsis is used.
table_lines <- c("| 排序 | 预测小时 | 排序 | 预测小时 | 排序 | 预测小时 | 排序 | 预测小时 |",
                 "|---:|---:|---:|---:|---:|---:|---:|---:|")
for (start in seq.int(1L, N_students, by = 4L)) {
  indices <- start:min(start + 3L, N_students)
  cells <- unlist(lapply(indices, function(i) c(as.character(i), sprintf("%.2f", prediction$hours_rounded[i]))))
  cells <- c(cells, rep("", 8L - length(cells)))
  table_lines <- c(table_lines, paste0("| ", paste(cells, collapse = " | "), " |"))
}
writeLines(table_lines, file.path(out_dir, "complete_prediction_table.md"), useBytes = TRUE)

stopifnot(nrow(history) == 37L, sum(history$ellipsis_after) == 5L,
          all(diff(history$hours) <= 0),
          N_students > 0, nrow(prediction) == N_students,
          all(is.finite(prediction$hours)), all(diff(prediction$hours) <= 1e-8),
          all(diff(prediction$hours_rounded) <= 0),
          max(abs(prediction$predictive_cdf_at_prediction - 0.5)) < 1e-7,
          summary$rounding_W1 <= 0.005 + 1e-10,
          abs(sum(scenarios$weight) - 1) < 1e-12,
          all(scenarios$p_near > 0))
for (s in seq_len(nrow(scenarios))) {
  probabilities <- seq(0.0001, 0.9999, length.out = 1001L)
  stopifnot(max(abs(scenario_cdf(scenario_quantile(probabilities, scenarios[s, ]), scenarios[s, ]) - probabilities)) < 1e-10,
            max(abs(scenario_cdf(c(0, 24), scenarios[s, ]) - c(scenarios$p_late[s], 1 - scenarios$p_early[s]))) < 1e-12)
}
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
capture.output({
  cat("N=146 is confirmed by the user. Model probabilities are assumptions.\n")
  cat("Historical values shown:", nrow(history), "; ellipses:", sum(history$ellipsis_after), "\n")
  cat("\nScenario assumptions\n"); print(scenarios)
  cat("\nPrediction summary\n"); print(summary)
  cat("\nSensitivity across hypothetical scenarios\n"); print(scenario_sensitivity)
  cat("\nModel-internal W1 comparison, not real forecasting accuracy\n"); print(model_risks)
  cat("\nWeighted model-internal risks\n"); print(weighted_risks)
  cat("\nAll predicted hours, descending\n"); cat(paste(sprintf("%.2f", prediction$hours_rounded), collapse = ", "), "\n")
  cat("\nVALIDATION PASSED: history transcription counts, N, order, finite values, median equations and rounding.\n")
}, file = file.path(out_dir, "calculation_summary.txt"))
cat("All Q7 calculations and validation passed; results saved to", out_dir, "\n")
