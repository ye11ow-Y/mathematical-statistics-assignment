# Run from the assignment root: Rscript analysis/q1_q3.R
# Base R only. No figures are produced; plot-ready data are saved instead.
options(digits = 12, warn = 2)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
out_dir <- file.path("output", "q1-q3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
master_seed <- 20261108L
seeds <- data.frame(
  experiment = c("q1_initial", "q1_repeated", "q1_path",
                 "q1_normal_n10", "q1_normal_n100", "q1_normal_n1000",
                 "q1_normal_n10000", "q2_shift", "q2_scale", "q3_rb"),
  seed = master_seed + seq_len(10L)
)
write.csv(seeds, file.path(out_dir, "seeds.csv"), row.names = FALSE)
config <- data.frame(
  experiment = seeds$experiment,
  n = c(100, 100, 10000, 10, 100, 1000, 10000, 10, 10, 10),
  repetitions = c(1, 1000, 1, rep(1000, 4), 10000, 10000, 10000),
  population_mean = c(rep(1, 7), 3, 3, 2),
  population_sd = c(rep(5, 7), 1, 3, 1),
  target = c("mu and sigma2", "mu and sigma2", "mu and sigma2",
             rep("mu", 4), "theta=3", "theta=3", "theta=1"),
  settings_source = c(rep("exam; seed chosen", 7),
                      rep("exam; B=10000 and seed chosen", 2),
                      "theta=1, n=10, B=10000 and seed chosen")
)
write.csv(config, file.path(out_dir, "configuration.csv"), row.names = FALSE)

# The empirical variance uses B; var() itself uses B-1.
summary_row <- function(est, truth, theory_mse, label) {
  b <- length(est)
  mse <- mean((est - truth)^2)
  variance_b <- mean((est - mean(est))^2)
  bias <- mean(est) - truth
  stopifnot(abs(mse - variance_b - bias^2) < 1e-10 * max(1, mse))
  data.frame(estimator = label, first_estimate = est[1],
             mean_estimate = mean(est), empirical_bias = bias,
             empirical_mse = mse, empirical_variance_B = variance_b,
             empirical_variance_Bminus1 = var(est), theoretical_mse = theory_mse,
             mse_mcse = sd((est - truth)^2) / sqrt(b), repetitions = b)
}
confidence_intervals <- function(xbar, s2, n, alpha = 0.05) {
  half <- qt(1 - alpha / 2, n - 1) * sqrt(s2 / n)
  data.frame(mu_lower = xbar - half, mu_upper = xbar + half,
             sigma2_lower = (n - 1) * s2 / qchisq(1 - alpha / 2, n - 1),
             sigma2_upper = (n - 1) * s2 / qchisq(alpha / 2, n - 1))
}
save_csv <- function(x, name) write.csv(x, file.path(out_dir, name), row.names = FALSE)

# Q1.1: both population parameters are treated as unknown in the intervals.
mu <- 1
sigma <- 5
n <- 100L
alpha <- 0.05
B <- 1000L
set.seed(seeds$seed[1])
x_initial <- rnorm(n, mean = mu, sd = sigma)
initial <- cbind(data.frame(n = n, mean = mean(x_initial),
                           variance_unbiased = var(x_initial),
                           variance_mle = mean((x_initial - mean(x_initial))^2)),
                 confidence_intervals(mean(x_initial), var(x_initial), n))
save_csv(initial, "q1_initial.csv")
save_csv(data.frame(i = seq_len(n), x = x_initial), "q1_initial_sample.csv")

# Q1.2(1)(2): each column is one independent sample of size 100.
set.seed(seeds$seed[2])
x_repeated <- matrix(rnorm(n * B, mean = mu, sd = sigma), nrow = n)
rep_means <- colMeans(x_repeated)
rep_vars <- apply(x_repeated, 2, var)
rep_intervals <- confidence_intervals(rep_means, rep_vars, n)
q1_repeated <- cbind(data.frame(replication = seq_len(B), mean = rep_means,
                               variance = rep_vars), rep_intervals)
save_csv(q1_repeated, "q1_repeated.csv")
q1_mse <- rbind(summary_row(rep_means, mu, sigma^2 / n, "sample_mean"),
                summary_row(rep_vars, sigma^2, 2 * sigma^4 / (n - 1), "sample_variance"))
save_csv(q1_mse, "q1_mse.csv")
coverage <- data.frame(
  parameter = c("mu", "sigma2"),
  coverage = c(mean(rep_intervals$mu_lower <= mu & mu <= rep_intervals$mu_upper),
               mean(rep_intervals$sigma2_lower <= sigma^2 & sigma^2 <= rep_intervals$sigma2_upper)),
  nominal = 1 - alpha, mcse_at_nominal = sqrt(alpha * (1 - alpha) / B))
save_csv(coverage, "q1_coverage.csv")

# Q1.2(3): a single length-10000 sequence; each n uses a prefix.
set.seed(seeds$seed[3])
x_path <- rnorm(10000L, mean = mu, sd = sigma)
n_grid <- c(seq.int(10L, 9990L, by = 20L), 10000L)
path_means <- vapply(n_grid, function(k) mean(x_path[seq_len(k)]), numeric(1))
path_vars <- vapply(n_grid, function(k) var(x_path[seq_len(k)]), numeric(1))
q1_path <- cbind(data.frame(n = n_grid, mean = path_means, variance = path_vars),
                 confidence_intervals(path_means, path_vars, n_grid))
# These are pointwise sampling reference bands around the true parameters,
# distinct from the confidence intervals above. They are NOT simultaneous bands.
q1_path$mean_reference_lower <- mu + qnorm(alpha / 2) * sigma / sqrt(n_grid)
q1_path$mean_reference_upper <- mu + qnorm(1 - alpha / 2) * sigma / sqrt(n_grid)
q1_path$variance_reference_lower <- sigma^2 * qchisq(alpha / 2, n_grid - 1) / (n_grid - 1)
q1_path$variance_reference_upper <- sigma^2 * qchisq(1 - alpha / 2, n_grid - 1) / (n_grid - 1)
save_csv(q1_path, "q1_consistency_path.csv")
save_csv(q1_path[match(c(10, 110, 1010, 10000), n_grid), ], "q1_path_checkpoints.csv")

# Q1.3: actually generate n observations in every replication.
# Do not substitute a direct draw of the sample mean's known distribution.
normal_sizes <- c(10L, 100L, 1000L, 10000L)
normal_raw <- vector("list", length(normal_sizes))
normal_summaries <- vector("list", length(normal_sizes))
for (j in seq_along(normal_sizes)) {
  k <- normal_sizes[j]
  set.seed(seeds$seed[j + 3L])
  estimates <- replicate(B, {
    x <- rnorm(k, mean = mu, sd = sigma)
    c(xbar = mean(x), s = sd(x))
  })
  z <- sqrt(k) * (estimates[1, ] - mu) / sigma
  t_stat <- sqrt(k) * (estimates[1, ] - mu) / estimates[2, ]
  normal_raw[[j]] <- data.frame(n = k, replication = seq_len(B),
                                mean = estimates[1, ], sd = estimates[2, ],
                                z_true_sigma = z, t_sample_sd = t_stat)
  normal_summaries[[j]] <- data.frame(
    n = k, z_mean = mean(z), z_variance_B = mean((z - mean(z))^2),
    z_q025 = unname(quantile(z, 0.025)), z_q975 = unname(quantile(z, 0.975)),
    z_fraction_within_196 = mean(abs(z) <= qnorm(0.975)),
    t_mean = mean(t_stat), t_variance_B = mean((t_stat - mean(t_stat))^2),
    t_theoretical_variance = (k - 1) / (k - 3))
}
q1_normality <- do.call(rbind, normal_raw)
q1_normality_summary <- do.call(rbind, normal_summaries)
save_csv(q1_normality, "q1_normality_samples.csv")
save_csv(q1_normality_summary, "q1_normality_summary.csv")

# Q2: use X = theta + Z for location, X = theta * Z for scale.
theta <- 3
n2 <- 10L
B2 <- 10000L
set.seed(seeds$seed[8])
z_shift <- matrix(rnorm(n2 * B2, mean = 0, sd = 1), nrow = n2)
x_shift <- theta + z_shift
shift_est <- colMeans(x_shift)
set.seed(seeds$seed[9])
z_scale <- matrix(rnorm(n2 * B2, mean = 1, sd = 1), nrow = n2)
x_scale <- theta * z_scale
scale_sample_mean <- colMeans(x_scale)
scale_est <- colSums(x_scale) / (n2 + 1)
q2_summary <- rbind(summary_row(shift_est, theta, 1 / n2, "shift_best_unbiased_linear"),
                    summary_row(scale_est, theta, theta^2 / (n2 + 1), "scale_best_linear"),
                    summary_row(scale_sample_mean, theta, theta^2 / n2, "scale_unbiased_mean"))
save_csv(q2_summary, "q2_summary.csv")
save_csv(data.frame(replication = seq_len(B2), shift = shift_est,
                    scale = scale_est, scale_sample_mean = scale_sample_mean), "q2_repeated.csv")
scale_gain <- (scale_sample_mean - theta)^2 - (scale_est - theta)^2
q2_gain <- data.frame(empirical_mse_reduction = mean(scale_gain),
                     theoretical_mse_reduction = theta^2 / n2 - theta^2 / (n2 + 1),
                     paired_mcse = sd(scale_gain) / sqrt(B2))
save_csv(q2_gain, "q2_scale_comparison.csv")

# Q3: paired comparison from the SAME sample in each replication.
theta3 <- 1
n3 <- 10L
B3 <- 10000L
set.seed(seeds$seed[10])
x_rb <- matrix(rnorm(n3 * B3, mean = 2 * theta3, sd = 1), nrow = n3)
original <- x_rb[1, ] / 2
improved <- colMeans(x_rb) / 2
q3_summary <- rbind(summary_row(original, theta3, 1 / 4, "X1_over_2"),
                    summary_row(improved, theta3, 1 / (4 * n3), "Xbar_over_2"))
save_csv(q3_summary, "q3_summary.csv")
save_csv(data.frame(replication = seq_len(B3), original = original,
                    rao_blackwell = improved), "q3_repeated.csv")
paired_gain <- (original - theta3)^2 - (improved - theta3)^2
q3_comparison <- data.frame(
  mse_ratio = q3_summary$empirical_mse[2] / q3_summary$empirical_mse[1],
  relative_reduction = 1 - q3_summary$empirical_mse[2] / q3_summary$empirical_mse[1],
  theoretical_ratio = 1 / n3, theoretical_reduction = 1 - 1 / n3,
  empirical_mse_gain = mean(paired_gain), theoretical_mse_gain = (n3 - 1) / (4 * n3),
  gain_mcse_paired = sd(paired_gain) / sqrt(B3))
save_csv(q3_comparison, "q3_comparison.csv")

# Validation: algebraic identities, counts, interval ordering and MC accuracy.
all_mse <- rbind(q1_mse, q2_summary, q3_summary)
all_mse$mcse_units_from_theory <- (all_mse$empirical_mse - all_mse$theoretical_mse) / all_mse$mse_mcse
stopifnot(all(abs(all_mse$mcse_units_from_theory) < 5),
          nrow(q1_path) == 501L, identical(tail(n_grid, 2), c(9990L, 10000L)),
          all(table(q1_normality$n) == 1000L),
          ncol(x_shift) == B2, ncol(x_scale) == B2, ncol(x_rb) == B3,
          all(rep_intervals$mu_lower < rep_intervals$mu_upper),
          all(rep_intervals$sigma2_lower < rep_intervals$sigma2_upper),
          all(abs(coverage$coverage - coverage$nominal) < 5 * coverage$mcse_at_nominal),
          max(abs(scale_est - n2 / (n2 + 1) * scale_sample_mean)) < 1e-12,
          max(abs(improved - colSums(x_rb) / (2 * n3))) < 1e-12,
          all(abs(q1_normality_summary$z_mean) < 5 / sqrt(B)),
          all(abs(q1_normality_summary$z_variance_B - (B - 1) / B) < 5 * sqrt(2 * (B - 1)) / B),
          abs(mean(paired_gain) - (n3 - 1) / (4 * n3)) < 5 * sd(paired_gain) / sqrt(B3))
save_csv(all_mse, "validation_mse.csv")
saveRDS(list(config = config, seeds = seeds, x_initial = x_initial,
             x_repeated = x_repeated, x_path = x_path,
             q1_normality = q1_normality, x_shift = x_shift,
             x_scale = x_scale, x_rb = x_rb), file.path(out_dir, "simulation_record.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
capture.output({
  cat("Q1 initial sample and intervals\n"); print(initial)
  cat("\nQ1 MSE and empirical variance\n"); print(q1_mse)
  cat("\nQ1 interval coverage (supplementary)\n"); print(coverage)
  cat("\nQ1 path checkpoints\n"); print(q1_path[match(c(10, 110, 1010, 10000), n_grid), 1:3])
  cat("\nQ1 standardisation\n"); print(q1_normality_summary)
  cat("\nQ2 summary\n"); print(q2_summary)
  cat("\nQ2 paired scale comparison\n"); print(q2_gain)
  cat("\nQ3 summary\n"); print(q3_summary)
  cat("\nQ3 paired comparison\n"); print(q3_comparison)
  cat("\nVALIDATION PASSED: identities, counts, intervals and Monte Carlo checks.\n")
}, file = file.path(out_dir, "calculation_summary.txt"))
cat("All calculations and validation passed; results saved to", out_dir, "\n")
