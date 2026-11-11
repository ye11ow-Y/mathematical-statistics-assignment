# Run from the assignment root: Rscript analysis/q4_q6.R
# Base R and stats. Plot-ready data are saved without creating figures.
options(digits = 12, warn = 2)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
out_dir <- file.path("output", "q4-q6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
master_seed <- 20261108L
normal_sizes <- c(4L, 16L, 64L, 256L, 1024L)
seeds <- data.frame(
  experiment = c("q5_bernoulli_path", "q6a_classification", "q6a_type1",
                 as.vector(rbind(paste0("q6b_H0_n", normal_sizes),
                                 paste0("q6b_H1_n", normal_sizes)))),
  seed = master_seed + c(51L, 61L, 62L, 71L:80L))
save_csv <- function(x, name) write.csv(x, file.path(out_dir, name), row.names = FALSE)
save_csv(seeds, "seeds.csv")
configuration <- data.frame(
  setting = c("q5_true_theta", "q5_n_grid", "q5_interval", "q6a_class_n_each",
              "q6a_class_rounds", "q6a_test_n", "q6a_test_B", "q6a_cutoff",
              "q6b_theta0", "q6b_theta1", "q6b_n", "q6b_B_each", "q6b_cutoff",
              "q6c_known", "q6c_symbolic"),
  value = c("0.3", "20,40,...,1000", "90% equal-tail", "10000", "1000",
            "100", "10000", "1.96", "0", "0.5", "4,16,64,256,1024",
            "10000", "1.65", "sample_mean=4; theta0=2.8", "n,sigma2,b2,theta1,model_weights"),
  source = c("exam", "exam", "level: exam; equal-tail: chosen", "exam", "chosen",
             "exam", "chosen", "exam", "chosen as planned", "chosen as planned",
             "exam", "chosen", "exam", "exam", "not specified by exam"))
save_csv(configuration, "configuration.csv")

# Q4 is analytic: eta = 2*sqrt(theta), eta' = 1/sqrt(theta),
# I_n(theta) = n/(2*theta^2), hence CRLB = 2*theta/n = eta^2/(2*n).

# Q5: all priors use the same prefix at each sample size.
set.seed(seeds$seed[1])
theta_true <- 0.3
bernoulli_path <- rbinom(1000L, size = 1, prob = theta_true)
successes <- cumsum(bernoulli_path)
save_csv(data.frame(i = seq_along(bernoulli_path), x = bernoulli_path,
                    cumulative_successes = successes), "q5_sample_path.csv")
priors <- data.frame(prior = c("Beta(3,7)", "Jeffreys", "Flat"),
                     a = c(3, 0.5, 1), b = c(7, 0.5, 1))
priors$prior_mean <- priors$a / (priors$a + priors$b)
priors$prior_strength <- priors$a + priors$b
save_csv(priors, "q5_priors.csv")
sample_sizes <- seq.int(20L, 1000L, by = 20L)
posterior_list <- lapply(seq_len(nrow(priors)), function(j) {
  a_post <- priors$a[j] + successes[sample_sizes]
  b_post <- priors$b[j] + sample_sizes - successes[sample_sizes]
  lower <- qbeta(0.05, a_post, b_post)
  upper <- qbeta(0.95, a_post, b_post)
  data.frame(prior = priors$prior[j], n = sample_sizes, successes = successes[sample_sizes],
             a_post = a_post, b_post = b_post,
             posterior_mean = a_post / (a_post + b_post),
             posterior_variance = a_post * b_post / ((a_post + b_post)^2 * (a_post + b_post + 1)),
             lower90 = lower, upper90 = upper, interval_width = upper - lower)
})
posterior <- do.call(rbind, posterior_list)
save_csv(posterior, "q5_posterior_path.csv")
save_csv(posterior[posterior$n %in% c(20, 100, 500, 1000), ], "q5_posterior_checkpoints.csv")
grid <- seq(0.00001, 0.99999, length.out = 1001L)
density_n20 <- do.call(rbind, lapply(seq_len(nrow(priors)), function(j) {
  row <- posterior_list[[j]][1, ]
  data.frame(prior = row$prior, n = 20, theta = grid,
             posterior_density = dbeta(grid, row$a_post, row$b_post),
             posterior_mean = row$posterior_mean, truth = theta_true)
}))
save_csv(density_n20, "q5_density_n20.csv")
bayes_risk <- data.frame(
  n = sample_sizes,
  posterior_mean_bayes_mse = 1 / (6 * (sample_sizes + 2)),
  classical_mean_bayes_mse = 1 / (6 * sample_sizes),
  bayesian_lower_bound = (1 / 6)^2 / (sample_sizes / 6 + 1 / 3),
  fixed_theta03_risk = (sample_sizes * theta_true * (1 - theta_true) + (1 - 2 * theta_true)^2) / (sample_sizes + 2)^2)
# Under a uniform prior the predictive law of T is uniform on 0,...,n.
bayes_risk$expected_posterior_variance <- vapply(sample_sizes, function(n) {
  t <- 0:n
  mean((t + 1) * (n - t + 1) / ((n + 2)^2 * (n + 3)))
}, numeric(1))
save_csv(bayes_risk, "q5_bayes_risk.csv")

# Probability summaries use Clopper-Pearson binomial intervals, including at 0 or 1.
probability_summary <- function(count, total, probability, label) {
  interval <- binom.test(count, total, conf.level = 0.95)$conf.int
  data.frame(experiment = label, count = count, total = total,
             empirical_probability = count / total, theoretical_probability = probability,
             mcse_theoretical = sqrt(probability * (1 - probability) / total),
             mc_lower95 = interval[1], mc_upper95 = interval[2])
}

# Q6A(1): actually generate both samples in each of 1000 rounds.
set.seed(seeds$seed[2])
rounds <- 1000L
class_n <- 10000L
counts <- matrix(0L, nrow = rounds, ncol = 2,
                 dimnames = list(NULL, c("positive_errors", "negative_errors")))
for (r in seq_len(rounds)) {
  counts[r, 1] <- sum(rnorm(class_n, mean = 4, sd = 1) <= 0)
  counts[r, 2] <- sum(rnorm(class_n, mean = -4, sd = 1) > 0)
}
class_probability <- pnorm(-4)
classification_rounds <- data.frame(round = seq_len(rounds), counts)
save_csv(classification_rounds, "q6a_classification_rounds.csv")
classification <- rbind(
  probability_summary(sum(counts[, 1]), rounds * class_n, class_probability, "N(4,1)"),
  probability_summary(sum(counts[, 2]), rounds * class_n, class_probability, "N(-4,1)"),
  probability_summary(sum(counts), 2 * rounds * class_n, class_probability, "pooled"))
save_csv(classification, "q6a_classification_summary.csv")
zero_error <- data.frame(
  event = c("one_group_zero_errors", "both_groups_zero_errors"),
  theoretical_probability = c((1 - class_probability)^class_n, (1 - class_probability)^(2 * class_n)),
  empirical_fraction_rounds = c(mean(counts[, 1] == 0), mean(rowSums(counts) == 0)))
save_csv(zero_error, "q6a_zero_error_rounds.csv")

# Q6A(2): one complete n=100 dataset corresponds to one decision.
set.seed(seeds$seed[3])
test_n <- 100L
test_B <- 10000L
cutoff_two_sided <- 1.96
test_data <- matrix(rnorm(test_n * test_B, mean = 4, sd = 1), nrow = test_n)
test_means <- colMeans(test_data)
test_z <- sqrt(test_n) * (test_means - 4)
test_reject <- abs(test_z) > cutoff_two_sided
type1 <- probability_summary(sum(test_reject), test_B,
                              2 * pnorm(-cutoff_two_sided), "two_sided_type1")
save_csv(data.frame(replication = seq_len(test_B), mean = test_means,
                    z = test_z, reject = as.integer(test_reject)), "q6a_type1_replications.csv")
save_csv(type1, "q6a_type1_summary.csv")
rm(test_data)
gc()

# Q6B: generate all n observations, under both H0 and H1, at every n.
theta0 <- 0
theta1 <- 0.5
test_B <- 10000L
cutoff_one_sided <- 1.65
np_replicates <- vector("list", 2L * length(normal_sizes))
np_summary <- vector("list", 2L * length(normal_sizes))
for (j in seq_along(normal_sizes)) {
  n <- normal_sizes[j]
  for (h in 0:1) {
    index <- 2L * (j - 1L) + h + 1L
    true_theta <- if (h == 0) theta0 else theta1
    set.seed(seeds$seed[index + 3L])
    x <- matrix(rnorm(n * test_B, mean = true_theta, sd = 1), nrow = n)
    xbar <- colMeans(x)
    z <- sqrt(n) * (xbar - theta0)
    reject <- z > cutoff_one_sided
    theoretical_probability <- pnorm(sqrt(n) * (true_theta - theta0) - cutoff_one_sided)
    np_summary[[index]] <- cbind(
      data.frame(n = n, hypothesis = paste0("H", h), true_theta = true_theta,
                 mean_threshold = theta0 + cutoff_one_sided / sqrt(n)),
      probability_summary(sum(reject), test_B, theoretical_probability,
                          if (h == 0) "size" else "power"))
    np_replicates[[index]] <- data.frame(n = n, hypothesis = paste0("H", h),
                                        replication = seq_len(test_B), mean = xbar,
                                        z = z, reject = as.integer(reject))
    rm(x)
    gc()
  }
}
np_summary <- do.call(rbind, np_summary)
np_replicates <- do.call(rbind, np_replicates)
save_csv(np_summary, "q6b_summary.csv")
save_csv(np_replicates, "q6b_replications.csv")
power_table <- data.frame(
  n = normal_sizes,
  theoretical_size = pnorm(-cutoff_one_sided),
  empirical_size = np_summary$empirical_probability[np_summary$hypothesis == "H0"],
  theoretical_power = pnorm(sqrt(normal_sizes) * (theta1 - theta0) - cutoff_one_sided),
  empirical_power = np_summary$empirical_probability[np_summary$hypothesis == "H1"],
  theoretical_type2 = pnorm(cutoff_one_sided - sqrt(normal_sizes) * (theta1 - theta0)))
save_csv(power_table, "q6b_power_table.csv")
curve_n <- 1:1024
save_csv(data.frame(n = curve_n,
                    power = pnorm(sqrt(curve_n) * (theta1 - theta0) - cutoff_one_sided),
                    type2 = pnorm(cutoff_one_sided - sqrt(curve_n) * (theta1 - theta0))), "q6b_theoretical_curve.csv")

# Numerical/algebraic validation and traceability.
stopifnot(nrow(posterior) == 150L, nrow(density_n20) == 3003L,
          all(posterior$lower90 < posterior$upper90),
          all(posterior$posterior_mean > 0 & posterior$posterior_mean < 1),
          max(abs(bayes_risk$posterior_mean_bayes_mse - bayes_risk$bayesian_lower_bound)) < 1e-14,
          max(abs(bayes_risk$posterior_mean_bayes_mse - bayes_risk$expected_posterior_variance)) < 1e-14,
          nrow(np_replicates) == 100000L,
          all(table(np_replicates$n, np_replicates$hypothesis) == 10000L))
for (j in seq_len(nrow(priors))) {
  rows <- posterior_list[[j]]
  stopifnot(max(abs(pbeta(rows$lower90, rows$a_post, rows$b_post) - 0.05)) < 1e-10,
            max(abs(pbeta(rows$upper90, rows$a_post, rows$b_post) - 0.95)) < 1e-10)
}
probability_checks <- rbind(classification, type1,
                            np_summary[, names(classification)])
stopifnot(all(abs(probability_checks$empirical_probability - probability_checks$theoretical_probability) <=
                5 * probability_checks$mcse_theoretical + 1 / probability_checks$total))
save_csv(probability_checks, "validation_probabilities.csv")
saveRDS(list(configuration = configuration, seeds = seeds, priors = priors,
             bernoulli_path = bernoulli_path, posterior = posterior,
             classification_rounds = classification_rounds,
             test_means = test_means, np_replicates = np_replicates),
        file.path(out_dir, "simulation_record.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
capture.output({
  cat("Q4: CRLB for eta=2*sigma is 2*sigma^2/n = eta^2/(2*n).\n")
  cat("\nQ5: first 20 observations\n"); print(bernoulli_path[1:20])
  cat("\nQ5: posterior checkpoints\n"); print(posterior[posterior$n %in% c(20, 100, 500, 1000), ])
  cat("\nQ5: risk checkpoints\n"); print(bayes_risk[bayes_risk$n %in% c(20, 100, 500, 1000), ])
  cat("\nQ6A: classification\n"); print(classification)
  cat("\nQ6A: zero-error rounds\n"); print(zero_error)
  cat("\nQ6A: Type I error\n"); print(type1)
  cat("\nQ6B: power table\n"); print(power_table)
  cat("\nQ6B: simulation uncertainty\n"); print(np_summary)
  cat("\nQ6C: symbolic model comparison; no omitted exam parameter has been assigned a value.\n")
  cat("\nVALIDATION PASSED: counts, Beta quantiles, exact Bayesian risk and simulation probabilities.\n")
}, file = file.path(out_dir, "calculation_summary.txt"))
cat("All calculations and validation passed; results saved to", out_dir, "\n")
