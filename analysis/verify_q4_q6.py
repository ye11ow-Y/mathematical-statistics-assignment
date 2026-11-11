"""Independently verify R output using NumPy and SciPy.

Run from the assignment root: python analysis/verify_q4_q6.py
"""

import csv
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.stats import beta, norm


RESULTS = Path(__file__).resolve().parents[1] / "output" / "q4-q6"
checks = []


def read(name):
    with (RESULTS / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def close(actual, expected):
    assert math.isclose(float(actual), float(expected), rel_tol=2e-9, abs_tol=2e-12), (
        actual, expected
    )


def whole(value):
    # R may export integer-valued totals using scientific notation.
    number = float(value)
    assert number.is_integer() and abs(number) < 2**53
    return int(number)


def exact_interval(count, total):
    lower = 0 if count == 0 else beta.ppf(0.025, count, total - count + 1)
    upper = 1 if count == total else beta.ppf(0.975, count + 1, total - count)
    return lower, upper


path = read("q5_sample_path.csv")
observations = np.array([whole(row["x"]) for row in path])
successes = observations.cumsum()
assert len(path) == 1000 and set(observations) == {0, 1}
assert successes.tolist() == [whole(row["cumulative_successes"]) for row in path]
posterior = read("q5_posterior_path.csv")
priors = {
    row["prior"]: (float(row["a"]), float(row["b"]))
    for row in read("q5_priors.csv")
}
assert len(posterior) == 150
assert Counter(whole(row["n"]) for row in posterior) == Counter(
    {n: 3 for n in range(20, 1001, 20)}
)
for row in posterior:
    n = whole(row["n"])
    count = whole(row["successes"])
    a, b = priors[row["prior"]]
    assert count == successes[n - 1]
    a_post, b_post = a + count, b + n - count
    close(row["a_post"], a_post)
    close(row["b_post"], b_post)
    close(row["posterior_mean"], a_post / (a_post + b_post))
    close(
        row["posterior_variance"],
        a_post * b_post / ((a_post + b_post)**2 * (a_post + b_post + 1)),
    )
    close(row["lower90"], beta.ppf(0.05, a_post, b_post))
    close(row["upper90"], beta.ppf(0.95, a_post, b_post))
    close(row["interval_width"], float(row["upper90"]) - float(row["lower90"]))
for row in read("q5_density_n20.csv"):
    a, b = priors[row["prior"]]
    close(
        row["posterior_density"],
        beta.pdf(float(row["theta"]), a + successes[19], b + 20 - successes[19]),
    )
checks.append("Q5: shared Bernoulli path, 150 posterior rows, 3003 density values and 90% Beta quantiles passed.")

J = quad(lambda theta: theta * (1 - theta), 0, 1)[0]
H = quad(lambda theta: (1 - 2 * theta)**2, 0, 1)[0]
for row in read("q5_bayes_risk.csv"):
    n = whole(row["n"])
    risk = quad(
        lambda theta: (n * theta * (1 - theta) + (1 - 2 * theta)**2) / (n + 2)**2,
        0, 1,
    )[0]
    close(row["posterior_mean_bayes_mse"], risk)
    close(row["bayesian_lower_bound"], J * J / (n * J + H))
    close(row["classical_mean_bayes_mse"], 1 / (6 * n))
    close(row["expected_posterior_variance"], risk)
    close(row["fixed_theta03_risk"], (0.21 * n + 0.16) / (n + 2)**2)
checks.append("Q5: Bayesian MSE, weighted lower bound, expected posterior variance and fixed-parameter risk passed.")

rounds = read("q6a_classification_rounds.csv")
assert len(rounds) == 1000
positive_errors = sum(whole(row["positive_errors"]) for row in rounds)
negative_errors = sum(whole(row["negative_errors"]) for row in rounds)
classification = read("q6a_classification_summary.csv")
assert [whole(row["count"]) for row in classification] == [
    positive_errors, negative_errors, positive_errors + negative_errors
]
for row in classification:
    close(row["theoretical_probability"], norm.cdf(-4))
type1_raw = read("q6a_type1_replications.csv")
assert len(type1_raw) == 10000
for row in type1_raw:
    close(row["z"], 10 * (float(row["mean"]) - 4))
    assert whole(row["reject"]) == int(abs(float(row["z"])) > 1.96)
type1_summary = read("q6a_type1_summary.csv")[0]
assert sum(whole(row["reject"]) for row in type1_raw) == whole(type1_summary["count"])
close(type1_summary["theoretical_probability"], 2 * norm.sf(1.96))
checks.append("Q6A: 1000 round error counts, pooled probability and 10000 complete test decisions passed.")

np_raw = read("q6b_replications.csv")
assert len(np_raw) == 100000
counts, totals = Counter(), Counter()
for row in np_raw:
    n = whole(row["n"])
    key = (n, row["hypothesis"])
    z = math.sqrt(n) * float(row["mean"])
    close(row["z"], z)
    assert whole(row["reject"]) == int(z > 1.65)
    counts[key] += whole(row["reject"])
    totals[key] += 1
for row in read("q6b_summary.csv"):
    n = whole(row["n"])
    key = (n, row["hypothesis"])
    assert totals[key] == 10000 and counts[key] == whole(row["count"])
    close(
        row["theoretical_probability"],
        norm.sf(1.65 - math.sqrt(n) * float(row["true_theta"])),
    )
for row in read("q6b_power_table.csv"):
    n = whole(row["n"])
    # A relative check retains accuracy in extremely small type-II tails.
    assert math.isclose(
        float(row["theoretical_type2"]),
        norm.cdf(1.65 - math.sqrt(n) * 0.5), rel_tol=2e-9,
    )
for row in read("q6b_theoretical_curve.csv"):
    n = whole(row["n"])
    close(row["power"], norm.sf(1.65 - math.sqrt(n) * 0.5))
    close(row["type2"], norm.cdf(1.65 - math.sqrt(n) * 0.5))
for row in read("validation_probabilities.csv"):
    count, total = whole(row["count"]), whole(row["total"])
    probability = float(row["theoretical_probability"])
    close(row["empirical_probability"], count / total)
    close(row["mcse_theoretical"], math.sqrt(probability * (1 - probability) / total))
    lower, upper = exact_interval(count, total)
    close(row["mc_lower95"], lower)
    close(row["mc_upper95"], upper)
checks.append("Q6B: 100000 decisions, theoretical power/type-II tails and all exact binomial intervals passed.")

# Internal cases check symbolic identities and are not assigned exam parameters.
for n, sigma2, b2 in [(1, 1, 0.25), (10, 4, 1), (100, 0.5, 4), (8, 3, 0.1)]:
    y, theta0 = 2, 2.8
    v = sigma2 / (4 * n)
    V = 1 / (1 / v + 1 / b2)
    m = V * (y / v + theta0 / b2)
    weight = b2 / (b2 + v)
    close(m, theta0 - 0.8 * weight)
    r0 = norm.cdf((theta0 - m) / math.sqrt(V))
    assert r0 > 0.5
    untruncated = norm.pdf(y, theta0, math.sqrt(b2 + v))
    m0 = quad(
        lambda theta: 2 * norm.pdf(y, theta, math.sqrt(v))
        * norm.pdf(theta, theta0, math.sqrt(b2)),
        -np.inf, theta0, epsabs=1e-12,
    )[0]
    close(m0, 2 * untruncated * r0)
    close(quad(lambda theta: 2 * norm.pdf(theta, theta0, math.sqrt(b2)), -np.inf, theta0)[0], 1)
checks.append("Q6C: posterior and truncated-normal marginal verified by quadrature on four internal cases.")
checks.append("PASSED. Numerical formulas independently verified using Python/NumPy/SciPy.")
report = "\n".join(checks) + "\n"
(RESULTS / "independent_verification.txt").write_text(report, encoding="utf-8")
print(report)
