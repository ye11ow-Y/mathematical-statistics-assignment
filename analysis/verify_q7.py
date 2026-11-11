"""Independent verification of the declared Q7 predictive model and its outputs."""

import csv
import math
import re
from pathlib import Path

import numpy as np
from pypdf import PdfReader
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import betainc, betaincc
from scipy.stats import beta


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "output" / "q7"
checks = []


def read(filename):
    with (RESULTS / filename).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def close(actual, expected, tolerance=2e-9):
    assert math.isclose(float(actual), float(expected), rel_tol=tolerance, abs_tol=tolerance), (
        actual, expected
    )


history = read("historical_fragments.csv")
page = PdfReader(ROOT / "examination.pdf").pages[8].extract_text()
display = page.split("2025:", 1)[1].split("Submission window:", 1)[0]
source_fragments = re.split(r"\.\s*\.\s*\.", display)
assert len(source_fragments) == 6
expected_values = []
for fragment_index, fragment in enumerate(source_fragments, start=1):
    values = [float(value) for value in re.findall(r"-?\d+\.\d+", fragment)]
    corresponding = [row for row in history if int(row["fragment"]) == fragment_index]
    assert len(values) == len(corresponding)
    for offset, (value, row) in enumerate(zip(values, corresponding), start=1):
        close(row["hours"], value)
        assert int(row["within_fragment"]) == offset
    expected_values.extend(values)
assert len(expected_values) == len(history) == 37
assert sum(row["ellipsis_after"] == "TRUE" for row in history) == 5
assert np.all(np.diff(expected_values) <= 0)
checks.append("Historical transcription matches the original PDF: 37 values, six fragments and five omissions, preserving repeated and negative values.")

scenarios = read("scenarios.csv")
numeric_columns = [key for key in scenarios[0] if key != "scenario"]
for scenario in scenarios:
    for key in numeric_columns:
        scenario[key] = float(scenario[key])
    close(scenario["p_near"], 1 - scenario["p_late"] - scenario["p_early"])
    assert 0 < scenario["a_near"] < 1
close(sum(row["weight"] for row in scenarios), 1)


def cdf(x, scenario):
    if x < 0:
        return scenario["p_late"] * math.exp(x / scenario["tau_late"])
    if x <= 24:
        return scenario["p_late"] + scenario["p_near"] * (x / 24)**scenario["a_near"]
    return 1 - scenario["p_early"] * math.exp(-(x - 24) / scenario["tau_early"])


def quantile(p, scenario):
    if p < scenario["p_late"]:
        return scenario["tau_late"] * math.log(p / scenario["p_late"])
    if p <= 1 - scenario["p_early"]:
        return 24 * ((p - scenario["p_late"]) / scenario["p_near"])**(1 / scenario["a_near"])
    return 24 - scenario["tau_early"] * math.log((1 - p) / scenario["p_early"])


prediction = read("predictions.csv")
N = int(prediction[0]["N_students"])
provenance = {row["item"]: row["value"] for row in read("provenance.csv")}
assert N == int(provenance["current_N"])
assert len(prediction) == N
assert [int(row["rank"]) for row in prediction] == list(range(1, N + 1))
hours = np.array([float(row["hours"]) for row in prediction])
rounded = np.array([float(row["hours_rounded"]) for row in prediction])
assert np.all(np.isfinite(hours))
assert np.all(np.diff(hours) <= 0) and np.all(np.diff(rounded) <= 0)


def predictive_order_cdf(x, descending_rank):
    j = N + 1 - descending_rank
    return sum(
        scenario["weight"] * betainc(j, descending_rank, cdf(x, scenario))
        for scenario in scenarios
    )


for row in prediction:
    i = int(row["rank"])
    x = float(row["hours"])
    beta_median = beta.ppf(0.5, N + 1 - i, i)
    close(row["beta_probability_median"], beta_median)
    close(predictive_order_cdf(x, i), 0.5, tolerance=1e-7)
    for scenario in scenarios:
        close(row[scenario["scenario"] + "_order_median"], quantile(beta_median, scenario))
    assert abs(float(row["hours_rounded"]) - x) <= 0.005 + 1e-10
for i in sorted({1, 2, 3, N // 2, N - 2, N - 1, N}):
    root = brentq(lambda x: predictive_order_cdf(x, i) - 0.5, -100, 2000, xtol=1e-10)
    close(hours[i - 1], root)
checks.append(f"All {N} predictive medians satisfy the independently recomputed cohort-mixture order-statistic equations; extreme and middle roots also independently solved.")

submission = read("predictions_for_submission.csv")
values = (RESULTS / "prediction_values.txt").read_text(encoding="utf-8").strip().split(", ")
assert len(submission) == len(values) == N
for i, (row, value) in enumerate(zip(submission, values), start=1):
    assert int(row["rank"]) == i
    assert row["hours"] == value == f"{rounded[i - 1]:.2f}"
manuscript = (RESULTS / "第七题正文.md").read_text(encoding="utf-8")
assert "<!-- PREDICTION_TABLE -->" not in manuscript and "\ufffd" not in manuscript
table_values = []
for line in manuscript.splitlines():
    cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
    if len(cells) == 8 and cells[0].isdigit():
        for offset in range(0, 8, 2):
            if cells[offset]:
                table_values.append((int(cells[offset]), cells[offset + 1]))
assert table_values == [(i, f"{rounded[i - 1]:.2f}") for i in range(1, N + 1)]
summary = read("prediction_summary.csv")[0]
close(summary["earliest_hours"], rounded[0])
close(summary["latest_hours"], rounded[-1])
assert int(summary["positive_count_rounded"]) == sum(rounded > 0)
assert int(summary["zero_count_rounded"]) == sum(rounded == 0)
assert int(summary["negative_count_rounded"]) == sum(rounded < 0)
close(summary["rounding_W1"], np.mean(abs(hours - rounded)))
checks.append("CSV, numeric-only submission sequence and complete manuscript table agree exactly; order, signs, counts and rounding passed.")

# Evaluate expected absolute loss directly from conditional order-statistic CDFs.
# This checks the hypothetical R simulation without sharing its random generator.
risk_rows = []
for scenario in scenarios:
    losses = []
    for i, value in enumerate(hours, start=1):
        j = N + 1 - i
        left_points = [point for point in [0, 24] if -100 < point < value]
        right_points = [point for point in [0, 24] if value < point < 2000]
        left = quad(
            lambda x: betainc(j, i, cdf(x, scenario)),
            -100, value, points=left_points, epsabs=1e-7,
        )[0]
        right = quad(
            lambda x: betaincc(j, i, cdf(x, scenario)),
            value, 2000, points=right_points, epsabs=1e-7,
        )[0]
        losses.append(left + right)
    exact_risk = float(np.mean(losses))
    r_result = next(
        row for row in read("hypothetical_risk_comparison.csv")
        if row["generating_scenario"] == scenario["scenario"]
        and row["candidate"] == "scenario_predictive_median"
    )
    simulated_risk, mcse = float(r_result["hypothetical_expected_W1"]), float(r_result["mcse"])
    assert abs(simulated_risk - exact_risk) <= 5 * mcse + 1e-6
    risk_rows.append({
        "scenario": scenario["scenario"],
        "hypothetical_expected_W1_by_quadrature": exact_risk,
        "R_simulated_hypothetical_W1": simulated_risk,
        "R_simulation_mcse": mcse,
    })
with (RESULTS / "independent_hypothetical_risk.csv").open("w", encoding="utf-8", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(risk_rows[0]))
    writer.writeheader()
    writer.writerows(risk_rows)
checks.append("Hypothetical W1 expectations independently integrated over all order-statistic CDFs; R simulations agree within Monte Carlo uncertainty in every scenario.")
checks.append("PASSED. This verifies model calculations and delivery integrity, not real forecast accuracy.")
report = "\n".join(checks) + "\n"
(RESULTS / "independent_verification.txt").write_text(report, encoding="utf-8")
print(report)
