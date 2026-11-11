# Mathematical Statistics Final Assignment

Reproducibility materials for a seven-question mathematical statistics assignment. The repository contains the calculation scripts, saved simulation records, numerical results, English figures, and the final report.

The report covers normal parameter estimation, optimal linear estimation in location and scale families, Rao–Blackwell improvement, a reparameterized Cramér–Rao bound, Bayesian estimation, hypothesis testing, and ordered submission-time prediction for 146 students.

## Report and figures

- [English report in Markdown](output/Mathematical_Statistics_Final_Assignment_EN.md)
- [Ten-page English Word report](output/Mathematical_Statistics_Final_Assignment_EN_Compact.docx)
- [English figures](output/figures_en)

The Word report includes all five figures and the complete set of 146 ordered predictions. Its page count was checked using Microsoft Word on A4 paper.

## Repository contents

| Location | Contents |
|---|---|
| `analysis/q1_q3.R` | Calculations and simulations for Questions 1–3 |
| `analysis/q4_q6.R` | Bayesian calculations and classification/testing simulations for Questions 4–6 |
| `analysis/q7.R` | Declared cohort scenarios, order-statistic predictive medians, and model-based risk comparisons |
| `analysis/verify_q4_q6.py` | Independent checks of posterior calculations, simulation summaries, exact intervals, and Bayesian testing identities |
| `analysis/verify_q7.py` | Checks of the historical transcription, predictive median equations, complete prediction table, and hypothetical risks |
| `analysis/plot_figures_en.py` | Entry point for regenerating all five English figures |
| `analysis/plot_figures.py` and `analysis/plot_figures_final.py` | Shared plotting functions and figure-signature checks used by the English entry point |
| `analysis/run_webr.mjs` | Alternative runner for the R scripts using webR |
| `output/q1-q3` | Initial samples, replication results, parameter tables, random seeds, R session information, and an RDS simulation record |
| `output/q4-q6` | Bernoulli sample path, posterior tables, classification counts, test statistics/decisions, random seeds, and an RDS simulation record |
| `output/q7` | Historical fragments, modelling assumptions, complete predictions, scenario reference values, risk comparisons, seeds, and verification records |
| `output/figures_en` | English PNG, PDF, and SVG figures, with provenance and validation records |
| `output/figures_final/figure_manifest.json` | Reference drawing signatures used to validate the English figures |
| `examination.pdf` | Source questions, also used by the Question 7 transcription check |
| `requirements.txt` | Python dependency versions used for verification and plotting |

## Data scope

The simulation outputs are generated from the stated models and fixed random seeds. The `.rds` files retain the simulation records saved by the R scripts and can be read using `readRDS()`. CSV files provide the corresponding inspectable numerical outputs and summaries.

For Question 7, the historical information comprises 37 values displayed in the assignment and five explicitly recorded omissions. The previous cohort size and the numbers of omitted observations are unknown. The three cohort scenarios, their weights, and their distribution parameters are explicit modelling assumptions. The supplied 146 predictions are model-based forecasts. Their actual forecasting error requires the observed submission sequence.

## Software requirements

- R with base packages and `stats`. The recorded calculations used R 4.6.0 through webR 0.6.0. Each result directory includes `sessionInfo.txt`.
- Python 3.12 with the versions listed in `requirements.txt`.
- Cambria and DengXian fonts for the established figure typography. The plotting code checks that these fonts are installed.

Run all commands from the repository root.

## Reproduce the R calculations

```sh
Rscript analysis/q1_q3.R
Rscript analysis/q4_q6.R
Rscript analysis/q7.R
```

The scripts save results in their corresponding `output` directories. Random seeds and replication counts are recorded in `seeds.csv`. The master seed is 20261108.

### Alternative R runner using webR

With Node.js and npm installed, prepare the runner dependency and execute each batch.

```sh
npm install --prefix tmp/r-runtime --ignore-scripts --no-audit --no-fund webr@0.6.0
node analysis/run_webr.mjs q1-q3
node analysis/run_webr.mjs q4-q6
node analysis/run_webr.mjs q7
```

The runtime directory is created locally by npm. The runner includes the Windows file-URL adaptation required by webR 0.6.0.

## Verify results and regenerate English figures

```sh
python -m pip install -r requirements.txt
python analysis/verify_q4_q6.py
python analysis/verify_q7.py
python analysis/plot_figures_en.py
```

The verification scripts write their check summaries to the corresponding result directories. The Question 7 verifier also reads the historical excerpt from `examination.pdf` and checks the complete table in `output/q7/第七题正文.md`. The English plotting entry point verifies the saved numerical inputs and the reference drawing signatures before exporting all five figures.

The English report currently contains a repository-link insertion point in Appendix B. The repository URL can be inserted once this repository has been created and uploaded.
