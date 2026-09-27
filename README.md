# Multi-Sensor vd-FPCA for Real-Time RUL Prediction

Companion code for the paper *"Variable-Domain Functional PCA for Real-Time
Remaining Useful Life Prediction in Multi-Sensor Systems"*
(Yildirim, Franco-Pereira, Kundu, Lillo). The pipeline estimates
variable-domain FPCA representations of C-MAPSS sensor trajectories,
builds real-time (per-cycle) score features, and predicts remaining useful
life (RUL) with a panel of supervised learners.

## Repository layout

| File | Purpose |
|------|---------|
| `run_cv_selection.R` | **Step 1 of the protocol in the paper**: steps 01–03, then the configuration selection on the training units (`07`) |
| `run_final_test.R` | **Step 2 of the protocol in the paper**: steps 04–05, then one evaluation of the selected configuration on the test units |
| `run_all.R` | Exploratory comparison only (steps 01–06, every configuration scored directly on the test units). Not the protocol reported in the paper |
| `01-vd-functions.R` | vd-FPCA helper functions (quadrature weights, covariance-slice prediction; adapted from Johns et al., 2019, JCGS) |
| `02-data-preparation.R` | Loads the per-sensor CSVs, min–max normalizes each sensor with constants estimated from the training rows only, builds `<s>traindata` (100×362), `<s>testdata` (100×303) for the 9 sensors and `RUL_TRUE` |
| `03-vdfpca-estimation.R` | Per-sensor vd-FPCA on training data: mean surface `gam(y ~ s(time, maxT))`, covariance surface, eigenfunctions per domain length (`pc_grid_<s>`, m = 2…303) |
| `04-realtime-scores.R` | Real-time score lists: for each test engine and each history step h (h = 1 is the full observed window), start-aligned scores of the test engine and of every eligible train engine (life ≥ OBS) at domain m = OBS−h+1. Mean surfaces are fit on eligible **training** curves only; the eigenbasis is training-only throughout |
| `05-feature-matrices.R` | Assembles the supervised feature matrices `binded_all_input` (K=3), `binded_all_input_PC1_2` (K=2), `binded_all_input_PC1` (K=1): one 101-row matrix per test engine (100 train rows + the test engine as row 101), last column = LIFE |
| `06-rul-prediction.R` | Definitions of the nine learners (SVR, elastic net, Gaussian process, boosted ridge, random forest, XGBoost, LightGBM, k-NN, neural network) and of the two input representations, with per-engine CV tuning. `07` and `run_final_test.R` reuse this definition block; running `06` on its own performs the exploratory comparison and writes `output/prediction_results.csv` |
| `07-training-cv-selection.R` | Configuration selection on the training units only: 200 pseudo-test cases, 9 learners x 2 input representations x K in {1,2,3}; writes `output/cv_selection_results.csv` |

## Requirements

R ≥ 4.5 with:

```
mgcv, reshape2, plyr, dplyr, caTools,
e1071, kernlab, glmnet, randomForest, xgboost, lightgbm, FNN, nnet
```

(`parallel`, used by `07-training-cv-selection.R`, ships with R.)

Results in the paper were produced with R 4.5.2, xgboost 3.2.0.1,
lightgbm 4.x, e1071 1.7-x on Windows 11 (single machine, no GPU).

## Data

The pipeline reads 10 CSV files derived from the public NASA C-MAPSS turbofan
benchmark (Saxena & Goebel, NASA Prognostics Data Repository):
`1-T24train_test_all.csv` … `9-W32train_test_all.csv` (per-sensor matrices,
rows 1–100 = training units, rows 101–200 = test units, columns = operating
cycles) and `1-T24test_RUL.csv` (true RUL of the 100 test units).

Place them in `./data`, or point the pipeline elsewhere with the
`CMAPSS_DATA_DIR` environment variable (or by defining `base_data_dir`
before sourcing).

## How to run

```r
setwd("<this folder>")
source("run_cv_selection.R")   # selection on the training units
source("run_final_test.R")     # one evaluation on the test units
```

Approximate runtime: step 03 ≈ 5 h on one core (the trivariate covariance
smoothers; peak memory ≈ 30 GB), the selection in step 07 ≈ 50 min on 7 cores,
step 04 ≈ 65 min, step 05 under a minute. `run_all.R`, the exploratory
comparison, takes a further 30–40 min for step 06. Checkpoints
(`output/vdfpca-estimates.RData`, `output/realtime-scores.RData`,
`output/feature-matrices.RData`) let you resume from any stage — e.g. to
re-run only the prediction step:

```r
setwd("<this folder>")
source("06-rul-prediction.R")   # loads output/feature-matrices.RData
```

## Evaluation protocol and expected results

Following standard practice on C-MAPSS, training targets use the
piecewise-linear RUL convention with a ceiling of 125 cycles. The same ceiling
is applied to the prediction and to the true RUL before the errors are taken,
so all three metrics are computed on capped values. RMSE / MAE / Score are
evaluated on the 100 test units at their last observed cycle (Score =
asymmetric penalty of Saxena et al., 2008; lower is better for all three).

The learner, its input representation and the number of retained components
K are selected on the **training units only** (`run_cv_selection.R` ->
`07-training-cv-selection.R`): every training unit is truncated at two
horizons drawn from the empirical distribution of observed-window / lifetime
ratios in the test set and treated as a pseudo-test unit, with the eligible
set, the centering mean, the scores and the learner rebuilt from the remaining
training units exactly as for a test unit (200 pseudo-test cases). All
9 learners x {current-domain scores, score-history summaries} x K in {1,2,3}
are scored, and the configuration with the lowest cross-validated RMSE is then
evaluated **once** on the test set (`run_final_test.R`). Per-engine seeds are
fixed (`set.seed(1000 + j)`), so a full run reproduces:

| Learner (best configuration by CV) | Input | K | CV RMSE | CV MAE | Test RMSE | Test MAE | Test Score |
|---|---|---|---|---|---|---|---|
| **Elastic net (selected)** | history summaries | 1 | **17.60** | **13.72** | **17.26** | **13.07** | **552** |
| Boosted ridge | history summaries | 2 | 17.70 | 13.72 | 17.80 | 13.79 | 785 |
| SVR (radial) | history summaries | 2 | 19.04 | 14.00 | 15.62 | 12.11 | 395 |
| Neural network (MLP) | history summaries | 1 | 19.37 | 14.56 | 20.47 | 14.95 | 1405 |
| Random forest | history summaries | 2 | 19.47 | 15.04 | 19.07 | 14.42 | 763 |
| XGBoost | current scores | 2 | 19.68 | 14.86 | 18.41 | 14.33 | 780 |
| LightGBM | history summaries | 2 | 19.70 | 14.90 | 18.54 | 14.01 | 758 |
| Gaussian process | history summaries | 2 | 20.04* | 15.05* | 17.13 | 13.09 | 500 |
| k-NN | current scores | 2 | 20.66 | 15.39 | 18.90 | 14.03 | 726 |

The test-set columns of the non-selected rows were computed after the selection,
for reference only (`output/test_all_configurations_posthoc.csv`), and play no
role in it. The full 54-configuration CV ranking, with the number of cases each
configuration was scored on, is written to `output/cv_selection_results.csv`.

\* The Gaussian-process fit failed on one pseudo-test case, so its
cross-validated values are over 199 cases rather than 200. Every other
configuration in the table was scored on all 200.

*Inputs:* "current scores" are the 9 x K vd-FPCA scores computed from the
unit's full observed window; "history summaries" condense each (sensor,
component) score history into six statistics (current value, mean, SD, slope,
mean of the last five values, current-minus-midpoint).

## Animations

The `animations/` folder contains animated illustrations of the real-time
mechanism (the supplementary videos of the paper). They are produced from
separate plotting scripts that are not part of this pipeline, and the first one
uses XGBoost rather than the selected elastic net, since it illustrates the
per-cycle updating mechanism rather than the reported result:

| File | Content |
|------|---------|
| `W32_engine_040_online_rul_XGB.gif` | Online RUL prediction for test engine 40 (sensor W32 shown): as the observation stream extends cycle by cycle, the eligible training set is re-formed and the predicted RUL trajectory is updated in place |
| `T30_test_engine_007_score_history.gif` | Evolution of the T30 score representation for a test engine: observed trajectory, score path over the growing domain, and PC1–PC2 phase plane advancing together |
| `T30_test_engine_007_PC1_distribution.gif` | The test engine's first-component score inside the historical score distribution of the eligible training units at each horizon |

## Citation

If you use this code, please cite the paper (reference to be added upon
publication) and Johns, Crainiceanu, Zipunnikov & Gellar (2019), *Variable-
Domain Functional Principal Component Analysis*, JCGS 28(4):993–1006.

## License

Creative Commons Attribution 4.0 International (CC BY 4.0) — see `LICENSE`.
The raw benchmark data belong to the public NASA C-MAPSS repository
(Saxena & Goebel, NASA Prognostics Center of Excellence).
