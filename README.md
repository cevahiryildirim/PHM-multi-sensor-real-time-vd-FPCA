# Multi-Sensor vd-FPCA for Real-Time RUL Prediction

Companion code for the paper *"Real-Time Remaining Useful Life Prediction with
Multi-Sensor Variable-Domain Functional Principal Component Analysis"*
(Yildirim, Franco-Pereira, Kundu, Lillo). The pipeline estimates
variable-domain FPCA representations of C-MAPSS sensor trajectories,
builds real-time (per-cycle) score features, and predicts remaining useful
life (RUL) with a panel of supervised learners.

## Repository layout

| File | Purpose |
|------|---------|
| `run_all.R` | End-to-end driver: sources steps 01–06 in order, writes checkpoints to `output/` |
| `01-vd-functions.R` | vd-FPCA helper functions (quadrature weights, covariance-slice prediction; adapted from Johns et al., 2019, JCGS) |
| `02-data-preparation.R` | Loads the per-sensor CSVs, min–max normalizes, builds `<s>traindata` (100×362), `<s>testdata` (100×303) for the 9 sensors and `RUL_TRUE` |
| `03-vdfpca-estimation.R` | Per-sensor vd-FPCA on training data: mean surface `gam(y ~ s(time, maxT))`, covariance surface, eigenfunctions per domain length (`pc_grid_<s>`, m = 3…303) |
| `04-realtime-scores.R` | Real-time score lists: for each test engine and each history step h (h = 1 is the full observed window), start-aligned scores of the test engine and of every eligible train engine (life ≥ OBS) at domain m = OBS−h+1. Mean surfaces are fit on eligible **training** curves only; the eigenbasis is training-only throughout |
| `05-feature-matrices.R` | Assembles the supervised feature matrices `binded_all_input` (K=3), `binded_all_input_PC1_2` (K=2), `binded_all_input_PC1` (K=1): one 101-row matrix per test engine (100 train rows + the test engine as row 101), last column = LIFE |
| `06-rul-prediction.R` | The 9 learner configurations reported in the paper (SVR, elastic net, Gaussian process, boosted ridge, random forest, XGBoost, LightGBM, k-NN) with per-engine CV tuning; evaluates RMSE / MAE / Saxena score and writes `output/prediction_results.csv` |

## Requirements

R ≥ 4.5 with:

```
mgcv, reshape2, plyr, dplyr, caTools, refund,
e1071, kernlab, glmnet, randomForest, xgboost, lightgbm, FNN
```

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
source("run_all.R")
```

Approximate single-threaded runtime: step 03 ≈ 15 min, step 04 ≈ 60–90 min,
step 06 ≈ 30–40 min; everything else is fast. Checkpoints
(`output/vdfpca-estimates.RData`, `output/realtime-scores.RData`,
`output/feature-matrices.RData`) let you resume from any stage — e.g. to
re-run only the prediction step:

```r
setwd("<this folder>")
source("06-rul-prediction.R")   # loads output/feature-matrices.RData
```

## Evaluation protocol and expected results

Following standard practice on C-MAPSS, training targets use the
piecewise-linear RUL convention with a ceiling of 125 cycles, and RMSE / MAE /
Score are computed on the 100 test units at their last observed cycle
(Score = asymmetric penalty of Saxena et al., 2008; lower is better for all
three metrics). Per-engine seeds are fixed (`set.seed(1000 + j)`), so a full
run reproduces:

| Learner | Feature set | K | RMSE | MAE | Score |
|---------|-------------|---|------|-----|-------|
| SVR (radial) | score-trajectory summaries | 2 | **15.97** | **12.35** | **406** |
| Elastic net | current-domain scores | 2 | 17.03 | 13.15 | 623 |
| Gaussian process | current-domain scores | 2 | 17.11 | 13.27 | 518 |
| Boosted ridge | score-trajectory summaries | 2 | 17.94 | 13.95 | 792 |
| Random forest | current-domain scores | 2 | 17.98 | 13.77 | 590 |
| XGBoost | current-domain scores | 3 | 18.01 | 13.91 | 790 |
| LightGBM | current-domain scores | 2 | 18.72 | 13.99 | 999 |
| k-NN | current-domain scores | 2 | 18.90 | 14.03 | 726 |
| Neural network (MLP) | current-domain scores | 2 | 19.55 | 14.76 | 1389 |

*Feature sets:* "current-domain scores" are the 9×K vd-FPCA scores computed
from the unit's full observed window; "score-trajectory summaries" condense
each (sensor, component) score history into six statistics (current value,
mean, SD, slope, mean of last five, current-minus-midpoint).

## Citation

If you use this code, please cite the paper (reference to be added upon
publication) and Johns, Crainiceanu, Zipunnikov & Gellar (2019), *Variable-
Domain Functional Principal Component Analysis*, JCGS 28(4):993–1006.

## License

Creative Commons Attribution 4.0 International (CC BY 4.0) — see `LICENSE`.
The raw benchmark data belong to the public NASA C-MAPSS repository
(Saxena & Goebel, NASA Prognostics Center of Excellence).
