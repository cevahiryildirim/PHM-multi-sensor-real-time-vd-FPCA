# Data

These CSV files are derived from the public NASA C-MAPSS turbofan degradation
benchmark (Saxena, A., Goebel, K., Simon, D., & Eklund, N.,
2008, "Damage propagation modeling for aircraft engine run-to-failure
simulation", IEEE PHM 2008; data from the NASA Prognostics Center of
Excellence Data Set Repository).

- `<n>-<sensor>train_test_all.csv` — one file per sensor
  (T24, T30, T50, P30, Ps30, phi, BPR, W31, W32). Rows 1-100 are the
  training units (run to failure), rows 101–200 are the test units
  (truncated before failure); columns are operating cycles.
- `1-T24test_RUL.csv` — true remaining useful life of the 100 test units
  (columns: Test.eng, OBS, RUL, LIFE).

Sensor values are the raw C-MAPSS measurements arranged per-unit per-cycle;
min–max normalization is performed by `02-data-preparation.R` at load time.
