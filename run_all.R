#######################################################
## run_all.R — exploratory comparison of the nine learners
##
## NOTE: this is NOT the protocol reported in the paper. It fits every
## configuration directly on the test engines and is kept only as the
## exploratory comparison. For the results in the paper use the two-step
## protocol instead:
##     source("run_cv_selection.R")   # selection on the training engines
##     source("run_final_test.R")     # one evaluation on the test engines
##
## Steps run here:
##   01 helper functions -> 02 data -> 03 vd-FPCA estimation ->
##   04 real-time score lists -> 05 feature matrices -> 06 prediction.
## Checkpoints are written to output/ after the expensive stages so a
## partial run can be resumed by sourcing the remaining steps manually.
##
## Usage:  setwd("<this folder>"); source("run_all.R")
## Data:   place the CSVs in ./data or set CMAPSS_DATA_DIR (see README).
## Runtime: roughly 1.5-2.5 h single-threaded (steps 04 and 06 dominate).
#######################################################

if (!dir.exists("output")) dir.create("output")

## source() with local = FALSE evaluates each step in the global workspace,
## so the steps see each other's objects exactly as in an interactive run.
step <- function(label, file) {
  t0 <- Sys.time()
  message("== ", label, " ==")
  source(file, local = FALSE)
  message(sprintf("   done in %.1f min",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

step("01 helper functions",      "01-vd-functions.R")
step("02 data preparation",      "02-data-preparation.R")
step("03 vd-FPCA estimation",    "03-vdfpca-estimation.R")
save(list = c(grep("^pc_grid_", ls(globalenv()), value = TRUE),
              grep("(traindata|testdata)$", ls(globalenv()), value = TRUE),
              "ndays_df", "RUL_TRUE"),
     file = "output/vdfpca-estimates.RData")

step("04 real-time score lists", "04-realtime-scores.R")
save(list = grep("^list_scores_.*_minus$", ls(globalenv()), value = TRUE),
     file = "output/realtime-scores.RData")

step("05 feature matrices",      "05-feature-matrices.R")
save(binded_all_input, binded_all_input_PC1, binded_all_input_PC1_2, RUL_TRUE,
     file = "output/feature-matrices.RData")

step("06 RUL prediction",        "06-rul-prediction.R")
message("Pipeline complete — see output/prediction_results.csv")
