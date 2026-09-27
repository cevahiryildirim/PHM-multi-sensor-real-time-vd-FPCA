#######################################################
## run_cv_selection.R — steps 01-03 with training-only normalization,
## then the training-only CV configuration selection (step 07).
## Usage: setwd("<this folder>"); source("run_cv_selection.R")
#######################################################
if (!dir.exists("output")) dir.create("output")
## Data location: ./data by default, or set the CMAPSS_DATA_DIR environment
## variable (02-data-preparation.R applies the same convention).
step <- function(label, file) {
  t0 <- Sys.time(); message("== ", label, " ==")
  source(file, local = FALSE)
  message(sprintf("   done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
step("01 helper functions",   "01-vd-functions.R")
step("02 data preparation",   "02-data-preparation.R")
step("03 vd-FPCA estimation", "03-vdfpca-estimation.R")
save(list = c(grep("^pc_grid_", ls(globalenv()), value = TRUE),
              grep("(traindata|testdata)$", ls(globalenv()), value = TRUE),
              "ndays_df", "RUL_TRUE"),
     file = "output/vdfpca-estimates.RData")
step("07 training-only CV selection", "07-training-cv-selection.R")
message("CV selection complete - see output/cv_selection_results.csv")
