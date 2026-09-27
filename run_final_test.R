#######################################################
## run_final_test.R — single evaluation of the CV-selected configuration
## on the official test set, after run_cv_selection.R has produced
## output/vdfpca-estimates.RData (train-only normalization) and
## output/cv_selection_results.csv (training-only CV ranking).
##
## Steps: 04 real-time scores (test engines) -> 05 feature matrices ->
## evaluate the top configuration of the CV ranking once. The other
## configurations are then evaluated as well, for reference only, and written
## to a separate file; those values play no role in the selection.
#######################################################
if (!dir.exists("output")) dir.create("output")
step <- function(label, file) {
  t0 <- Sys.time(); message("== ", label, " ==")
  source(file, local = FALSE)
  message(sprintf("   done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
load("output/vdfpca-estimates.RData")
step("04 real-time score lists", "04-realtime-scores.R")
save(list = grep("^list_scores_.*_minus$", ls(globalenv()), value = TRUE),
     file = "output/realtime-scores.RData")
step("05 feature matrices", "05-feature-matrices.R")
save(binded_all_input, binded_all_input_PC1, binded_all_input_PC1_2, RUL_TRUE,
     file = "output/feature-matrices.RData")

## ---- learners (definition block of step 06) ----------------------------
## Only the definition block of our own 06-rul-prediction.R (feature builders
## and fit_* functions) is evaluated here, so the learners stay in one place.
src06 <- readLines("06-rul-prediction.R")
cut   <- grep("^## ---- CONFIG-BLOCK-START", src06)[1]
stopifnot(!is.na(cut))
eval(parse(text = src06[seq_len(cut - 1)]), envir = globalenv())   # our own definitions
learners <- list(svr = fit_svr, enet = fit_enet, gp = fit_gp, xgblin = fit_xgblin,
                 rf = fit_rf, xgb = fit_xgb, lgbm = fit_lgbm, knn = fit_knn, mlp = fit_mlp)
std_of   <- c(svr = TRUE, enet = TRUE, gp = TRUE, xgblin = TRUE, rf = FALSE,
              xgb = FALSE, lgbm = FALSE, knn = TRUE, mlp = TRUE)
DATS <- list(binded_all_input_PC1, binded_all_input_PC1_2, binded_all_input)

run_cfg <- function(learner, mode, K) {
  DAT <- DATS[[K]]; pred <- rep(NA_real_, 100)
  for (j in 1:100) {
    set.seed(1000 + j); m <- DAT[[j]]; q <- (ncol(m) - 1) / (N_SENSORS * K)
    Xf <- if (mode == "cur") build_current(m, q, K) else build_summaries(m, q, K)
    colnames(Xf) <- paste0("f", seq_len(ncol(Xf)))
    y_life <- m[1:100, ncol(m)]; ok <- complete.cases(Xf[1:100, , drop = FALSE])
    Xtr <- Xf[1:100, , drop = FALSE][ok, , drop = FALSE]; Xte <- Xf[101, , drop = FALSE]
    if (std_of[[learner]]) { s <- standardize(Xtr, Xte); Xtr <- s$tr; Xte <- s$te }
    ytr <- pmin(y_life[ok] - RUL_TRUE[j, "OBS"], RUL_CAP)
    pred[j] <- tryCatch(learners[[learner]](Xtr, ytr, Xte, sum(ok)), error = function(e) NA_real_)
  }
  pred
}
## Both the prediction and the true RUL are capped at RUL_CAP before the
## errors are taken, so all three metrics are computed on capped values.
## n is the number of engines actually evaluated (a learner that fails on an
## engine leaves NA, which would otherwise be dropped silently).
metrics <- function(pred) {
  d <- pmin(pred, RUL_CAP) - pmin(RUL_TRUE[, "RUL"], RUL_CAP)
  c(RMSE = sqrt(mean(d^2, na.rm = TRUE)), MAE = mean(abs(d), na.rm = TRUE),
    Score = sum(ifelse(d < 0, exp(-d / 13) - 1, exp(d / 10) - 1), na.rm = TRUE),
    n = sum(!is.na(d)))
}

## ---- the ONE pre-specified evaluation ----------------------------------
cv <- read.csv("output/cv_selection_results.csv", stringsAsFactors = FALSE)
sel <- cv[1, ]
message(sprintf("CV-selected configuration: %s / %s / K=%d (CV RMSE %.2f)", sel$learner, sel$mode, sel$K, sel$RMSE))
pred_sel <- run_cfg(sel$learner, sel$mode, sel$K)
saveRDS(pred_sel, "output/predictions_selected_configuration.rds")
m_sel <- metrics(pred_sel)
if (m_sel["n"] != 100) warning("selected configuration evaluated on ", m_sel["n"], " of 100 test engines")
write.csv(data.frame(learner = sel$learner, mode = sel$mode, K = sel$K,
                     RMSE = round(m_sel["RMSE"], 2), MAE = round(m_sel["MAE"], 2),
                     Score = round(m_sel["Score"]), n = m_sel["n"]),
          "output/test_selected_configuration.csv", row.names = FALSE)
cat(sprintf("\nTEST (single evaluation) %s/%s/K=%d: RMSE=%.2f MAE=%.2f Score=%.0f (n=%d)\n",
            sel$learner, sel$mode, sel$K, m_sel["RMSE"], m_sel["MAE"], m_sel["Score"], m_sel["n"]))

## ---- reference: every configuration on the test set ---------------------
## Computed after the selection, for reference only; these values play no
## role in choosing the configuration evaluated above.
ref <- data.frame()
for (i in seq_len(nrow(cv))) {
  p <- run_cfg(cv$learner[i], cv$mode[i], cv$K[i]); mm <- metrics(p)
  ref <- rbind(ref, data.frame(learner = cv$learner[i], mode = cv$mode[i], K = cv$K[i],
                               CV_RMSE = cv$RMSE[i], test_RMSE = round(mm["RMSE"], 2),
                               test_MAE = round(mm["MAE"], 2), test_Score = round(mm["Score"]),
                               n = mm["n"]))
  cat(sprintf("  [reference] %-7s %-4s K=%d  CV %.2f  test %.2f / %.2f / %.0f\n",
              cv$learner[i], cv$mode[i], cv$K[i], cv$RMSE[i], mm["RMSE"], mm["MAE"], mm["Score"]))
}
write.csv(ref, "output/test_all_configurations_posthoc.csv", row.names = FALSE)
message("final test evaluation complete")
