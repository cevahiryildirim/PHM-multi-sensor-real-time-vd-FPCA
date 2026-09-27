#######################################################
## 07-training-cv-selection.R
##
## Training-only selection of the learner / feature representation / K.
## The 100 TRAINING engines are turned into pseudo-test cases by truncating
## each of them at horizons drawn from the empirical distribution of
## (observed window / lifetime) in the official test set. For a pseudo-test
## case (i, h):
##   - the eligible set is every OTHER training engine with life >= h,
##   - the centering mean is re-fitted on that eligible set (risk-set mean),
##   - scores are computed exactly as in 04-realtime-scores.R
##     (fixed training eigenfunctions, start-aligned first-m cycles),
##   - feature matrices are assembled exactly as in 05-feature-matrices.R,
##   - every learner configuration is fitted on the eligible rows with the
##     capped target min(L - h, 125) and predicts the case.
## Configurations are ranked by cross-validated RMSE over the pseudo-test
## cases. No test engine is fitted, predicted or scored here.
##
## The one thing taken from the test set is the empirical distribution of the
## observed fraction of life, RUL_TRUE$OBS / RUL_TRUE$LIFE, used only to draw
## the truncation points so that the pseudo-test cases have the same censoring
## pattern as the evaluation set. No individual test response enters a model,
## a target or an error.
##
## Expected upstream objects (steps 01-03): <s>traindata, pc_grid_<s>,
## ndays_df, RUL_TRUE (its OBS and LIFE columns only).
## Output: output/cv_selection_results.csv, output/cv_selection_cases.csv,
##         output/cv_selection_predictions.rds
#######################################################

suppressMessages({ library(mgcv); library(reshape2); library(dplyr); library(parallel) })

cv_sensors   <- c("T24", "T30", "T50", "P30", "ps30", "phi", "BPR", "W31", "W32")
CV_TRUNC_PER_ENGINE <- 2      # truncation points per training engine
CV_MAX_HIST  <- 250           # history steps kept (as in step 05)
CV_MIN_H     <- 20            # never truncate before cycle 20
CV_MIN_RUL   <- 5             # keep at least 5 cycles of remaining life
CV_WORKERS   <- if (exists("cv_workers")) cv_workers else max(1, min(7, detectCores() - 1))
RUL_CAP      <- 125
N_SENSORS    <- 9
if (!dir.exists("output")) dir.create("output")

## ---- 1. pseudo-test cases ------------------------------------------------
train_life <- as.numeric(ndays_df[order(as.numeric(ndays_df$id)), "maxT"])
stopifnot(length(train_life) == 100)
test_frac  <- RUL_TRUE[, "OBS"] / RUL_TRUE[, "LIFE"]          # empirical truncation fractions
## the reference eigenbasis is estimated up to the longest observed test window
## (grid 2..303), so a pseudo-test cannot be truncated beyond that domain length
CV_MAX_H   <- length(get(paste0("pc_grid_", cv_sensors[1]))[["efunctions"]]) + 1
set.seed(20260920)
cases <- do.call(rbind, lapply(1:100, function(i) {
  fr <- sample(test_frac, CV_TRUNC_PER_ENGINE, replace = TRUE)
  h  <- pmin(pmax(round(fr * train_life[i]), CV_MIN_H), train_life[i] - CV_MIN_RUL, CV_MAX_H)
  data.frame(case = NA_integer_, engine = i, life = train_life[i], h = h,
             rul_true = train_life[i] - h)
}))
cases$case <- seq_len(nrow(cases))
write.csv(cases, "output/cv_selection_cases.csv", row.names = FALSE)
if (exists("cv_case_limit")) cases <- cases[seq_len(cv_case_limit), ]   # smoke-test hook
message(sprintf("%d pseudo-test cases; h range %d-%d; median frac %.2f",
                nrow(cases), min(cases$h), max(cases$h), median(cases$h / cases$life)))

## ---- 2. scores for one case (mirrors compute_minus_scores in step 04) ----
score_case <- function(cs, traindata, pc_grid, train_life, max_hist, n_pc = 3) {
  i <- cs$engine; h <- cs$h
  q <- h - 3; if (q < 1) return(NULL)
  elig <- setdiff(which(train_life >= h), i)               # risk set, engine i excluded
  if (!length(elig)) return(NULL)
  train_long <- reshape2::melt(traindata[elig, , drop = FALSE], na.rm = TRUE) %>%
    dplyr::rename(id = Var1, time = Var2, y = value) %>%
    dplyr::filter(!is.na(y)) %>% dplyr::group_by(id) %>%
    dplyr::mutate(maxT = dplyr::n()) %>% dplyr::ungroup()
  fit_mean <- mgcv::gam(y ~ s(time, maxT), data = train_long, method = "REML")
  out <- vector("list", min(max_hist, q))
  for (hs in seq_len(min(max_hist, q))) {
    m  <- h - hs + 1
    ef <- pc_grid[["efunctions"]][[m - 1]][, seq_len(n_pc), drop = FALSE]
    mu_m <- as.vector(mgcv::predict.gam(fit_mean, newdata = data.frame(time = seq_len(m), maxT = m)))
    slot <- vector("list", 101)
    for (id in elig) slot[[id]] <- as.numeric(crossprod(ef, as.numeric(traindata[id, seq_len(m)]) - mu_m))
    slot[[101]] <- as.numeric(crossprod(ef, as.numeric(traindata[i, seq_len(m)]) - mu_m))
    out[[hs]] <- slot
  }
  out
}

## ---- 3. feature matrix for one case (mirrors step 05 layout) -------------
assemble_case <- function(score_by_sensor, K, train_life, cs) {
  blocks <- lapply(score_by_sensor, function(sl) {
    do.call(cbind, lapply(sl, function(slot) {
      b <- matrix(NA_real_, 101, K)
      for (r in seq_len(101)) if (length(slot[[r]]) >= K) b[r, ] <- slot[[r]][seq_len(K)]
      b
    }))
  })
  cbind(do.call(cbind, blocks), c(train_life, cs$life))   # last column = LIFE
}

## ---- 4. learners: reuse the definitions of step 06 without running it ----
## The feature builders and fit_* functions are defined in our own
## 06-rul-prediction.R; only its definition block (everything before the
## configuration list) is evaluated here so the learners stay in one place.
binded_all_input <- NULL; binded_all_input_PC1_2 <- NULL   # prevents 06 from loading a checkpoint
src06 <- readLines("06-rul-prediction.R")
cut   <- grep("^## ---- CONFIG-BLOCK-START", src06)[1]
stopifnot(!is.na(cut))
eval(parse(text = src06[seq_len(cut - 1)]), envir = globalenv())
learners <- list(svr = fit_svr, enet = fit_enet, gp = fit_gp, xgblin = fit_xgblin,
                 rf = fit_rf, xgb = fit_xgb, lgbm = fit_lgbm, knn = fit_knn, mlp = fit_mlp)
std_of   <- c(svr = TRUE, enet = TRUE, gp = TRUE, xgblin = TRUE, rf = FALSE,
              xgb = FALSE, lgbm = FALSE, knn = TRUE, mlp = TRUE)
grid <- expand.grid(learner = names(learners), mode = c("cur", "summ"), K = 1:3,
                    stringsAsFactors = FALSE)

predict_case <- function(feat_by_K, cs) {
  res <- numeric(nrow(grid))
  for (g in seq_len(nrow(grid))) {
    K <- grid$K[g]; m <- feat_by_K[[K]]
    q <- (ncol(m) - 1) / (N_SENSORS * K)
    Xf <- if (grid$mode[g] == "cur") build_current(m, q, K) else build_summaries(m, q, K)
    colnames(Xf) <- paste0("f", seq_len(ncol(Xf)))
    y_life <- m[1:100, ncol(m)]
    ok  <- complete.cases(Xf[1:100, , drop = FALSE])
    Xtr <- Xf[1:100, , drop = FALSE][ok, , drop = FALSE]; Xte <- Xf[101, , drop = FALSE]
    if (std_of[[grid$learner[g]]]) { s <- standardize(Xtr, Xte); Xtr <- s$tr; Xte <- s$te }
    ytr <- pmin(y_life[ok] - cs$h, RUL_CAP)
    set.seed(1000 + cs$case)
    res[g] <- tryCatch(learners[[grid$learner[g]]](Xtr, ytr, Xte, sum(ok)),
                       error = function(e) NA_real_)
  }
  res
}

## ---- 5. run: scoring + learners per case, in parallel over cases ---------
worker_fun <- function(c_idx) tryCatch({
  cs <- cases[c_idx, ]
  t0 <- Sys.time()
  sc <- lapply(cv_sensors, function(s)
    score_case(cs, get(paste0(s, "traindata")), get(paste0("pc_grid_", s)), train_life, CV_MAX_HIST))
  t_score <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (any(vapply(sc, is.null, logical(1)))) return(list(case = c_idx, pred = rep(NA_real_, nrow(grid)),
                                                        t_score = t_score, t_fit = NA))
  feat_by_K <- lapply(1:3, function(K) assemble_case(sc, K, train_life, cs))
  t1 <- Sys.time()
  pred <- predict_case(feat_by_K, cs)
  list(case = c_idx, pred = pred, t_score = t_score,
       t_fit = as.numeric(difftime(Sys.time(), t1, units = "secs")))
}, error = function(e) list(case = c_idx, pred = rep(NA_real_, nrow(grid)), t_score = NA, t_fit = NA,
                            error = conditionMessage(e)))

export_objs <- c("cases", "train_life", "cv_sensors", "CV_MAX_HIST", "CV_MAX_H", "RUL_CAP", "N_SENSORS",
                 "score_case", "assemble_case", "predict_case", "grid", "learners", "std_of",
                 "pos_of", "build_current", "build_summaries", "standardize", "cvfolds",
                 paste0("fit_", names(learners)),
                 paste0(cv_sensors, "traindata"), paste0("pc_grid_", cv_sensors))
cl <- makeCluster(CV_WORKERS)
invisible(clusterEvalQ(cl, suppressMessages({
  library(mgcv); library(reshape2); library(dplyr); library(e1071); library(kernlab)
  library(glmnet); library(randomForest); library(xgboost); library(lightgbm); library(FNN); library(nnet) })))
clusterExport(cl, export_objs, envir = globalenv())
message("running ", nrow(cases), " cases x ", nrow(grid), " configurations on ", CV_WORKERS, " workers ...")
t_all <- Sys.time()
res <- parLapply(cl, seq_len(nrow(cases)), worker_fun)
stopCluster(cl)
message(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))

## ---- 6. aggregate ---------------------------------------------------------
P <- do.call(rbind, lapply(res, `[[`, "pred"))                 # cases x configs
truth <- pmin(cases$rul_true, RUL_CAP)
metrics <- t(apply(P, 2, function(p) {
  d <- pmin(p, RUL_CAP) - truth; ok <- !is.na(d)
  c(RMSE = sqrt(mean(d[ok]^2)), MAE = mean(abs(d[ok])),
    Score = sum(ifelse(d[ok] < 0, exp(-d[ok] / 13) - 1, exp(d[ok] / 10) - 1)),
    n = sum(ok))
}))
results <- cbind(grid, as.data.frame(metrics))
results <- results[order(results$RMSE), ]
results$RMSE <- round(results$RMSE, 2); results$MAE <- round(results$MAE, 2); results$Score <- round(results$Score)
write.csv(results, "output/cv_selection_results.csv", row.names = FALSE)
saveRDS(list(cases = cases, grid = grid, pred = P,
             t_score = sapply(res, `[[`, "t_score"), t_fit = sapply(res, `[[`, "t_fit")),
        "output/cv_selection_predictions.rds")
cat("\nTraining-only CV selection (", nrow(cases), " pseudo-test cases):\n", sep = "")
print(head(results, 15), row.names = FALSE)
cat("\nselected configuration:", results$learner[1], "/", results$mode[1], "/ K =", results$K[1], "\n")
