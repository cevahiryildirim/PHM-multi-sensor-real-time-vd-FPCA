#######################################################
## 06-rul-prediction.R
##
## Definitions of the nine learners and of the two input representations,
## reused by 07-training-cv-selection.R and run_final_test.R.
##
## Run on its own, this file performs the EXPLORATORY comparison: it fits a
## fixed set of configurations directly on the test engines. That is neither
## the protocol nor the table reported in the paper, which come from
## run_cv_selection.R followed by run_final_test.R.
##
## In the exploratory comparison, for every test engine j:
##   - take its feature matrix (100 train rows + row 101 = test engine),
##   - keep eligible train rows (complete cases = lives >= OBS_j),
##   - train on the capped RUL target  y_i = min(LIFE_i - OBS_j, 125)
##     (piecewise-linear RUL convention, ceiling 125),
##   - predict row 101 and convert to RUL.
## Metrics over the 100 test engines (both predicted and true RUL capped
## at 125): RMSE, MAE, and the asymmetric score of Saxena et al. (2008).
##
## Reproducibility: the per-engine seed set.seed(1000 + j) fixes the CV folds
## and the stochastic learners.
##
## Expected upstream objects (step 05): binded_all_input (K=3),
## binded_all_input_PC1_2 (K=2), RUL_TRUE. If absent, the checkpoint
## output/feature-matrices.RData is loaded.
#######################################################

suppressMessages({
  library(e1071); library(kernlab); library(glmnet)
  library(randomForest); library(xgboost); library(lightgbm); library(FNN)
  library(nnet)
})

if (!exists("binded_all_input") || !exists("binded_all_input_PC1_2")) {
  load("output/feature-matrices.RData")
}
if (!dir.exists("output")) dir.create("output")

RUL_CAP <- 125
N_SENSORS <- 9
prediction_engines <- if (exists("prediction_engines")) prediction_engines else 1:100

## ---- feature builders --------------------------------------------------
## Column layout (fixed by construction in step 05): per sensor, K PC columns
## per history step, history steps in DESCENDING domain order (step 1 = the
## test engine's current domain m = OBS_j, last step = m = 4).
pos_of <- function(p, mi, q) (p - 1) * q + mi

build_current <- function(m, q, K) {
  idx <- unlist(lapply(1:N_SENSORS, function(p)
    ((pos_of(p, 1, q) - 1) * K + 1):((pos_of(p, 1, q) - 1) * K + K)))
  m[, idx, drop = FALSE]
}

build_summaries <- function(m, q, K) {
  ## per (sensor, PC) score trajectory, oriented so index q = current domain:
  ## last value, mean, sd, slope, mean of last 5, current-minus-midpoint
  nf <- N_SENSORS * K * q
  Xall <- m[, 1:nf, drop = FALSE]
  feats <- list()
  tt <- 1:q; ttc <- tt - mean(tt); den <- sum(ttc^2)
  for (p in 1:N_SENSORS) for (k in 1:K) {
    cidx <- rev(sapply(1:q, function(mi) (pos_of(p, mi, q) - 1) * K + k))
    Z <- Xall[, cidx, drop = FALSE]
    l5 <- max(1, q - 4); hf <- max(1, ceiling(q / 2))
    feats[[length(feats) + 1]] <- cbind(
      Z[, q], rowMeans(Z), apply(Z, 1, sd),
      as.numeric(Z %*% ttc) / den,
      rowMeans(Z[, l5:q, drop = FALSE]), Z[, q] - Z[, hf])
  }
  do.call(cbind, feats)
}

standardize <- function(Xtr, Xte) {
  mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); sdv[sdv < 1e-12] <- 1
  list(tr = sweep(sweep(Xtr, 2, mu), 2, sdv, "/"),
       te = sweep(sweep(Xte, 2, mu), 2, sdv, "/"))
}

## ---- learners (each returns a scalar prediction of capped RUL) ---------
cvfolds <- function(n, k) sample(rep(1:k, length.out = n))

fit_svr <- function(X, y, Xte, n) {
  d <- ncol(X)
  grid <- expand.grid(cost = c(1, 10, 100), gamma = c(0.1, 0.5, 2) / d)
  if (n >= 10) {
    fo <- cvfolds(n, min(5, n)); best <- NULL
    for (g in seq_len(nrow(grid))) {
      errs <- c()
      for (fold in unique(fo)) {
        tr <- fo != fold
        f <- svm(x = X[tr, , drop = FALSE], y = y[tr], type = "eps-regression",
                 kernel = "radial", cost = grid$cost[g], gamma = grid$gamma[g],
                 scale = FALSE)
        errs <- c(errs, (predict(f, X[!tr, , drop = FALSE]) - y[!tr])^2)
      }
      sc <- sqrt(mean(errs))
      if (is.null(best) || sc < best$sc) best <- list(g = g, sc = sc)
    }
    gb <- grid[best$g, ]
  } else gb <- data.frame(cost = 10, gamma = 0.5 / d)
  f <- svm(x = X, y = y, type = "eps-regression", kernel = "radial",
           cost = gb$cost, gamma = gb$gamma, scale = FALSE)
  as.numeric(predict(f, Xte))
}

fit_enet <- function(X, y, Xte, n) {
  best <- NULL; nf <- max(3, min(10, n))
  for (a in c(0, 0.5, 1)) {
    cv <- tryCatch(cv.glmnet(X, y, alpha = a, nfolds = nf, grouped = FALSE,
                             standardize = FALSE), error = function(e) NULL)
    if (is.null(cv)) next
    if (is.null(best) || min(cv$cvm) < best$sc) best <- list(cv = cv, sc = min(cv$cvm))
  }
  if (is.null(best)) {
    as.numeric(predict(glmnet(X, y, alpha = 0, lambda = 1, standardize = FALSE), Xte))
  } else as.numeric(predict(best$cv, Xte, s = "lambda.min"))
}

fit_gp <- function(X, y, Xte, n) {
  f <- gausspr(x = X, y = y, kernel = "rbfdot", kpar = "automatic",
               variance.model = FALSE)
  as.numeric(predict(f, Xte))
}

fit_xgblin <- function(X, y, Xte, n) {
  ## nthread must stay 1: gblinear's parallel (shotgun) updater is
  ## non-deterministic with more threads, which would break reproducibility.
  dtr <- xgb.DMatrix(X, label = y); best <- NULL
  for (lam in c(0.1, 1, 10, 100)) {
    prm <- list(booster = "gblinear", objective = "reg:squarederror",
                lambda = lam, alpha = 0, eta = 0.5, nthread = 1)
    if (n >= 8) {
      cv <- xgb.cv(params = prm, data = dtr, nrounds = 300, nfold = min(5, n),
                   early_stopping_rounds = 20, verbose = 0)
      el <- cv$evaluation_log
      tcol <- grep("test.*rmse.*mean", names(el), value = TRUE)[1]
      bi <- which.min(el[[tcol]]); sc <- min(el[[tcol]])
    } else { bi <- 100; sc <- lam }
    if (is.null(best) || sc < best$sc) best <- list(prm = prm, ni = bi, sc = sc)
  }
  f <- xgb.train(params = best$prm, data = dtr, nrounds = best$ni, verbose = 0)
  predict(f, xgb.DMatrix(Xte))
}

fit_rf <- function(X, y, Xte, n) {
  f <- randomForest(x = X, y = y, ntree = 500)
  as.numeric(predict(f, Xte))
}

fit_xgb <- function(X, y, Xte, n) {
  dtr <- xgb.DMatrix(X, label = y)
  grid <- expand.grid(eta = c(0.05, 0.1), md = c(2, 3)); best <- NULL
  if (n >= 8) {
    for (g in seq_len(nrow(grid))) {
      prm <- list(objective = "reg:squarederror", eta = grid$eta[g],
                  max_depth = grid$md[g], subsample = 0.9,
                  colsample_bytree = 0.7, nthread = 2)
      cv <- xgb.cv(params = prm, data = dtr, nrounds = 800, nfold = min(5, n),
                   early_stopping_rounds = 30, verbose = 0)
      el <- cv$evaluation_log
      tcol <- grep("test.*rmse.*mean", names(el), value = TRUE)[1]
      bi <- which.min(el[[tcol]]); sc <- min(el[[tcol]])
      if (is.null(best) || sc < best$sc) best <- list(prm = prm, ni = bi, sc = sc)
    }
  } else {
    best <- list(prm = list(objective = "reg:squarederror", eta = 0.05,
                            max_depth = 2, subsample = 0.9,
                            colsample_bytree = 0.7, nthread = 2), ni = 300)
  }
  f <- xgb.train(params = best$prm, data = dtr, nrounds = best$ni, verbose = 0)
  predict(f, xgb.DMatrix(Xte))
}

fit_lgbm <- function(X, y, Xte, n) {
  best <- NULL
  for (nl in c(3, 7)) {
    prm <- list(objective = "regression", learning_rate = 0.05, num_leaves = nl,
                min_data_in_leaf = max(2, floor(n / 10)), feature_fraction = 0.7,
                bagging_fraction = 0.9, bagging_freq = 1,
                num_threads = 2, verbosity = -1)
    if (n >= 8) {
      cv <- lgb.cv(params = prm, data = lgb.Dataset(X, label = y), nrounds = 800,
                   nfold = min(5, n), early_stopping_rounds = 30, verbose = -1)
      bi <- cv$best_iter; sc <- cv$best_score
    } else { bi <- 200; sc <- nl }
    if (is.null(best) || sc < best$sc) best <- list(prm = prm, ni = bi, sc = sc)
  }
  f <- lgb.train(params = best$prm, data = lgb.Dataset(X, label = y),
                 nrounds = best$ni, verbose = -1)
  predict(f, Xte)
}

fit_mlp <- function(X, y, Xte, n) {
  ## single-hidden-layer neural network (nnet); size/decay chosen by CV.
  ## Deeper architectures are not meaningful here: the per-horizon training
  ## sets contain at most 100 units (as few as 4 for long-history engines).
  grid <- expand.grid(size = c(3, 5, 8), decay = c(0.01, 0.1, 1))
  if (n >= 10) {
    fo <- cvfolds(n, min(5, n)); best <- NULL
    for (g in seq_len(nrow(grid))) {
      errs <- c()
      for (fold in unique(fo)) {
        tr <- fo != fold
        f <- nnet(x = X[tr, , drop = FALSE], y = y[tr], size = grid$size[g],
                  decay = grid$decay[g], linout = TRUE, maxit = 500, trace = FALSE)
        errs <- c(errs, (predict(f, X[!tr, , drop = FALSE]) - y[!tr])^2)
      }
      sc <- sqrt(mean(errs))
      if (is.null(best) || sc < best$sc) best <- list(g = g, sc = sc)
    }
    gb <- grid[best$g, ]
  } else gb <- data.frame(size = 3, decay = 1)
  f <- nnet(x = X, y = y, size = gb$size, decay = gb$decay,
            linout = TRUE, maxit = 500, trace = FALSE)
  as.numeric(predict(f, Xte))
}

fit_knn <- function(X, y, Xte, n) {
  ks <- 1:min(20, n - 1); best <- NULL
  for (k in ks) {
    sc <- sqrt(mean((knn.reg(train = X, y = y, k = k)$pred - y)^2))  # LOOCV
    if (is.null(best) || sc < best$sc) best <- list(k = k, sc = sc)
  }
  knn.reg(train = X, test = Xte, y = y, k = best$k)$pred
}

## ---- CONFIG-BLOCK-START ------------------------------------------------
## Everything above this marker is the shared definition block (feature
## builders and fit_* functions). 07-training-cv-selection.R and
## run_final_test.R read this file and evaluate only that block, so do not
## change the marker text without updating them.
configs <- list(
  list(name = "SVR (radial), score-trajectory summaries, K=2",
       input = "K2", mode = "summ", fit = fit_svr,    std = TRUE),
  list(name = "Elastic net, current-domain scores, K=2",
       input = "K2", mode = "cur",  fit = fit_enet,   std = TRUE),
  list(name = "Gaussian process, current-domain scores, K=2",
       input = "K2", mode = "cur",  fit = fit_gp,     std = TRUE),
  list(name = "Boosted ridge, score-trajectory summaries, K=2",
       input = "K2", mode = "summ", fit = fit_xgblin, std = TRUE),
  list(name = "Random forest, current-domain scores, K=2",
       input = "K2", mode = "cur",  fit = fit_rf,     std = FALSE),
  list(name = "XGBoost, current-domain scores, K=3",
       input = "K3", mode = "cur",  fit = fit_xgb,    std = FALSE),
  list(name = "LightGBM, current-domain scores, K=2",
       input = "K2", mode = "cur",  fit = fit_lgbm,   std = FALSE),
  list(name = "k-NN, current-domain scores, K=2",
       input = "K2", mode = "cur",  fit = fit_knn,    std = TRUE),
  list(name = "Neural network (MLP), current-domain scores, K=2",
       input = "K2", mode = "cur",  fit = fit_mlp,    std = TRUE)
)

## ---- run ---------------------------------------------------------------
run_config <- function(cfg) {
  DAT <- if (cfg$input == "K2") binded_all_input_PC1_2 else binded_all_input
  K <- if (cfg$input == "K2") 2L else 3L
  pred_rul <- rep(NA_real_, 100)
  for (j in prediction_engines) {
    set.seed(1000 + j)
    m <- DAT[[j]]
    q <- (ncol(m) - 1) / (N_SENSORS * K)
    Xf <- if (cfg$mode == "cur") build_current(m, q, K) else build_summaries(m, q, K)
    colnames(Xf) <- paste0("f", seq_len(ncol(Xf)))
    y_life <- m[1:100, ncol(m)]
    ok <- complete.cases(Xf[1:100, , drop = FALSE])
    Xtr <- Xf[1:100, , drop = FALSE][ok, , drop = FALSE]
    Xte <- Xf[101, , drop = FALSE]
    if (cfg$std) { s <- standardize(Xtr, Xte); Xtr <- s$tr; Xte <- s$te }
    OBSj <- RUL_TRUE[j, "OBS"]
    ytr <- pmin(y_life[ok] - OBSj, RUL_CAP)     # capped RUL target
    pred_rul[j] <- cfg$fit(Xtr, ytr, Xte, sum(ok))
  }
  pred_rul
}

eval_metrics <- function(pred_rul, engines) {
  d <- pmin(pred_rul[engines], RUL_CAP) - pmin(RUL_TRUE[engines, "RUL"], RUL_CAP)
  c(RMSE = sqrt(mean(d^2)), MAE = mean(abs(d)),
    Score = sum(ifelse(d < 0, exp(-d / 13) - 1, exp(d / 10) - 1)))
}

results <- data.frame()
for (cfg in configs) {
  t0 <- Sys.time()
  pred_rul <- run_config(cfg)
  mets <- eval_metrics(pred_rul, prediction_engines)
  results <- rbind(results, data.frame(
    Learner = cfg$name, RMSE = round(mets["RMSE"], 2), MAE = round(mets["MAE"], 2),
    Score = round(mets["Score"]), secs = round(as.numeric(
      difftime(Sys.time(), t0, units = "secs")))))
  saveRDS(pred_rul, file.path("output", paste0(
    "predictions_", gsub("[^A-Za-z0-9]+", "_", cfg$name), ".rds")))
  cat(sprintf("%-52s RMSE=%6.2f  MAE=%6.2f  Score=%6.0f  (%ds)\n",
              cfg$name, mets["RMSE"], mets["MAE"], mets["Score"],
              round(as.numeric(difftime(Sys.time(), t0, units = "secs")))))
}
write.csv(results, "output/prediction_results.csv", row.names = FALSE)
cat("\nresults written to output/prediction_results.csv\n")
