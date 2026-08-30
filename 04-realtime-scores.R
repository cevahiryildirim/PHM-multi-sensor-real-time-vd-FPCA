#######################################################
## 04-realtime-scores.R
##
## Real-time vd-FPCA score lists for the 9 sensors.
## For each test engine i and each history step h
## (h = 1 -> full observed window; larger h -> h-1 cycles removed from the
## END of the test engine's window), the working domain is
##     m = OBS_i - h + 1.
## At each (i, h) the scores of every eligible train engine
## (life >= OBS_i) and of the test engine itself (slot 101) are computed
## from their FIRST m cycles (start-aligned), using the train-estimated
## eigenfunctions at domain m.
##
## Output (names/structure identical to the original script):
##   list_scores_<sensor>_minus[[h]][[testeng]][[trainno]]  = c(PC1, PC2, PC3)
##   Empty/NULL slots mean "not available" (h beyond the test window, or
##   ineligible train engine) — downstream 05-feature-matrices.R checks
##   length()==0.
##
## FIXES vs the original exploratory script (results differ slightly):
##  (1) The mean surface is fit on ELIGIBLE TRAIN curves only; the test
##      engine's curve no longer enters mean estimation (matches the paper's
##      "no test data enter the estimation step").
##  (2) The mean used at domain m is mu(t; maxT = m) — the vd-FPCA
##      conditional mean — instead of mu(t; maxT = OBS_i) truncated to 1:m.
##  (3) The mean GAM is fit ONCE per (sensor, test engine) instead of once
##      per (sensor, test engine, h): ~q_i times fewer GAM fits (days -> ~1 h).
##
## Index convention (VERIFIED): pc_grid_<s>$efunctions[[k]] holds the
## eigenfunctions for domain length m = k + 1 (a (k+1) x npc matrix), so the
## entry for domain m is efunctions[[m - 1]] with exactly m rows.
##
## Expected upstream objects (steps 1-3):
##   <s>traindata (100x362), <s>testdata (100x303) for the 9 sensors,
##   pc_grid_<s>, ndays_df (columns: id, maxT, ...)
#######################################################

library(mgcv)
library(reshape2)
library(dplyr)

compute_minus_scores <- function(traindata, testdata, pc_grid, ndays_df,
                                 engines = seq_len(nrow(testdata)),
                                 max_hist = 300, n_pc = 3, verbose = TRUE) {
  list_scores <- vector("list", max_hist)
  for (h in seq_len(max_hist)) list_scores[[h]] <- vector("list", nrow(testdata))

  for (i in engines) {
    obs_i <- length(na.omit(testdata[i, ]))
    q_i <- obs_i - 3
    if (q_i < 1) next

    elig_ids <- as.numeric(ndays_df[ndays_df[, "maxT"] >= obs_i, "id"])
    if (!length(elig_ids)) next

    ## ---- mean surface: eligible TRAIN curves only, fit ONCE per test engine
    train_long <- reshape2::melt(traindata[elig_ids, , drop = FALSE], na.rm = TRUE) %>%
      dplyr::rename(id = Var1, time = Var2, y = value) %>%
      dplyr::filter(!is.na(y)) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(maxT = dplyr::n()) %>%
      dplyr::ungroup()
    fit_mean <- mgcv::gam(y ~ s(time, maxT), data = train_long, method = "REML")

    for (h in seq_len(min(max_hist, q_i))) {
      m <- obs_i - h + 1
      ef <- pc_grid[["efunctions"]][[m - 1]]      # (m x npc), see convention above
      ef <- ef[, seq_len(n_pc), drop = FALSE]
      mu_m <- as.vector(mgcv::predict.gam(
        fit_mean, newdata = data.frame(time = seq_len(m), maxT = m)))

      slot <- list()
      for (id in elig_ids) {
        resid <- as.numeric(traindata[id, seq_len(m)]) - mu_m
        slot[[id]] <- as.numeric(crossprod(ef, resid))
      }
      resid <- as.numeric(testdata[i, seq_len(m)]) - mu_m
      slot[[101]] <- as.numeric(crossprod(ef, resid))
      list_scores[[h]][[i]] <- slot
    }
    if (verbose) message("test engine ", i, " done (OBS=", obs_i,
                         ", eligible=", length(elig_ids), ")")
  }
  list_scores
}

## ---- driver: run for the 9 sensors (station order as in the manuscript) ----
## Set `minus_sensors` before sourcing to run a subset;
## set `minus_skip_driver <- TRUE` to only load the function.
if (!exists("minus_skip_driver") || !isTRUE(minus_skip_driver)) {
  if (!exists("minus_sensors")) {
    minus_sensors <- c("T24", "T30", "T50", "P30", "ps30", "phi", "BPR", "W31", "W32")
  }
  for (s in minus_sensors) {
    assign(paste0("list_scores_", s, "_minus"),
           compute_minus_scores(get(paste0(s, "traindata")),
                                get(paste0(s, "testdata")),
                                get(paste0("pc_grid_", s)),
                                ndays_df))
  }
}
