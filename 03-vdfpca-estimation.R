#######################################################
## 03-vdfpca-estimation.R
## Variable-domain FPCA workflow for sensor trajectories
##
## Expected upstream objects:
## - Sensor matrices from 02-data-preparation.R
## - Helper functions from 01-vd-functions.R
#######################################################

lowest <- 3
highest <- 303
gridM <- seq(lowest, highest, by = 1)

sensor_order <- c("T30", "T24", "T50", "P30", "ps30", "phi", "BPR", "W31", "W32")

build_sensor_long_data <- function(sensor_matrix) {
  sensor_long <- reshape2::melt(sensor_matrix, na.rm = TRUE) %>%
    dplyr::rename(id = Var1, time = Var2, y = value) %>%
    dplyr::arrange(id, time) %>%
    dplyr::filter(!is.na(y))

  ndays <- plyr::ddply(sensor_long, ~id, plyr::summarize, maxT = length(y)) %>%
    dplyr::arrange(maxT) %>%
    dplyr::filter(maxT != 1)

  ndays$newid <- seq_len(nrow(ndays))

  sensor_l2 <- dplyr::inner_join(sensor_long, ndays, by = "id") %>%
    dplyr::select(-id) %>%
    dplyr::rename(id = newid) %>%
    dplyr::select(id, maxT, time, y) %>%
    dplyr::arrange(id, time)

  list(
    long_data = sensor_long,
    ndays_df = ndays,
    l2_data = sensor_l2
  )
}

fit_mean_surface <- function(sensor_l2) {
  mgcv::gam(y ~ s(time, maxT), data = sensor_l2, method = "REML")
}

fit_covariance_surface <- function(sensor_cov_data) {
  mgcv::gam(kprod ~ s(stime, ttime, maxT), data = sensor_cov_data)
}

run_sensor_vdfpca <- function(sensor_prefix, target_env = parent.frame()) {
  train_matrix_name <- paste0(sensor_prefix, "traindata")
  test_matrix_name <- paste0(sensor_prefix, "testdata")

  train_long_name <- paste0(sensor_prefix, "traindata_long")
  train_l2_name <- paste0(sensor_prefix, "traindata_l2")
  train_fit_name <- paste0("fit_", sensor_prefix, "traindata")
  train_cov_name <- paste0(sensor_prefix, "traindata_cov")
  cov_fit_name <- paste0("cov_", sensor_prefix, "traindata")
  pc_grid_name <- paste0("pc_grid_", sensor_prefix)

  test_long_name <- paste0(sensor_prefix, "testdata_long")
  test_l2_name <- paste0(sensor_prefix, "testdata_l2")
  test_fit_name <- paste0("fit_", sensor_prefix, "traindata_new")

  train_matrix <- get(train_matrix_name, envir = target_env, inherits = TRUE)
  test_matrix <- get(test_matrix_name, envir = target_env, inherits = TRUE)

  train_data <- build_sensor_long_data(train_matrix)
  assign(train_long_name, train_data$long_data, envir = target_env)
  assign("ndays_df", train_data$ndays_df, envir = target_env)

  train_l2 <- train_data$l2_data
  assign(train_l2_name, train_l2, envir = target_env)

  train_fit <- fit_mean_surface(train_l2)
  assign(train_fit_name, train_fit, envir = target_env)

  train_l2$mean <- as.vector(stats::predict(train_fit))
  train_l2$ydiff <- train_l2$y - train_l2$mean
  assign(train_l2_name, train_l2, envir = target_env)

  train_cov <- datsetup_cov(train_l2)
  assign(train_cov_name, train_cov, envir = target_env)

  cov_fit <- fit_covariance_surface(train_cov)
  assign(cov_fit_name, cov_fit, envir = target_env)

  pc_grid <- get_pcs_M(gridM, cov_fit, Hz = 1, includezero = FALSE, npcs = 3)
  assign(pc_grid_name, pc_grid, envir = target_env)

  test_data <- build_sensor_long_data(test_matrix)
  assign(test_long_name, test_data$long_data, envir = target_env)
  assign("ndays_df_test", test_data$ndays_df, envir = target_env)

  test_l2 <- test_data$l2_data
  assign(test_l2_name, test_l2, envir = target_env)

  test_fit <- fit_mean_surface(test_l2)
  assign(test_fit_name, test_fit, envir = target_env)

  invisible(
    list(
      train_long_name = train_long_name,
      train_l2_name = train_l2_name,
      train_fit_name = train_fit_name,
      train_cov_name = train_cov_name,
      cov_fit_name = cov_fit_name,
      pc_grid_name = pc_grid_name,
      test_long_name = test_long_name,
      test_l2_name = test_l2_name,
      test_fit_name = test_fit_name
    )
  )
}

for (sensor_prefix in sensor_order) {
  run_sensor_vdfpca(sensor_prefix)
}

## Plotting helpers and example PC visuals live in:
## source("visual vdFPCA.R")
