#######################################################
## 02-data-preparation.R
## Loads the per-sensor C-MAPSS CSVs (see README, "Data") and builds
## the normalized sensor matrices <s>traindata (100x362) / <s>testdata
## (100x303) for the 9 sensors, plus RUL_TRUE (Test.eng, OBS, RUL, LIFE).
#######################################################

## Directory holding the CSVs (1-T24train_test_all.csv, ..., 1-T24test_RUL.csv).
## Override by defining base_data_dir before sourcing, or via the
## CMAPSS_DATA_DIR environment variable; defaults to ./data.
if (!exists("base_data_dir")) {
  base_data_dir <- Sys.getenv("CMAPSS_DATA_DIR", unset = "data")
}

RUL_TRUE <- read.csv(
  file.path(base_data_dir, "1-T24test_RUL.csv"),
  header = TRUE
)
RUL_TRUE <- as.matrix(RUL_TRUE)

prepare_sensor_dataset <- function(file_name,
                                   sensor_label,
                                   train_rows = 1:100,
                                   test_rows = 101:200,
                                   total_cols = 2:363,
                                   test_cols = 303,
                                   transpose_test_cols = 303) {
  all_train_test <- read.csv(
    file.path(base_data_dir, file_name),
    header = TRUE,
    row.names = 1
  )

  raw_matrix <- as.matrix(all_train_test[, total_cols])
  ## Min-max constants are estimated from the TRAINING rows only and then
  ## applied unchanged to the test rows (inductive normalization).
  observed_values <- na.omit(as.vector(raw_matrix[train_rows, , drop = FALSE]))
  min_value <- min(observed_values)
  max_value <- max(observed_values)

  data_matrix <- matrix(as.numeric(raw_matrix), 200, 362)
  data_matrix <- (data_matrix - min_value) / (max_value - min_value)

  train_data <- data_matrix[train_rows, , drop = FALSE]
  test_data <- data_matrix[test_rows, 1:test_cols, drop = FALSE]

  t_train_data <- t(data_matrix[train_rows, , drop = FALSE])
  t_test_data <- t(data_matrix[test_rows, 1:transpose_test_cols, drop = FALSE])

  colnames(train_data) <- seq_len(362)
  rownames(train_data) <- seq_len(length(train_rows))
  colnames(test_data) <- seq_len(test_cols)
  rownames(test_data) <- seq_len(length(test_rows))

  list(
    all_train_test = all_train_test,
    data = data_matrix,
    min_value = min_value,
    max_value = max_value,
    traindata = train_data,
    t_traindata = t_train_data,
    testdata = test_data,
    t_testdata = t_test_data
  )
}

sensor_specs <- list(
  list(prefix = "T24", display = "T24", file = "1-T24train_test_all.csv", transpose_test_cols = 362),
  list(prefix = "T30", display = "T30", file = "2-T30train_test_all.csv", transpose_test_cols = 303),
  list(prefix = "T50", display = "T50", file = "3-T50train_test_all.csv", transpose_test_cols = 303),
  list(prefix = "P30", display = "P30", file = "4-P30train_test_all.csv", transpose_test_cols = 303),
  list(prefix = "ps30", display = "ps30", file = "5-ps30train_test_all.csv", transpose_test_cols = 303),
  list(prefix = "phi", display = "phi", file = "6-phitrain_test_all.csv", transpose_test_cols = 303),
  list(prefix = "BPR", display = "BPR", file = "7-BPRtrain_test_all.csv", transpose_test_cols = 303),
  list(prefix = "W31", display = "W31", file = "8-W31train_test_all.csv", transpose_test_cols = 303),
  list(prefix = "W32", display = "W32", file = "9-W32train_test_all.csv", transpose_test_cols = 303)
)

for (spec in sensor_specs) {
  sensor_data <- prepare_sensor_dataset(
    file_name = spec$file,
    sensor_label = spec$display,
    transpose_test_cols = spec$transpose_test_cols
  )

  assign(paste0(spec$prefix, "_all_train_test"), sensor_data$all_train_test)
  assign(paste0(spec$prefix, "data"), sensor_data$data)
  assign(paste0("min", spec$prefix), sensor_data$min_value)
  assign(paste0("max", spec$prefix), sensor_data$max_value)
  assign(paste0(spec$prefix, "traindata"), sensor_data$traindata)
  assign(paste0("t", spec$prefix, "traindata"), sensor_data$t_traindata)
  assign(paste0(spec$prefix, "testdata"), sensor_data$testdata)
  assign(paste0("t", spec$prefix, "testdata"), sensor_data$t_testdata)
}

sensor_prefixes <- vapply(sensor_specs, function(spec) spec$prefix, character(1))

train_matrix_names <- paste0(sensor_prefixes, "traindata")
test_matrix_names <- paste0(sensor_prefixes, "testdata")
t_train_matrix_names <- paste0("t", sensor_prefixes, "traindata")
t_test_matrix_names <- paste0("t", sensor_prefixes, "testdata")

sensor_train_matrices <- mget(train_matrix_names, inherits = FALSE)
sensor_test_matrices <- mget(test_matrix_names, inherits = FALSE)
sensor_train_matrices_t <- mget(t_train_matrix_names, inherits = FALSE)
sensor_test_matrices_t <- mget(t_test_matrix_names, inherits = FALSE)

names(sensor_train_matrices) <- sensor_prefixes
names(sensor_test_matrices) <- sensor_prefixes
names(sensor_train_matrices_t) <- sensor_prefixes
names(sensor_test_matrices_t) <- sensor_prefixes

stopifnot(
  all(vapply(sensor_train_matrices, is.matrix, logical(1))),
  all(vapply(sensor_test_matrices, is.matrix, logical(1))),
  all(vapply(sensor_train_matrices_t, is.matrix, logical(1))),
  all(vapply(sensor_test_matrices_t, is.matrix, logical(1)))
)
