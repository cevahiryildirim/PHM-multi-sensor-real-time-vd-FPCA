#######################################################
## 05-feature-matrices.R
##
## Builds the supervised feature matrices from the real-time score lists
## of step 5. For each test engine i (101-row matrix: 100 train engines +
## the test engine as row 101):
##   columns = 9 sensors x K PCs x q_i history steps, in DESCENDING domain
##   order within each sensor block (first block = current domain m = OBS_i,
##   last = m = 4); the final column is LIFE (train lives; row 101 holds the
##   test engine's true life, used for evaluation only — never as input).
## Ineligible train engines (life < OBS_i) appear as all-NA rows; drop them
## with complete.cases()/na.omit to obtain the paper's eligible set I_j.
##
## Outputs:
## - binded_all_input        (K = 3)
## - binded_all_input_PC1_2  (K = 2)
## - binded_all_input_PC1    (K = 1)
## (plus the per-sensor lists list_all_test_engs_history_<sensor><suffix>)
##
## Model comparison lives in 06-rul-prediction.R; this script only builds
## the inputs.
##
## Expected upstream objects:
## - list_scores_<sensor>_minus  (step 04)
## - T24traindata                (step 02; for train lives)
## - RUL_TRUE                    (step 02)
#######################################################

history_sensors <- c("T24", "T30", "T50", "P30", "ps30", "phi", "BPR", "W31", "W32")
history_num_test_engines <- 100
history_num_engine_rows <- 101
history_num_history_steps <- 250
history_test_row <- history_num_engine_rows

.require_workspace_object <- function(object_name, target_env = parent.frame()) {
  if (!exists(object_name, envir = target_env, inherits = TRUE)) {
    stop(
      sprintf("Missing required object '%s'. Please source the upstream script first.", object_name),
      call. = FALSE
    )
  }
}

.score_list_name <- function(sensor) {
  paste0("list_scores_", sensor, "_minus")
}

.history_output_name <- function(sensor, suffix = "") {
  paste0("list_all_test_engs_history_", sensor, suffix)
}

.build_history_matrix_list <- function(sensor,
                                       component_count = 3,
                                       target_env = parent.frame()) {
  score_list_object <- .score_list_name(sensor)
  .require_workspace_object(score_list_object, target_env = target_env)

  sensor_score_list <- get(score_list_object, envir = target_env, inherits = TRUE)
  history_list <- vector("list", history_num_test_engines)

  for (testeng in seq_len(history_num_test_engines)) {
    history_matrix <- matrix(NA_real_, history_num_engine_rows, 3)

    for (history in seq_len(history_num_history_steps)) {
      current_history <- sensor_score_list[[history]][[testeng]]
      if (length(current_history) == 0) next   # h beyond the test window

      block_matrix <- matrix(NA_real_, history_num_engine_rows, component_count)
      for (trainno in seq_len(history_num_engine_rows)) {
        current_scores <- current_history[[trainno]]
        if (length(current_scores) >= component_count) {
          block_matrix[trainno, ] <- current_scores[seq_len(component_count)]
        }                                       # else: ineligible engine -> NA
      }
      history_matrix <- cbind(history_matrix, block_matrix)
    }
    history_list[[testeng]] <- history_matrix[, -c(1:3), drop = FALSE]
  }

  history_list
}

.assign_history_outputs <- function(component_count = 3,
                                    suffix = "",
                                    target_env = parent.frame()) {
  for (sensor in history_sensors) {
    assign(
      .history_output_name(sensor, suffix = suffix),
      .build_history_matrix_list(
        sensor = sensor,
        component_count = component_count,
        target_env = target_env
      ),
      envir = target_env
    )
  }
}

.combine_history_inputs <- function(suffix = "",
                                    response_column,
                                    output_name,
                                    target_env = parent.frame()) {
  combined_inputs <- vector("list", history_num_test_engines)

  for (i in seq_len(history_num_test_engines)) {
    sensor_blocks <- lapply(
      history_sensors,
      function(sensor) {
        get(.history_output_name(sensor, suffix = suffix), envir = target_env, inherits = TRUE)[[i]]
      }
    )
    combined_inputs[[i]] <- do.call(
      cbind,
      c(sensor_blocks, list(response_column[[i]]))
    )
  }

  assign(output_name, combined_inputs, envir = target_env)
  invisible(combined_inputs)
}

.require_workspace_object("RUL_TRUE")
.require_workspace_object("T24traindata")

TRAIN_LIFES <- matrix(NA_real_, history_num_test_engines, 1)
for (i in seq_len(history_num_test_engines)) {
  TRAIN_LIFES[i, ] <- length(stats::na.omit(T24traindata[i, ]))
}

Response_column <- vector("list", history_num_test_engines)
for (i in seq_len(history_num_test_engines)) {
  Response_column[[i]] <- rbind(TRAIN_LIFES, RUL_TRUE[i, 4])
}

.assign_history_outputs(component_count = 3, suffix = "")
.assign_history_outputs(component_count = 1, suffix = "_PC1")
.assign_history_outputs(component_count = 2, suffix = "_PC1_2")

.combine_history_inputs(
  suffix = "",
  response_column = Response_column,
  output_name = "binded_all_input"
)
.combine_history_inputs(
  suffix = "_PC1",
  response_column = Response_column,
  output_name = "binded_all_input_PC1"
)
.combine_history_inputs(
  suffix = "_PC1_2",
  response_column = Response_column,
  output_name = "binded_all_input_PC1_2"
)

## Optional export: set `history_save_path <- "<file>.RData"` before sourcing.
if (exists("history_save_path")) {
  save(binded_all_input, binded_all_input_PC1, binded_all_input_PC1_2,
       RUL_TRUE, file = history_save_path)
  message("binded inputs saved to ", history_save_path)
}
