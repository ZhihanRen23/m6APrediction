#' One-hot style DNA base encoding
#'
#' Convert equal-length DNA sequences into a data frame,
#' where each position is represented as a factor column (A/T/C/G).
#' Useful for merging with model features or one-hot expansion.
#'
#' @param dna_vec A character vector, each element is a DNA sequence
#' of equal length containing only A/T/C/G.
#' @return A data.frame with n rows and L columns (L is the sequence length),
#' column names are nt_pos1..nt_posL, and factor levels are fixed as A/T/C/G.
#' @examples
#' dna_encoding(c("ATCG", "AACC"))
#' @export
dna_encoding <- function(dna_vec) {
  stopifnot(is.character(dna_vec), length(dna_vec) > 0)
  L <- nchar(dna_vec[1])
  if (any(nchar(dna_vec) != L)) {
    stop("All DNA sequences must have the same length.")
  }
  # Split characters into a matrix -> data frame
  m <- matrix(unlist(strsplit(dna_vec, "")), ncol = L, byrow = TRUE)
  colnames(m) <- paste0("nt_pos", seq_len(L))
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  # Fix factor level order
  df[] <- lapply(df, function(x) factor(x, levels = c("A", "T", "C", "G")))
  df
}

#' Batch prediction of m6A probability and status
#'
#' Given a trained random forest model and feature data,
#' output the predicted positive-class probability and classification
#' result for each row of input.
#'
#' @param ml_fit A trained random forest model (e.g. randomForest::randomForest object)
#' @param feature_df A data.frame containing the same feature columns
#' used during training (types/levels must match).
#' @param positive_threshold Numeric value between 0 and 1; the probability
#' threshold to define a Positive sample. Default is 0.5.
#' @return A data.frame with two additional columns:
#' \itemize{
#'   \item \code{predicted_m6A_prob}: numeric, predicted positive probability
#'   \item \code{predicted_m6A_status}: character, "Positive" or "Negative"
#' }
#' @details
#' This function assumes the model is binary classification and supports
#' \code{predict(..., type = "prob")}. If DNA sequence columns are present,
#' you can use \code{\link{dna_encoding}} to generate and merge them before prediction.
#' @examples
#' \dontrun{
#'   # Load model and example data from package extdata
#'   miniDB    <- readRDS(system.file("extdata", "miniDB.rds", package = "m6APrediction"))
#'   m6A_model <- read.csv(system.file("extdata", "m6A_model.csv", package = "m6APrediction"))
#'
#'   # Use a smaller subset for quick demonstration
#'   m6A_small <- head(m6A_model, 10)
#'
#'   out <- prediction_multiple(miniDB, m6A_small, positive_threshold = 0.6)
#'   head(out)
#' }
#' @import randomForest
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  stopifnot(is.data.frame(feature_df))
  if (!is.numeric(positive_threshold) || length(positive_threshold) != 1L ||
      is.na(positive_threshold) || positive_threshold <= 0 || positive_threshold >= 1) {
    stop("positive_threshold must be a single numeric value between 0 and 1.")
  }

  # Predict probabilities (assuming binary classification)
  prob_mat <- stats::predict(ml_fit, newdata = feature_df, type = "prob")
  if (!is.matrix(prob_mat) || ncol(prob_mat) != 2) {
    stop("The model did not return a two-column probability matrix. Check the randomForest model and type='prob'.")
  }

  # Identify positive-class column (prefer column named 'Positive', otherwise use the second column)
  pos_col <- if ("Positive" %in% colnames(prob_mat)) "Positive" else colnames(prob_mat)[2]
  pos_prob <- as.numeric(prob_mat[, pos_col])

  status <- ifelse(pos_prob >= positive_threshold, "Positive", "Negative")

  out <- feature_df
  out$predicted_m6A_prob <- pos_prob
  out$predicted_m6A_status <- status
  out
}

#' Single-sample m6A prediction
#'
#' Construct a one-row feature data.frame and call
#' \code{\link{prediction_multiple}} internally to return the
#' predicted probability and classification.
#'
#' @param ml_fit A trained random forest model.
#' @param features A named list or a one-row data.frame of input features
#' with column names matching the training data.
#' @param positive_threshold Numeric threshold for Positive classification.
#' Default is 0.5.
#' @return A named list containing \code{predicted_m6A_prob} and \code{predicted_m6A_status}.
#' @examples
#' \dontrun{
#'   miniDB    <- readRDS(system.file("extdata", "miniDB.rds", package = "m6APrediction"))
#'   m6A_model <- read.csv(system.file("extdata", "m6A_model.csv", package = "m6APrediction"))
#'
#'   # Take the first row as an example
#'   one_row <- m6A_model[1, , drop = FALSE]
#'
#'   prediction_single(miniDB, one_row, positive_threshold = 0.5)
#' }
#' @export
prediction_single <- function(ml_fit, features, positive_threshold = 0.5) {
  # Ensure a one-row data.frame
  if (is.data.frame(features)) {
    if (nrow(features) != 1L) stop("features must be a one-row data.frame or a named list.")
    newdf <- features
  } else if (is.list(features)) {
    newdf <- as.data.frame(features, stringsAsFactors = FALSE, optional = TRUE)
    if (nrow(newdf) != 1L) {
      # Force to use the first row if list expansion produces multiple rows
      newdf <- newdf[1, , drop = FALSE]
    }
  } else {
    stop("features must be either a named list or a one-row data.frame.")
  }

  outdf <- prediction_multiple(ml_fit, newdf, positive_threshold = positive_threshold)
  list(
    predicted_m6A_prob   = outdf$predicted_m6A_prob[1],
    predicted_m6A_status = outdf$predicted_m6A_status[1]
  )
}
