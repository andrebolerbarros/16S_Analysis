# trim_prediction.R
# Core trim parameter prediction algorithm

#' Calculate minimum trim base for paired reads
#'
#' Given forward and reverse read lengths and a required minimum combined length,
#' compute the minimum acceptable trim position for each direction.
#'
#' @param forward_length Integer length of forward reads (after primer trimming)
#' @param reverse_length Integer length of reverse reads (after primer trimming)
#' @param min_combined_length Integer minimum combined read length
#'   (amplicon_length + minimum_overlap)
#' @return Integer vector c(min_forward, min_reverse)
#' @keywords internal
calc_min_trim <- function(forward_length, reverse_length, min_combined_length) {
  if (forward_length + reverse_length < min_combined_length) {
    warning("Combined read lengths (", forward_length + reverse_length,
            ") are less than required combined length (", min_combined_length, ").")
    return(c(forward_length, reverse_length))
  }
  min_forward <- min_combined_length - reverse_length
  min_reverse <- min_combined_length - forward_length
  c(min_forward, min_reverse)
}

#' Generate all possible trim position combinations
#'
#' Enumerates every possible (forward_trim, reverse_trim) pair such that the
#' combined trimmed length equals min_combined_length. Positions are 0-indexed
#' (matching the expected error matrix indices).
#'
#' @param forward_length Forward read length (after primer trimming)
#' @param reverse_length Reverse read length (after primer trimming)
#' @param min_combined_length Minimum combined length
#' @return A matrix with 2 columns (forward_pos, reverse_pos), each row a trim combo
#' @keywords internal
make_all_trim_positions <- function(forward_length, reverse_length, min_combined_length) {
  mins <- calc_min_trim(forward_length, reverse_length, min_combined_length)
  min_fwd <- mins[1]
  min_rev <- mins[2]

  # 0-indexed positions
  fwd_positions <- seq(min_fwd - 1L, forward_length - 1L)
  rev_positions <- seq(reverse_length - 1L, min_rev - 1L)

  # The trim positions trade off: as forward gets longer, reverse gets shorter
  n <- length(fwd_positions)
  if (n != length(rev_positions)) {
    # Adjust to have matching lengths
    n <- min(n, length(rev_positions))
    fwd_positions <- fwd_positions[1:n]
    rev_positions <- rev_positions[1:n]
  }

  cbind(forward_pos = fwd_positions, reverse_pos = rev_positions)
}

#' Pad a raw expected error value to a max EE threshold
#'
#' Ceiling of the value plus 1, matching FIGARO's padMaxExpectedError.
#'
#' @param raw_value Numeric expected error value
#' @return Integer padded max expected error
#' @keywords internal
pad_max_ee <- function(raw_value) {
  ceiling(raw_value) + 1L
}

#' Calculate the score for a trim parameter set
#'
#' Score = (readRetention * 100) - ((fwdMaxEE - 1)^2 + (revMaxEE - 1)^2)
#'
#' @param read_retention Fraction of reads retained (0-1)
#' @param fwd_max_ee Forward max expected error
#' @param rev_max_ee Reverse max expected error
#' @return Numeric score
#' @keywords internal
calc_score <- function(read_retention, fwd_max_ee, rev_max_ee) {
  (read_retention * 100) - ((fwd_max_ee - 1)^2 + (rev_max_ee - 1)^2)
}

#' Run trim parameter testing (lite version)
#'
#' For each candidate trim position pair, determines the max expected error
#' threshold from the fitted curves, then counts what fraction of reads would
#' pass that filter. Returns a ranked table of results.
#'
#' @param fwd_ee_matrix Forward expected error matrix (positions x reads)
#' @param rev_ee_matrix Reverse expected error matrix (positions x reads)
#' @param trim_positions Matrix from make_all_trim_positions
#' @param min_trim_positions Integer vector c(fwd_min, rev_min)
#' @param fwd_curve Forward figaro_curve object (or NULL for default model)
#' @param rev_curve Reverse figaro_curve object (or NULL for default model)
#' @param forward_primer_length Forward primer length (for adjusting output positions)
#' @param reverse_primer_length Reverse primer length (for adjusting output positions)
#' @return A data.frame with columns: forward_trim, reverse_trim, forward_max_ee,
#'   reverse_max_ee, read_retention, score, sorted by score descending
#' @keywords internal
run_trim_test <- function(fwd_ee_matrix, rev_ee_matrix,
                           trim_positions,
                           min_trim_positions = c(0L, 0L),
                           fwd_curve = NULL, rev_curve = NULL,
                           forward_primer_length = 0L,
                           reverse_primer_length = 0L) {

  fwd_min <- min_trim_positions[1]
  rev_min <- min_trim_positions[2]

  n_combos <- nrow(trim_positions)
  results <- vector("list", n_combos)

  for (i in seq_len(n_combos)) {
    fwd_pos <- trim_positions[i, 1]
    rev_pos <- trim_positions[i, 2]

    # Get max expected error from curve or default model
    if (!is.null(fwd_curve)) {
      fwd_max_ee <- pad_max_ee(fwd_curve$predict_fn(fwd_pos))
    } else {
      fwd_max_ee <- pad_max_ee(0.0356 * exp(0.015 * fwd_pos))
    }

    if (!is.null(rev_curve)) {
      rev_max_ee <- pad_max_ee(rev_curve$predict_fn(rev_pos))
    } else {
      rev_max_ee <- pad_max_ee(0.0289 * exp(0.0203 * rev_pos))
    }

    # Get expected errors for all reads at these positions
    fwd_row_idx <- fwd_pos - fwd_min + 1L
    rev_row_idx <- rev_pos - rev_min + 1L

    # Bounds checking
    if (fwd_row_idx < 1 || fwd_row_idx > nrow(fwd_ee_matrix)) next
    if (rev_row_idx < 1 || rev_row_idx > nrow(rev_ee_matrix)) next

    fwd_errors <- fwd_ee_matrix[fwd_row_idx, ]
    rev_errors <- rev_ee_matrix[rev_row_idx, ]

    # Handle different numbers of reads between forward and reverse
    n_reads <- min(length(fwd_errors), length(rev_errors))
    fwd_errors <- fwd_errors[1:n_reads]
    rev_errors <- rev_errors[1:n_reads]

    # Count kept reads: those where both forward and reverse EE < maxEE
    kept <- sum(fwd_errors < fwd_max_ee & rev_errors < rev_max_ee)
    retention <- kept / n_reads

    score <- calc_score(retention, fwd_max_ee, rev_max_ee)

    # Output positions are 1-indexed and include primer length offset
    results[[i]] <- data.frame(
      forward_trim = fwd_pos + 1L + forward_primer_length,
      reverse_trim = rev_pos + 1L + reverse_primer_length,
      forward_max_ee = fwd_max_ee,
      reverse_max_ee = rev_max_ee,
      read_retention_percent = round(100 * retention, 2),
      score = round(score, 4),
      stringsAsFactors = FALSE
    )
  }

  result_df <- do.call(rbind, results[!vapply(results, is.null, logical(1))])
  result_df <- result_df[order(-result_df$score), ]
  rownames(result_df) <- NULL
  result_df
}
