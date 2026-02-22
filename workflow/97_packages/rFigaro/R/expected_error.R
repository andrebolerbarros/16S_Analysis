# expected_error.R
# Build expected error matrices from FASTQ files

#' Build expected error matrix for a single FASTQ file
#'
#' For each read, computes the cumulative expected error at every position.
#' Returns a matrix where rows = positions and columns = reads.
#'
#' @param path Path to FASTQ file
#' @param subsample Integer subsampling factor (keep every Nth read). 1 = all reads.
#' @param left_trim Number of bases to trim from the left (e.g., primer length)
#' @param start_position Start position index (0-based) for the output matrix.
#'   Cumulative EE is computed from position 0, but only positions >= start_position
#'   are returned.
#' @param compact If TRUE, store as integer (floor of EE). Saves memory.
#' @return A numeric matrix: rows = positions (from start_position), cols = reads
#' @keywords internal
build_expected_error_matrix <- function(path, subsample = 1L, left_trim = 0L,
                                        start_position = 0L, compact = TRUE) {
  # Detect encoding before opening the main read connection
  base <- detect_quality_encoding(path)

  con <- open_fastq(path)
  on.exit(close(con))
  ee_list <- list()
  read_count <- 0L
  line_count <- 0L

  repeat {
    lines <- readLines(con, n = 4L)
    if (length(lines) < 4L) break
    line_count <- line_count + 1L

    # Subsampling: keep every Nth read (1-indexed modular check)
    if ((line_count - 1L) %% subsample != 0L) next

    qual_string <- lines[4]
    if (left_trim > 0L) {
      qual_string <- substring(qual_string, left_trim + 1L)
    }

    ee <- cumulative_expected_error(qual_string, base)

    # Trim to start_position
    if (start_position > 0L && length(ee) > start_position) {
      ee <- ee[(start_position + 1L):length(ee)]
    }

    if (compact) {
      ee <- as.integer(floor(ee))
    }

    ee_list[[length(ee_list) + 1L]] <- ee
    read_count <- read_count + 1L
  }

  if (read_count == 0L) {
    stop("No reads found in: ", path)
  }

  # Ensure all vectors have the same length (pad shorter ones)
  max_len <- max(vapply(ee_list, length, integer(1)))
  ee_mat <- vapply(ee_list, function(x) {
    if (length(x) < max_len) {
      c(x, rep(if (compact) as.integer(max(x) + 1L) else max(x) + 1, max_len - length(x)))
    } else {
      x
    }
  }, if (compact) integer(max_len) else numeric(max_len))

  # ee_mat is now positions x reads
  ee_mat
}

#' Build combined expected error matrices for one direction across multiple files
#'
#' @param file_paths Character vector of FASTQ file paths
#' @param sample_ids Character vector of sample IDs (for ordering)
#' @param subsample Integer subsampling factor
#' @param start_position Start position for output
#' @param primer_length Bases to left-trim
#' @param n_cores Number of cores for parallel processing
#' @return Matrix: rows = positions, columns = reads (concatenated across files)
#' @keywords internal
build_combined_ee_matrix <- function(file_paths, sample_ids = NULL, subsample = 1L,
                                     start_position = 0L, primer_length = 0L,
                                     n_cores = 1L) {
  worker_fn <- function(fp) {
    build_expected_error_matrix(fp, subsample = subsample, left_trim = primer_length,
                                start_position = start_position, compact = TRUE)
  }

  if (n_cores > 1L && .Platform$OS.type != "windows") {
    matrices <- parallel::mclapply(file_paths, worker_fn, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    # Export required functions to cluster
    parallel::clusterExport(cl, c("build_expected_error_matrix", "open_fastq",
                                   "detect_quality_encoding", "cumulative_expected_error",
                                   "quality_to_phred", "phred_to_perror"),
                            envir = environment())
    matrices <- parallel::parLapply(cl, file_paths, worker_fn)
  } else {
    matrices <- lapply(file_paths, worker_fn)
  }

  # Ensure all matrices have the same number of rows (pad if needed)
  max_rows <- max(vapply(matrices, nrow, integer(1)))
  matrices <- lapply(matrices, function(m) {
    if (nrow(m) < max_rows) {
      pad <- matrix(as.integer(max(m) + 1L), nrow = max_rows - nrow(m), ncol = ncol(m))
      rbind(m, pad)
    } else {
      m
    }
  })

  # Column-bind: positions x (all reads)
  do.call(cbind, matrices)
}

#' Build combined expected error matrices for both forward and reverse reads
#'
#' @param fastq_table data.frame from find_fastq_files
#' @param subsample Integer subsampling factor
#' @param min_trim_positions Tuple of (forward_min, reverse_min) start positions
#' @param forward_primer_length Forward primer length to left-trim
#' @param reverse_primer_length Reverse primer length to left-trim
#' @param n_cores Number of cores
#' @return A list with elements: forward (matrix), reverse (matrix)
#' @keywords internal
build_both_ee_matrices <- function(fastq_table, subsample = 1L,
                                    min_trim_positions = c(0L, 0L),
                                    forward_primer_length = 0L,
                                    reverse_primer_length = 0L,
                                    n_cores = 1L) {
  fwd_files <- fastq_table$file_path[fastq_table$direction == 1]
  rev_files <- fastq_table$file_path[fastq_table$direction == 2]

  fwd_matrix <- build_combined_ee_matrix(
    fwd_files, subsample = subsample,
    start_position = min_trim_positions[1],
    primer_length = forward_primer_length,
    n_cores = n_cores
  )

  rev_matrix <- build_combined_ee_matrix(
    rev_files, subsample = subsample,
    start_position = min_trim_positions[2],
    primer_length = reverse_primer_length,
    n_cores = n_cores
  )

  list(forward = fwd_matrix, reverse = rev_matrix)
}

#' Compute the expected error percentile array for a single FASTQ file
#'
#' For each position, computes the Nth percentile of cumulative expected error
#' across all reads.
#'
#' @param path Path to FASTQ file
#' @param subsample Subsampling factor
#' @param percentile Percentile to compute (default 83)
#' @param primer_length Bases to left-trim
#' @return Numeric vector of percentile expected errors, one per position
#' @keywords internal
expected_error_percentile_array <- function(path, subsample = 1L, percentile = 83L,
                                             primer_length = 0L) {
  ee_mat <- build_expected_error_matrix(path, subsample = subsample,
                                         left_trim = primer_length,
                                         start_position = 0L, compact = FALSE)
  # ee_mat: rows = positions, cols = reads
  apply(ee_mat, 1, function(row) quantile(row, probs = percentile / 100, names = FALSE))
}

#' Compute mean expected error percentile arrays across multiple FASTQ files
#'
#' @param file_paths Character vector of FASTQ file paths
#' @param subsample Subsampling factor
#' @param percentile Percentile to compute
#' @param primer_length Bases to left-trim
#' @param n_cores Number of cores
#' @return Numeric vector of mean percentile expected errors
#' @keywords internal
mean_ee_percentile_array <- function(file_paths, subsample = 1L, percentile = 83L,
                                      primer_length = 0L, n_cores = 1L) {
  worker_fn <- function(fp) {
    expected_error_percentile_array(fp, subsample = subsample,
                                    percentile = percentile,
                                    primer_length = primer_length)
  }

  if (n_cores > 1L && .Platform$OS.type != "windows") {
    arrays <- parallel::mclapply(file_paths, worker_fn, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("expected_error_percentile_array",
                                   "build_expected_error_matrix", "open_fastq",
                                   "detect_quality_encoding", "cumulative_expected_error",
                                   "quality_to_phred", "phred_to_perror"),
                            envir = environment())
    arrays <- parallel::parLapply(cl, file_paths, worker_fn)
  } else {
    arrays <- lapply(file_paths, worker_fn)
  }

  # Ensure all arrays are the same length
  max_len <- max(vapply(arrays, length, integer(1)))
  arrays <- lapply(arrays, function(a) {
    if (length(a) < max_len) c(a, rep(NA_real_, max_len - length(a)))
    else a
  })

  ee_matrix <- do.call(rbind, arrays)
  colMeans(ee_matrix, na.rm = TRUE)
}
