#' @title rFigaro: Optimize Trimming Parameters for Paired-End Amplicon Sequencing
#' @description R implementation of FIGARO. Analyzes paired-end FASTQ files to
#'   determine optimal trimming parameters for DADA2 and other denoising pipelines.
#' @docType package
#' @name rFigaro-package
"_PACKAGE"

#' Run FIGARO analysis to determine optimal trim parameters
#'
#' Analyzes paired-end FASTQ files from amplicon sequencing to determine the
#' optimal forward and reverse trimming positions for use with DADA2 or similar
#' denoising tools. The algorithm fits exponential expected error curves to
#' quality profiles and scores all possible trim combinations based on read
#' retention and expected error thresholds.
#'
#' @param input_directory Path to directory containing paired-end FASTQ files
#' @param amplicon_length Integer length of the target amplicon (excluding primers)
#' @param forward_primer_length Integer length of the forward primer
#' @param reverse_primer_length Integer length of the reverse primer
#' @param min_overlap Integer minimum overlap for merging paired reads (default 20)
#' @param naming_standard File naming convention. One of: "nononsense" (default),
#'   "illumina", "zymo", "keriksson", "fvieira", "yzhang"
#' @param subsample Integer subsampling factor. Use -1 for auto-calculation based
#'   on total file size (default), or a positive integer (1 = no subsampling)
#' @param percentile Integer percentile for expected error model (default 83)
#' @param n_cores Number of CPU cores for parallel processing (default 1)
#' @param generate_plots Whether to generate expected error plots (default TRUE,
#'   requires ggplot2)
#' @param output_directory Optional directory to save results (JSON + plots).
#'   If NULL (default), results are returned but not saved to disk.
#' @param output_filename Output JSON filename (default "trimParameters.json")
#' @param top_n Number of top results to return (default 20, use Inf for all)
#'
#' @return A list with elements:
#'   \describe{
#'     \item{parameters}{data.frame of ranked trim parameter candidates with
#'       columns: forward_trim, reverse_trim, forward_max_ee, reverse_max_ee,
#'       read_retention_percent, score}
#'     \item{forward_curve}{figaro_curve object for forward reads}
#'     \item{reverse_curve}{figaro_curve object for reverse reads}
#'     \item{forward_read_length}{Detected forward read length}
#'     \item{reverse_read_length}{Detected reverse read length}
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- figaro(
#'   input_directory = "path/to/fastq_files",
#'   amplicon_length = 253,
#'   forward_primer_length = 19,
#'   reverse_primer_length = 20
#' )
#'
#' # View top results
#' head(result$parameters)
#'
#' # Access the best trim positions for DADA2
#' best <- result$parameters[1, ]
#' cat("truncLen = c(", best$forward_trim, ", ", best$reverse_trim, ")\n")
#' cat("maxEE = c(", best$forward_max_ee, ", ", best$reverse_max_ee, ")\n")
#'
#' # Save results to disk
#' result <- figaro(
#'   input_directory = "path/to/fastq_files",
#'   amplicon_length = 253,
#'   forward_primer_length = 19,
#'   reverse_primer_length = 20,
#'   output_directory = "path/to/output"
#' )
#'
#' # View expected error curves
#' if (!is.null(result$forward_curve$plot)) {
#'   print(result$forward_curve$plot)
#' }
#' }
#'
#' @export
figaro <- function(input_directory,
                   amplicon_length,
                   forward_primer_length,
                   reverse_primer_length,
                   min_overlap = 20L,
                   naming_standard = "nononsense",
                   subsample = -1L,
                   percentile = 83L,
                   n_cores = 1L,
                   generate_plots = TRUE,
                   output_directory = NULL,
                   output_filename = "trimParameters.json",
                   top_n = 20L) {

  # Validate inputs
  stopifnot(is.character(input_directory), length(input_directory) == 1)
  stopifnot(dir.exists(input_directory))
  stopifnot(is.numeric(amplicon_length), amplicon_length > 0)
  stopifnot(is.numeric(forward_primer_length), forward_primer_length >= 0)
  stopifnot(is.numeric(reverse_primer_length), reverse_primer_length >= 0)
  stopifnot(is.numeric(min_overlap), min_overlap > 0)
  stopifnot(is.numeric(percentile), percentile >= 1, percentile <= 100)

  amplicon_length <- as.integer(amplicon_length)
  forward_primer_length <- as.integer(forward_primer_length)
  reverse_primer_length <- as.integer(reverse_primer_length)
  min_overlap <- as.integer(min_overlap)
  percentile <- as.integer(percentile)
  subsample <- as.integer(subsample)
  n_cores <- as.integer(max(1L, n_cores))
  top_n <- as.integer(top_n)

  # Auto-calculate subsample if -1
  if (subsample <= 0L) {
    fastq_table_temp <- find_fastq_files(input_directory, naming_standard)
    total_size <- estimate_fastq_size(fastq_table_temp)
    gigabytes <- total_size / 1e9
    subsample <- max(1L, round(gigabytes * 10))
    message("Auto-calculated subsample factor: ", subsample)
  }

  min_combined_length <- amplicon_length + min_overlap

  message("Finding FASTQ files...")
  fastq_table <- find_fastq_files(input_directory, naming_standard)
  n_fwd <- sum(fastq_table$direction == 1)
  n_rev <- sum(fastq_table$direction == 2)
  message(sprintf("Found %d forward and %d reverse FASTQ files.", n_fwd, n_rev))

  message("Checking read lengths...")
  read_lengths <- check_read_lengths(fastq_table)
  fwd_read_length <- read_lengths$forward_length
  rev_read_length <- read_lengths$reverse_length
  message(sprintf("Forward read length: %d", fwd_read_length))
  message(sprintf("Reverse read length: %d", rev_read_length))

  # Effective read lengths after primer trimming
  eff_fwd_length <- fwd_read_length - forward_primer_length
  eff_rev_length <- rev_read_length - reverse_primer_length

  if (eff_fwd_length + eff_rev_length < min_combined_length) {
    stop("Combined effective read lengths (", eff_fwd_length + eff_rev_length,
         ") are less than required minimum (", min_combined_length,
         "). Check your amplicon_length and min_overlap settings.")
  }

  message("Fitting expected error curves...")
  curves <- calculate_ee_curves(
    fastq_table, subsample = subsample, percentile = percentile,
    generate_plots = generate_plots,
    forward_primer_length = forward_primer_length,
    reverse_primer_length = reverse_primer_length,
    n_cores = n_cores
  )

  message("Building expected error matrices...")
  min_trim_pos <- calc_min_trim(eff_fwd_length, eff_rev_length, min_combined_length)

  ee_matrices <- build_both_ee_matrices(
    fastq_table, subsample = subsample,
    min_trim_positions = min_trim_pos,
    forward_primer_length = forward_primer_length,
    reverse_primer_length = reverse_primer_length,
    n_cores = n_cores
  )

  message("Generating trim position candidates...")
  trim_positions <- make_all_trim_positions(eff_fwd_length, eff_rev_length,
                                             min_combined_length)

  message("Evaluating ", nrow(trim_positions), " trim position combinations...")
  result_table <- run_trim_test(
    fwd_ee_matrix = ee_matrices$forward,
    rev_ee_matrix = ee_matrices$reverse,
    trim_positions = trim_positions,
    min_trim_positions = min_trim_pos,
    fwd_curve = curves$forward_curve,
    rev_curve = curves$reverse_curve,
    forward_primer_length = forward_primer_length,
    reverse_primer_length = reverse_primer_length
  )

  # Limit to top N
  if (is.finite(top_n) && nrow(result_table) > top_n) {
    result_table <- result_table[1:top_n, ]
  }

  message(sprintf("Analysis complete. Top score: %.2f", result_table$score[1]))
  message(sprintf("Best trim positions: forward=%d, reverse=%d",
                  result_table$forward_trim[1], result_table$reverse_trim[1]))

  result <- list(
    parameters = result_table,
    forward_curve = curves$forward_curve,
    reverse_curve = curves$reverse_curve,
    forward_read_length = fwd_read_length,
    reverse_read_length = rev_read_length
  )
  class(result) <- "figaro_result"

  # Save outputs if output_directory is specified
  if (!is.null(output_directory)) {
    save_figaro_results(result, output_directory, output_filename)
  }

  result
}

#' Run FIGARO from pre-computed error profiles
#'
#' Alternative entry point when you already have expected error data (e.g., from
#' ShortRead or dada2 quality profiles) and want to use FIGARO's scoring algorithm
#' without re-reading FASTQ files.
#'
#' @param fwd_ee_matrix Forward expected error matrix (positions x reads)
#' @param rev_ee_matrix Reverse expected error matrix (positions x reads)
#' @param forward_read_length Forward read length
#' @param reverse_read_length Reverse read length
#' @param amplicon_length Target amplicon length (excluding primers)
#' @param forward_primer_length Forward primer length
#' @param reverse_primer_length Reverse primer length
#' @param min_overlap Minimum overlap for merging (default 20)
#' @param fwd_ee_percentile Optional forward EE percentile array for curve fitting
#' @param rev_ee_percentile Optional reverse EE percentile array for curve fitting
#' @param top_n Top N results to return
#'
#' @return A figaro_result list (same structure as figaro())
#' @export
figaro_from_error_profiles <- function(fwd_ee_matrix, rev_ee_matrix,
                                        forward_read_length, reverse_read_length,
                                        amplicon_length,
                                        forward_primer_length = 0L,
                                        reverse_primer_length = 0L,
                                        min_overlap = 20L,
                                        fwd_ee_percentile = NULL,
                                        rev_ee_percentile = NULL,
                                        top_n = 20L) {

  eff_fwd <- forward_read_length - forward_primer_length
  eff_rev <- reverse_read_length - reverse_primer_length
  min_combined <- amplicon_length + min_overlap

  min_trim_pos <- calc_min_trim(eff_fwd, eff_rev, min_combined)
  trim_positions <- make_all_trim_positions(eff_fwd, eff_rev, min_combined)

  # Fit curves if percentile data provided
  fwd_curve <- NULL
  rev_curve <- NULL
  if (!is.null(fwd_ee_percentile)) {
    fwd_curve <- fit_exponential_curve(
      seq(0, length(fwd_ee_percentile) - 1), fwd_ee_percentile,
      generate_plot = FALSE
    )
  }
  if (!is.null(rev_ee_percentile)) {
    rev_curve <- fit_exponential_curve(
      seq(0, length(rev_ee_percentile) - 1), rev_ee_percentile,
      generate_plot = FALSE
    )
  }

  result_table <- run_trim_test(
    fwd_ee_matrix, rev_ee_matrix,
    trim_positions, min_trim_pos,
    fwd_curve, rev_curve,
    forward_primer_length, reverse_primer_length
  )

  if (is.finite(top_n) && nrow(result_table) > top_n) {
    result_table <- result_table[1:top_n, ]
  }

  result <- list(
    parameters = result_table,
    forward_curve = fwd_curve,
    reverse_curve = rev_curve,
    forward_read_length = forward_read_length,
    reverse_read_length = reverse_read_length
  )
  class(result) <- "figaro_result"
  result
}

#' Save FIGARO results to disk
#'
#' @param result A figaro_result object
#' @param output_directory Directory to save files
#' @param output_filename JSON output filename
#' @keywords internal
save_figaro_results <- function(result, output_directory, output_filename = "trimParameters.json") {
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  # Save JSON
  json_path <- file.path(output_directory, output_filename)
  # Build JSON-compatible list matching FIGARO's output format
  json_list <- lapply(seq_len(nrow(result$parameters)), function(i) {
    row <- result$parameters[i, ]
    list(
      trimPosition = c(row$forward_trim, row$reverse_trim),
      maxExpectedError = c(row$forward_max_ee, row$reverse_max_ee),
      readRetentionPercent = row$read_retention_percent,
      score = row$score
    )
  })
  jsonlite::write_json(json_list, json_path, pretty = TRUE, auto_unbox = TRUE)
  message("Results saved to: ", json_path)

  # Save plots if available
  if (!is.null(result$forward_curve$plot)) {
    fwd_plot_path <- file.path(output_directory, "forwardExpectedError.png")
    tryCatch({
      ggplot2::ggsave(fwd_plot_path, result$forward_curve$plot,
                      width = 8, height = 5, dpi = 150)
      message("Forward error plot saved to: ", fwd_plot_path)
    }, error = function(e) {
      warning("Could not save forward plot: ", e$message)
    })
  }

  if (!is.null(result$reverse_curve$plot)) {
    rev_plot_path <- file.path(output_directory, "reverseExpectedError.png")
    tryCatch({
      ggplot2::ggsave(rev_plot_path, result$reverse_curve$plot,
                      width = 8, height = 5, dpi = 150)
      message("Reverse error plot saved to: ", rev_plot_path)
    }, error = function(e) {
      warning("Could not save reverse plot: ", e$message)
    })
  }

  invisible(NULL)
}

#' Print method for figaro_result objects
#' @param x A figaro_result object
#' @param n Number of top results to display (default 10)
#' @param ... Additional arguments (ignored)
#' @export
print.figaro_result <- function(x, n = 10, ...) {
  cat("FIGARO Trim Parameter Optimization Results\n")
  cat("============================================\n")
  cat(sprintf("Forward read length: %d\n", x$forward_read_length))
  cat(sprintf("Reverse read length: %d\n", x$reverse_read_length))
  cat("\n")

  if (!is.null(x$forward_curve)) {
    cat("Forward expected error model:\n  ")
    print(x$forward_curve)
  }
  if (!is.null(x$reverse_curve)) {
    cat("Reverse expected error model:\n  ")
    print(x$reverse_curve)
  }
  cat("\n")

  n_show <- min(n, nrow(x$parameters))
  cat(sprintf("Top %d trim parameter candidates (of %d total):\n", n_show, nrow(x$parameters)))
  print(x$parameters[1:n_show, ], row.names = FALSE)

  cat("\n--- DADA2 recommended parameters (from top result) ---\n")
  best <- x$parameters[1, ]
  cat(sprintf("truncLen = c(%d, %d)\n", best$forward_trim, best$reverse_trim))
  cat(sprintf("maxEE = c(%d, %d)\n", best$forward_max_ee, best$reverse_max_ee))
  cat(sprintf("Expected read retention: %.1f%%\n", best$read_retention_percent))

  invisible(x)
}
