# quality_scores.R
# Quality score encoding and expected error calculation utilities

#' Convert a quality character to a numeric Phred score
#' @param char Single character from a quality string
#' @param base Integer encoding base (default 33 for Sanger/Illumina 1.8+)
#' @return Integer Phred score
#' @keywords internal
char_to_phred <- function(char, base = 33L) {
  utf8ToInt(char) - base
}

#' Convert an entire quality string to a numeric Phred score vector
#' @param quality_string A quality string from a FASTQ record
#' @param base Integer encoding base (default 33)
#' @return Integer vector of Phred scores
#' @keywords internal
quality_to_phred <- function(quality_string, base = 33L) {
  utf8ToInt(quality_string) - base
}

#' Convert Phred score to error probability
#' @param phred Numeric Phred score(s)
#' @return Numeric error probability
#' @keywords internal
phred_to_perror <- function(phred) {
  10^(-phred / 10)
}

#' Compute the cumulative expected error array for a quality string
#' @param quality_string A quality string from a FASTQ record
#' @param base Integer encoding base (default 33)
#' @return Numeric vector of cumulative expected errors at each position
#' @keywords internal
cumulative_expected_error <- function(quality_string, base = 33L) {
  phred <- quality_to_phred(quality_string, base)
  perrors <- phred_to_perror(phred)
  cumsum(perrors)
}

#' Detect quality score encoding from a FASTQ file
#' @param path Path to FASTQ file
#' @param n_lines Number of quality lines to sample (default 100)
#' @return Integer base value (33 or 64)
#' @keywords internal
detect_quality_encoding <- function(path, n_lines = 100L) {
  con <- open_fastq(path)
  on.exit(close(con))

  min_ord <- 126L
  max_ord <- 33L
  count <- 0L

  while (count < n_lines) {
    # Read 4 lines (one FASTQ record)
    lines <- readLines(con, n = 4L)
    if (length(lines) < 4L) break
    qual <- lines[4L]
    ords <- utf8ToInt(qual)
    min_ord <- min(min_ord, min(ords))
    max_ord <- max(max_ord, max(ords))
    count <- count + 1L
  }

  # Sanger/Illumina 1.8+ uses base 33 (chars ! to J, ASCII 33-74)
  # Illumina 1.3-1.7 uses base 64 (chars @ to h, ASCII 64-104)
  if (min_ord < 59L) {
    return(33L)
  } else {
    return(64L)
  }
}
