#!/usr/bin/env Rscript
# run_figaro.R - Command-line interface for rFigaro
#
# Usage:
#   Rscript run_figaro.R -a AMPLICON_LENGTH -f FWD_PRIMER -r REV_PRIMER [options]
#
# Examples:
#   Rscript run_figaro.R -i /path/to/fastqs -a 253 -f 19 -r 20
#   Rscript run_figaro.R -i /path/to/fastqs -a 253 -f 19 -r 20 -o /path/to/output -m 20 -p 83

# Parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Defaults
  params <- list(
    input_directory = getwd(),
    output_directory = getwd(),
    amplicon_length = NULL,
    forward_primer_length = NULL,
    reverse_primer_length = NULL,
    min_overlap = 20L,
    subsample = -1L,
    percentile = 83L,
    naming_standard = "nononsense",
    output_filename = "trimParameters.json",
    n_cores = 1L
  )

  if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
    cat("
rFigaro - Optimize trimming parameters for paired-end amplicon sequencing
R implementation of FIGARO (Zymo Research)

Usage:
  Rscript run_figaro.R -a AMPLICON_LENGTH -f FWD_PRIMER -r REV_PRIMER [options]

Required arguments:
  -a, --ampliconLength        Length of amplicon (excluding primers)
  -f, --forwardPrimerLength   Length of forward primer
  -r, --reversePrimerLength   Length of reverse primer

Optional arguments:
  -i, --inputDirectory        Directory with FASTQ files (default: current dir)
  -o, --outputDirectory       Output directory (default: current dir)
  -n, --outputFileName        Output JSON filename (default: trimParameters.json)
  -m, --minimumOverlap        Minimum read overlap (default: 20)
  -s, --subsample             Subsampling factor; -1 for auto (default: -1)
  -p, --percentile            Percentile for EE model (default: 83)
  -F, --fileNamingStandard    Naming standard (default: nononsense)
                              Options: nononsense, illumina, zymo, keriksson,
                                       fvieira, yzhang
  -c, --cores                 Number of CPU cores (default: 1)
  -h, --help                  Show this help message
")
    quit(status = 0)
  }

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    val <- if (i < length(args)) args[i + 1] else NULL

    if (arg %in% c("-a", "--ampliconLength")) {
      params$amplicon_length <- as.integer(val); i <- i + 2
    } else if (arg %in% c("-f", "--forwardPrimerLength")) {
      params$forward_primer_length <- as.integer(val); i <- i + 2
    } else if (arg %in% c("-r", "--reversePrimerLength")) {
      params$reverse_primer_length <- as.integer(val); i <- i + 2
    } else if (arg %in% c("-i", "--inputDirectory")) {
      params$input_directory <- val; i <- i + 2
    } else if (arg %in% c("-o", "--outputDirectory")) {
      params$output_directory <- val; i <- i + 2
    } else if (arg %in% c("-n", "--outputFileName")) {
      params$output_filename <- val; i <- i + 2
    } else if (arg %in% c("-m", "--minimumOverlap")) {
      params$min_overlap <- as.integer(val); i <- i + 2
    } else if (arg %in% c("-s", "--subsample")) {
      params$subsample <- as.integer(val); i <- i + 2
    } else if (arg %in% c("-p", "--percentile")) {
      params$percentile <- as.integer(val); i <- i + 2
    } else if (arg %in% c("-F", "--fileNamingStandard")) {
      params$naming_standard <- val; i <- i + 2
    } else if (arg %in% c("-c", "--cores")) {
      params$n_cores <- as.integer(val); i <- i + 2
    } else {
      warning("Unknown argument: ", arg)
      i <- i + 1
    }
  }

  # Validate required
  if (is.null(params$amplicon_length)) stop("Amplicon length (-a) is required.")
  if (is.null(params$forward_primer_length)) stop("Forward primer length (-f) is required.")
  if (is.null(params$reverse_primer_length)) stop("Reverse primer length (-r) is required.")

  params
}

# Main execution
main <- function() {
  start_time <- Sys.time()

  params <- parse_args()

  # Check if rFigaro is available
  if (!requireNamespace("rFigaro", quietly = TRUE)) {
    # Try sourcing files directly
    script_dir <- dirname(sys.frame(1)$ofile)
    pkg_r_dir <- file.path(dirname(script_dir), "R")
    if (dir.exists(pkg_r_dir)) {
      for (f in list.files(pkg_r_dir, pattern = "\\.R$", full.names = TRUE)) {
        source(f)
      }
    } else {
      stop("rFigaro package not found. Install it or run from the package directory.")
    }
  } else {
    library(rFigaro)
  }

  result <- figaro(
    input_directory = params$input_directory,
    amplicon_length = params$amplicon_length,
    forward_primer_length = params$forward_primer_length,
    reverse_primer_length = params$reverse_primer_length,
    min_overlap = params$min_overlap,
    naming_standard = params$naming_standard,
    subsample = params$subsample,
    percentile = params$percentile,
    n_cores = params$n_cores,
    generate_plots = TRUE,
    output_directory = params$output_directory,
    output_filename = params$output_filename
  )

  print(result)

  elapsed <- Sys.time() - start_time
  message(sprintf("\nRun time: %s", format(elapsed)))
}

main()
