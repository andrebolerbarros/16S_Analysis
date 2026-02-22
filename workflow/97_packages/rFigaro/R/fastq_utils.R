# fastq_utils.R
# FASTQ file handling, naming standard parsing, and file discovery

# Expected FASTQ file extensions
.fastq_extensions <- c(".fastq", ".fq", ".fastq.gz", ".fq.gz")

#' Open a FASTQ file (gzipped or plain)
#' @param path Path to FASTQ file
#' @return A connection object
#' @keywords internal
open_fastq <- function(path) {
  if (!file.exists(path)) {
    stop("FASTQ file not found: ", path)
  }
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, open = "r")
  } else {
    con <- file(path, open = "r")
  }
  con
}

#' Parse FASTQ file name to extract sample info
#'
#' Supports naming standards: "nononsense" (default), "illumina", "zymo",
#' "keriksson", "fvieira", "yzhang"
#'
#' @param filename The file name (basename only)
#' @param naming_standard Character string naming standard to use
#' @return A list with elements: group, sample, direction (1 or 2)
#' @keywords internal
parse_fastq_name <- function(filename, naming_standard = "nononsense") {
  # Strip path, keep only file name
  basename_only <- basename(filename)

  switch(tolower(naming_standard),
    "nononsense" = .parse_nononsense(basename_only),
    "illumina"   = .parse_illumina(basename_only),
    "zymo" = , "zymoservices" = , "zymoservicesnamingstandard" = .parse_zymo(basename_only),
    "keriksson"  = .parse_keriksson(basename_only),
    "fvieira"    = .parse_fvieira(basename_only),
    "yzhang"     = .parse_yzhang(basename_only),
    stop("Unknown naming standard: ", naming_standard,
         ". Supported: nononsense, illumina, zymo, keriksson, fvieira, yzhang")
  )
}

.parse_nononsense <- function(filename) {
  base_name <- sub("\\.(fq|fastq)(\\.gz)?$", "", filename)
  m <- regmatches(base_name, regexpr("_R?([12])(_\\d\\d\\d)?$", base_name))
  if (length(m) == 0 || nchar(m) == 0) {
    stop("Could not infer read direction from filename: ", filename)
  }
  direction <- as.integer(sub(".*_R?([12])(_\\d\\d\\d)?$", "\\1", base_name))
  sample <- sub("_R?[12](_\\d\\d\\d)?$", "", base_name)
  list(group = sample, sample = sample, direction = direction)
}

.parse_illumina <- function(filename) {
  base_name <- sub("\\..*$", "", filename)
  parts <- strsplit(base_name, "_")[[1]]
  if (length(parts) < 5) {
    stop(filename, " does not appear to be a valid Illumina file name.")
  }
  n <- length(parts)
  group <- paste(parts[1:(n - 4)], collapse = "_")
  sample_num <- as.integer(sub("S", "", parts[n - 3]))
  direction <- as.integer(sub("R", "", parts[n - 1]))
  list(group = group, sample = sample_num, direction = direction)
}

.parse_zymo <- function(filename) {
  base_name <- sub("\\..*$", "", filename)
  parts <- strsplit(base_name, "_")[[1]]
  if (length(parts) != 3) {
    stop(filename, " does not appear to be a valid Zymo Services file name.")
  }
  direction <- as.integer(sub("R", "", parts[3]))
  list(group = parts[1], sample = parts[2], direction = direction)
}

.parse_keriksson <- function(filename) {
  parts <- strsplit(filename, "\\.")[[1]]
  group <- parts[1]
  sample_dir <- strsplit(parts[2], "_")[[1]]
  sample <- sample_dir[1]
  direction <- as.integer(gsub("[Rr]", "", sample_dir[2]))
  list(group = group, sample = sample, direction = direction)
}

.parse_fvieira <- function(filename) {
  base_name <- sub("\\..*$", "", filename)
  parts <- strsplit(base_name, "_")[[1]]
  if (length(parts) != 2) {
    stop(filename, " does not appear to be a valid fvieira file name.")
  }
  direction <- as.integer(gsub("[Rr]", "", parts[2]))
  list(group = "default", sample = parts[1], direction = direction)
}

.parse_yzhang <- function(filename) {
  base_name <- sub("\\..*$", "", filename)
  parts <- strsplit(base_name, "_")[[1]]
  if (length(parts) != 3) {
    stop(filename, " does not appear to be a valid yzhang file name.")
  }
  direction <- as.integer(gsub("[Rr]", "", parts[3]))
  list(group = "default", sample = parts[1], direction = direction)
}

#' Find all FASTQ files in a directory and parse their naming info
#'
#' @param directory Path to directory containing FASTQ files
#' @param naming_standard Naming standard to use for parsing
#' @return A data.frame with columns: file_path, file_name, group, sample,
#'   direction, sample_id
#' @keywords internal
find_fastq_files <- function(directory, naming_standard = "nononsense") {
  if (!dir.exists(directory)) {
    stop("Directory not found: ", directory)
  }

  all_files <- list.files(directory, full.names = FALSE)
  fastq_files <- character(0)
  for (ext in .fastq_extensions) {
    fastq_files <- c(fastq_files, all_files[grepl(paste0(gsub("\\.", "\\\\.", ext), "$"), all_files)])
  }

  if (length(fastq_files) == 0) {
    stop("No FASTQ files found in: ", directory)
  }

  # Remove duplicates (a file could match both .fastq.gz and .gz patterns)
  fastq_files <- unique(fastq_files)

  result <- do.call(rbind, lapply(fastq_files, function(fn) {
    info <- parse_fastq_name(fn, naming_standard)
    data.frame(
      file_path = file.path(directory, fn),
      file_name = fn,
      group = as.character(info$group),
      sample = as.character(info$sample),
      direction = info$direction,
      stringsAsFactors = FALSE
    )
  }))

  result$sample_id <- paste(result$group, result$sample, sep = "_")
  result
}

#' Estimate read length from a FASTQ file
#' @param path Path to FASTQ file
#' @param n_reads Number of reads to sample (default 100)
#' @return A list with elements: length (integer), variance (numeric)
#' @keywords internal
estimate_read_length <- function(path, n_reads = 100L) {
  con <- open_fastq(path)
  on.exit(close(con))

  lengths <- integer(0)
  count <- 0L

  while (count < n_reads) {
    lines <- readLines(con, n = 4L)
    if (length(lines) < 4L) break
    lengths <- c(lengths, nchar(lines[2]))
    count <- count + 1L
  }

  if (length(lengths) == 0) {
    stop("No reads found in: ", path)
  }

  mean_len <- round(mean(lengths))
  var_len <- if (length(lengths) > 1) var(lengths) else 0
  list(length = mean_len, variance = var_len)
}

#' Check that read lengths are consistent across all FASTQ files
#' @param fastq_table A data.frame from find_fastq_files
#' @return A list with elements: forward_length, reverse_length
#' @keywords internal
check_read_lengths <- function(fastq_table) {
  fwd <- fastq_table[fastq_table$direction == 1, ]
  rev <- fastq_table[fastq_table$direction == 2, ]

  if (nrow(fwd) == 0) stop("No forward (R1) FASTQ files found.")
  if (nrow(rev) == 0) stop("No reverse (R2) FASTQ files found.")
  if (nrow(fwd) != nrow(rev)) {
    warning("Different number of forward (", nrow(fwd), ") and reverse (",
            nrow(rev), ") files.")
  }

  fwd_info <- lapply(fwd$file_path, estimate_read_length)
  rev_info <- lapply(rev$file_path, estimate_read_length)

  fwd_lengths <- sapply(fwd_info, `[[`, "length")
  rev_lengths <- sapply(rev_info, `[[`, "length")
  fwd_vars <- sapply(fwd_info, `[[`, "variance")
  rev_vars <- sapply(rev_info, `[[`, "variance")

  if (length(unique(fwd_lengths)) != 1) {
    stop("Forward read files have different lengths: ",
         paste(unique(fwd_lengths), collapse = ", "))
  }
  if (length(unique(rev_lengths)) != 1) {
    stop("Reverse read files have different lengths: ",
         paste(unique(rev_lengths), collapse = ", "))
  }
  if (any(fwd_vars > 0)) {
    warning("Forward reads have variable lengths within files.")
  }
  if (any(rev_vars > 0)) {
    warning("Reverse reads have variable lengths within files.")
  }

  list(forward_length = fwd_lengths[1], reverse_length = rev_lengths[1])
}

#' Estimate total FASTQ file sizes in bytes (accounting for gzip)
#' @param fastq_table A data.frame from find_fastq_files
#' @return Total estimated uncompressed size in bytes
#' @keywords internal
estimate_fastq_size <- function(fastq_table) {
  sizes <- file.info(fastq_table$file_path)$size
  gzipped <- grepl("\\.gz$", fastq_table$file_path)
  # Estimate gzipped FASTQ as ~3.5x compressed size
  sizes[gzipped] <- sizes[gzipped] * 3.5
  sum(sizes)
}
