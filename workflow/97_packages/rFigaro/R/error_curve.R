# error_curve.R
# Exponential curve fitting for expected error profiles

#' Fit an exponential curve y = a * exp(b * x) + c
#'
#' Uses nonlinear least squares (nls) to fit the expected error profile.
#'
#' @param x_values Numeric vector of x values (positions)
#' @param y_values Numeric vector of y values (expected errors)
#' @param plot_name Title for the plot (if generated)
#' @param generate_plot If TRUE, generate a ggplot2 plot object
#' @return A list with elements: a, b, c (coefficients), r_squared,
#'   predict_fn (function), plot (ggplot or NULL)
#' @keywords internal
fit_exponential_curve <- function(x_values, y_values, plot_name = "Expected error by position",
                                   generate_plot = TRUE) {
  if (length(x_values) < 3) {
    stop("Need at least 3 data points to fit an exponential curve.")
  }

  # Initial parameters matching FIGARO: p0=(0.03, 0.015, 0)
  start_a <- 0.03
  start_b <- 0.015
  start_c <- 0

  df <- data.frame(x = x_values, y = y_values)

  fit <- tryCatch({
    nls(y ~ a * exp(b * x) + c, data = df,
        start = list(a = start_a, b = start_b, c = start_c),
        algorithm = "port",
        lower = c(a = -2, b = -1, c = -8),
        upper = c(a = 2, b = 1, c = 8),
        control = nls.control(maxiter = 200, warnOnly = TRUE))
  }, error = function(e) {
    # Fallback: try with different starting values
    tryCatch({
      nls(y ~ a * exp(b * x) + c, data = df,
          start = list(a = 0.05, b = 0.02, c = -0.5),
          algorithm = "port",
          lower = c(a = -2, b = -1, c = -8),
          upper = c(a = 2, b = 1, c = 8),
          control = nls.control(maxiter = 500, warnOnly = TRUE))
    }, error = function(e2) {
      warning("Exponential curve fitting failed. Using default error model.")
      NULL
    })
  })

  if (is.null(fit)) {
    # Return a default model using the FIGARO hardcoded forward model
    predict_fn <- function(x) 0.0356 * exp(0.015 * x)
    result <- list(
      a = 0.0356, b = 0.015, c = 0,
      r_squared = NA_real_,
      predict_fn = predict_fn,
      plot = NULL
    )
    class(result) <- "figaro_curve"
    return(result)
  }

  coefs <- coef(fit)
  a_val <- coefs["a"]
  b_val <- coefs["b"]
  c_val <- coefs["c"]

  predict_fn <- function(x) a_val * exp(b_val * x) + c_val

  # Calculate R-squared
  predicted <- predict_fn(x_values)
  r_squared <- cor(y_values, predicted)^2

  plot_obj <- NULL
  if (generate_plot) {
    plot_obj <- tryCatch({
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        pred_df <- data.frame(x = x_values, observed = y_values, predicted = predicted)
        sign_c <- if (c_val >= 0) "+" else "-"
        eq_label <- sprintf("%.4fe^(%.4fx) %s %.4f\nr^2 = %.6f",
                            a_val, b_val, sign_c, abs(c_val), r_squared)

        ggplot2::ggplot(pred_df, ggplot2::aes(x = x)) +
          ggplot2::geom_line(ggplot2::aes(y = observed), colour = "black") +
          ggplot2::geom_line(ggplot2::aes(y = predicted), colour = "blue",
                             linetype = "dashed") +
          ggplot2::labs(x = "Position in Read", y = "Expected Error",
                        title = plot_name) +
          ggplot2::annotate("text", x = 0, y = max(predicted) * 0.45,
                            label = eq_label, hjust = 0, size = 3) +
          ggplot2::theme_minimal()
      } else {
        NULL
      }
    }, error = function(e) NULL)
  }

  result <- list(
    a = unname(a_val),
    b = unname(b_val),
    c = unname(c_val),
    r_squared = r_squared,
    predict_fn = predict_fn,
    plot = plot_obj
  )
  class(result) <- "figaro_curve"
  result
}

#' Calculate expected error curves for both forward and reverse reads
#'
#' @param fastq_table data.frame from find_fastq_files
#' @param subsample Subsampling factor
#' @param percentile Percentile for expected error (default 83)
#' @param generate_plots Whether to generate plots
#' @param forward_primer_length Forward primer length
#' @param reverse_primer_length Reverse primer length
#' @param n_cores Number of cores
#' @return A list with elements: forward_curve, reverse_curve (each a figaro_curve)
#' @keywords internal
calculate_ee_curves <- function(fastq_table, subsample = 1L, percentile = 83L,
                                 generate_plots = TRUE,
                                 forward_primer_length = 0L,
                                 reverse_primer_length = 0L,
                                 n_cores = 1L) {
  fwd_files <- fastq_table$file_path[fastq_table$direction == 1]
  rev_files <- fastq_table$file_path[fastq_table$direction == 2]

  # Get the sample group name for plot titles
  group_name <- fastq_table$group[1]

  # Calculate percentile EE arrays
  fwd_ee <- mean_ee_percentile_array(fwd_files, subsample = subsample,
                                      percentile = percentile,
                                      primer_length = forward_primer_length,
                                      n_cores = n_cores)
  rev_ee <- mean_ee_percentile_array(rev_files, subsample = subsample,
                                      percentile = percentile,
                                      primer_length = reverse_primer_length,
                                      n_cores = n_cores)

  ord_suffix <- function(n) {
    tens <- n %% 100
    if (tens %in% c(11, 12, 13)) return(paste0(n, "th"))
    ones <- n %% 10
    suffix <- switch(as.character(ones), "1" = "st", "2" = "nd", "3" = "rd", "th")
    paste0(n, suffix)
  }

  fwd_curve <- fit_exponential_curve(
    seq(0, length(fwd_ee) - 1), fwd_ee,
    plot_name = sprintf("%s forward reads. %s percentile", group_name, ord_suffix(percentile)),
    generate_plot = generate_plots
  )

  rev_curve <- fit_exponential_curve(
    seq(0, length(rev_ee) - 1), rev_ee,
    plot_name = sprintf("%s reverse reads. %s percentile", group_name, ord_suffix(percentile)),
    generate_plot = generate_plots
  )

  list(forward_curve = fwd_curve, reverse_curve = rev_curve)
}

#' Print method for figaro_curve objects
#' @param x A figaro_curve object
#' @param ... Additional arguments (ignored)
#' @export
print.figaro_curve <- function(x, ...) {
  sign_c <- if (x$c >= 0) "+" else "-"
  cat(sprintf("Exponential fit: %.4fe^(%.4fx) %s %.4f\n",
              x$a, x$b, sign_c, abs(x$c)))
  if (!is.na(x$r_squared)) {
    cat(sprintf("R-squared: %.6f\n", x$r_squared))
  }
  invisible(x)
}
