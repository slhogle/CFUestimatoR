#' Poisson cutoff estimator of CFU
#'
#' Pools all observations with counts below `N` (counts at or above `N` are
#' distrusted) and estimates the rate as total count divided by total dilution
#' times volume.
#'
#' Setting `N = Inf` excludes nothing and gives the naive Poisson estimate,
#' which pools every observation.
#'
#' @param counts Numeric vector of colony counts.
#' @param dil Numeric vector of dilution factors (e.g. `1e-5`), same length as
#'   `counts`.
#' @param N Counts greater than or equal to `N` are excluded. Defaults to 300,
#'   the value recommended for 100 mm petri dishes. Use `Inf` to keep all counts
#'   (naive Poisson estimate).
#' @param V Sampling volume. The default of 1 returns absolute CFU in the
#'   amount plated or spotted.
#'
#' @return A one-row data frame with columns `r`, `stderr` and `type`. `r` is
#'   `NA` if no observations fall below `N`.
#' @export
#' @examples
#' counts <- c(2, 4, 3, 25, 32)
#' dil <- c(1e-6, 1e-6, 1e-6, 1e-5, 1e-5)
#' find_poisson_cutoff(counts, dil, N = 30)
#'
#' # No cutoff: naive Poisson estimate
#' find_poisson_cutoff(counts, dil, N = Inf)
find_poisson_cutoff <- function(counts, dil, N = 300, V = 1) {
  check_counts_dil(counts, dil)
  mask <- !is.na(counts) & !is.na(dil) & counts < N
  total <- sum(counts[mask])

  if (!any(mask)) {
    r <- NA_real_
    se <- NA_real_
  } else if (total == 0) {
    r <- 0
    se <- NA_real_
  } else {
    r <- total / (sum(dil[mask]) * V)
    se <- sqrt(r * r / total)
  }
  data.frame(r = r, stderr = se, type = "poisson_cutoff")
}

#' Maximum likelihood estimator of CFU
#'
#' Maximum likelihood estimator that accounts for saturation of a spot or plate
#' at roughly `N` countable colonies (a most probable number style estimate).
#'
#' @inheritParams find_poisson_cutoff
#' @param N Estimated maximum number of countable colonies, e.g. the area of
#'   the plate or spot divided by the average colony size. Defaults to 5000, the
#'   value recommended for 100 mm petri dishes (about 50 for a 2 ul spot).
#'
#' @return A one-row data frame with columns `r`, `stderr` and `type`. If any
#'   count exceeds `N` the estimate is `Inf`. If all counts are zero the
#'   estimate is 0 with an `NA` standard error.
#' @export
#' @examples
#' find_mle(c(2, 4, 3, 25, 32), c(1e-6, 1e-6, 1e-6, 1e-5, 1e-5), N = 100)
find_mle <- function(counts, dil, N = 5000, V = 1) {
  check_counts_dil(counts, dil)
  keep <- !is.na(counts) & !is.na(dil)
  counts <- counts[keep]
  dil <- dil[keep]

  if (length(counts) == 0) {
    return(data.frame(r = NA_real_, stderr = NA_real_, type = "ML_estimator"))
  }
  if (any(counts > N)) {
    return(data.frame(r = Inf, stderr = Inf, type = "ML_estimator"))
  }
  if (all(counts == 0)) {
    return(data.frame(r = 0, stderr = NA_real_, type = "ML_estimator"))
  }

  f <- function(x) {
    sum(counts * dil / (N * -expm1(-x * dil * V / N))) - sum(dil)
  }
  r <- uniroot(f, interval = c(1e-6, 1e12), tol = 1e-10 * max(counts))$root
  p0 <- exp(-r * dil * V / N)
  inv_var <- sum(dil^2 * V^2 * counts * p0 / (N^2 * (1 - p0)^2))

  data.frame(r = r, stderr = sqrt(1 / inv_var), type = "ML_estimator")
}

#' Estimate CFU from a dilution series
#'
#' Applies either [find_mle()] (the default) or [find_poisson_cutoff()] to each
#' group of a data frame of colony counts and converts the result to CFU/ml.
#'
#' Absolute CFU in the amount plated or spotted is divided by `sampling_vol`,
#' which is always in milliliters, so `CFU_ml` is in CFU/ml. The default of 0.1
#' (100 ul) is the usual volume plated on a 100 mm petri dish (use 0.002 for a
#' 2 ul spot). The same volume must be used for all dilutions.
#'
#' @param data A data frame with one row per counted spot or plate. Input data
#'   is expected to have a header row, e.g. as read by [utils::read.csv()] or
#'   `readr::read_csv()`.
#' @param count,dilution <[`tidy-select`][tidyselect::language]> The single
#'   column holding colony counts and the single column holding dilution
#'   factors, e.g. `count = n_colonies`. Strings such as `"n_colonies"` also
#'   work. Default to columns named `count` and `dilution`. Dilutions must be
#'   numeric (e.g. `1e-6`).
#' @param group <[`tidy-select`][tidyselect::language]> Columns to group by.
#'   Any combination of tidyselect expressions can be used, e.g. `c(strain,
#'   rep)`, `starts_with("treat")`, `where(is.character)` or `all_of(vars)`.
#'   Leave as `NULL` for no grouping. Strings such as `"strain"` also work.
#' @param method Estimator to use: `"mle"` for [find_mle()] (default) or
#'   `"poisson"` for [find_poisson_cutoff()].
#' @param N Threshold passed to the chosen estimator. For `"mle"`, the maximum
#'   number of countable colonies (see [find_mle()]); defaults to 5000, as
#'   recommended for 100 mm petri dishes. For `"poisson"`, counts at or above
#'   `N` are excluded (see [find_poisson_cutoff()]); defaults to 300, also as
#'   recommended for 100 mm petri dishes.
#' @param sampling_vol Volume plated or spotted, in milliliters (ml), used to
#'   convert absolute CFU to CFU/ml. Either a single positive number applied to
#'   all rows (e.g. `0.002` for a 2 ul spot), or the name of a column holding the volume for
#'   each row (e.g. `vol_plated_ml`). Defaults to 0.1 (100 ul). The volume must
#'   be constant within each group.
#'
#' @return A data frame with one row per group: the grouping columns,
#'   `estimator` (the estimator used), `CFU_ml` (CFU/ml) and `CFU_ml_stderr`
#'   (its standard error, also per ml).
#' @export
#' @examples
#' path <- system.file("extdata", "cfu_example.tsv", package = "CFUestimatoR")
#' cfu <- read.delim(path)
#'
#' # Plate counts on 100 mm dishes; replicates are pooled within each group
#' estimate_cfu(
#'   cfu,
#'   count = colony_count,
#'   group = c(measurement, time_hours, growth_media, sample_or_strain),
#'   N = 5000,
#'   sampling_vol = vol_plated_ml
#' )
#'
#' # Poisson cutoff estimator, grouping by any tidyselect expression
#' estimate_cfu(
#'   cfu,
#'   count = colony_count,
#'   group = c(sample_or_strain, growth_media),
#'   method = "poisson",
#'   N = 500,
#'   sampling_vol = vol_plated_ml
#' )
estimate_cfu <- function(
  data,
  count = count,
  dilution = dilution,
  group = NULL,
  method = c("mle", "poisson"),
  N = NULL,
  sampling_vol = 0.1
) {
  stopifnot(is.data.frame(data))
  count <- select_one(rlang::enquo(count), data, "count")
  dilution <- select_one(rlang::enquo(dilution), data, "dilution")
  group <- names(tidyselect::eval_select(rlang::enquo(group), data))
  method <- match.arg(method)
  if (is.null(N)) {
    N <- switch(method, mle = 5000, poisson = 300)
  }
  estimator <- switch(method, mle = find_mle, poisson = find_poisson_cutoff)
  vol <- resolve_volume(rlang::enquo(sampling_vol), data)

  reserved <- intersect(group, c("estimator", "CFU_ml", "CFU_ml_stderr"))
  if (length(reserved) > 0) {
    stop(
      "`group` cannot include columns named ",
      paste0("`", reserved, "`", collapse = ", "),
      "; rename them first.",
      call. = FALSE
    )
  }

  if (length(group) == 0) {
    idx <- list(seq_len(nrow(data)))
    keys <- data.frame(row.names = 1L)
  } else {
    f <- interaction(data[group], drop = TRUE, lex.order = TRUE)
    idx <- split(seq_len(nrow(data)), f)
    keys <- data[vapply(idx, `[`, integer(1), 1), group, drop = FALSE]
  }

  est <- lapply(idx, function(i) {
    estimator(data[[count]][i], data[[dilution]][i], N = N)
  })
  est <- do.call(rbind, unname(est))

  group_vol <- vapply(
    idx,
    function(i) {
      v <- unique(vol[i])
      if (length(v) != 1) {
        stop(
          "`sampling_vol` must be constant within each group, but varies in ",
          "at least one group.",
          call. = FALSE
        )
      }
      v
    },
    numeric(1)
  )

  res <- as.data.frame(keys, stringsAsFactors = FALSE)
  rownames(res) <- NULL
  res$estimator <- est$type
  res$CFU_ml <- est$r / group_vol
  res$CFU_ml_stderr <- est$stderr / group_vol
  res
}

# Returns a numeric vector of length nrow(data): a column of `data` if the
# expression names one, otherwise a recycled scalar.
resolve_volume <- function(quo, data) {
  expr <- rlang::quo_get_expr(quo)
  is_column <- (is.symbol(expr) && as.character(expr) %in% names(data)) ||
    (is.character(expr) && length(expr) == 1 && expr %in% names(data))

  if (is_column) {
    vol <- data[[as.character(expr)]]
  } else {
    vol <- rlang::eval_tidy(quo)
    if (is.numeric(vol) && length(vol) == 1) {
      vol <- rep(vol, nrow(data))
    }
  }

  if (
    !is.numeric(vol) || length(vol) != nrow(data) || anyNA(vol) || any(vol <= 0)
  ) {
    stop(
      "`sampling_vol` must be a single positive number or the name of a ",
      "numeric column of positive volumes in ml.",
      call. = FALSE
    )
  }
  vol
}

select_one <- function(quo, data, arg) {
  loc <- tidyselect::eval_select(
    quo,
    data,
    allow_rename = FALSE,
    error_call = rlang::caller_env()
  )
  if (length(loc) != 1) {
    stop(
      "`",
      arg,
      "` must select exactly one column, not ",
      length(loc),
      ".",
      call. = FALSE
    )
  }
  names(loc)
}

check_counts_dil <- function(counts, dil) {
  if (!is.numeric(counts) || !is.numeric(dil)) {
    stop("`counts` and `dil` must be numeric.", call. = FALSE)
  }
  if (length(counts) != length(dil)) {
    stop("`counts` and `dil` must be the same length.", call. = FALSE)
  }
  if (any(counts < 0, na.rm = TRUE) || any(dil <= 0, na.rm = TRUE)) {
    stop("`counts` must be >= 0 and `dil` must be > 0.", call. = FALSE)
  }
  invisible(TRUE)
}
