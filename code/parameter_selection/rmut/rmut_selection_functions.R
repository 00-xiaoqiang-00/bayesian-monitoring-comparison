# RMUT calibration-cache post-processing functions.
# Statistical operations follow the formal selector.
# The caller supplies cache_path(case, kind); no fitting functions are loaded.
options(stringsAsFactors = FALSE)
designs <- list(
  IMpower = list(n_max = 680L, interim = c(272L, 408L, 544L)),
  KEYNOTE598 = list(n_max = 568L, interim = c(228L, 340L, 454L)),
  OAK = list(n_max = 850L, interim = c(340L, 510L, 680L)),
  survival_colon = list(n_max = 620L, interim = c(248L, 372L, 496L))
)

pp_grid <- round(seq(0.80, 0.99, by = 0.0005), 4)
lambda_grid <- round(seq(0, 5, by = 0.005), 4)

prepare_cache <- function(path) {
  x <- readRDS(path)
  attempted <- length(x)
  if (attempted != 10000L) stop("Expected 10,000 records: ", path)
  valid <- vapply(x, function(z) {
    flag <- if (!is.null(z$fit_ok)) isTRUE(z$fit_ok) else if (!is.null(z$ok)) isTRUE(z$ok) else TRUE
    pp <- if (!is.null(z$PP)) z$PP else z$PP_interims
    flag && length(pp) == 3L && all(is.finite(pp)) &&
      length(z$PP_final) == 1L && is.finite(z$PP_final)
  }, logical(1L))
  if (mean(valid) < 0.95) stop("Valid-fit rate below 95%: ", path)
  pp <- matrix(NA_real_, attempted, 3L)
  final <- rep(NA_real_, attempted)
  for (i in which(valid)) {
    z <- x[[i]]
    pp[i, ] <- as.numeric(if (!is.null(z$PP)) z$PP else z$PP_interims)
    final[i] <- as.numeric(z$PP_final)
  }
  list(
    pp = pp,
    final = final,
    valid = valid,
    attempted = attempted,
    failed = attempted - sum(valid)
  )
}

scenario_metrics <- function(cache, pp_threshold, lambda_prime, interim_counts, n_max) {
  n_used <- n_used_vector(cache, pp_threshold, lambda_prime, interim_counts, n_max)
  valid_idx <- which(cache$valid)
  boundary <- pp_threshold * (interim_counts / n_max)^lambda_prime
  below <- sweep(cache$pp[valid_idx, , drop = FALSE], 2L, boundary, FUN = "<=")
  stopped <- rowSums(below) > 0L
  success <- !stopped & cache$final[valid_idx] > pp_threshold
  success_n <- sum(success)
  c(
    success_n = success_n,
    rate = success_n / cache$attempted,
    pet = sum(stopped) / cache$attempted,
    ess = sum(n_used) / cache$attempted
  )
}

find_exact_target_pp <- function(cache, lambda_prime, interim_counts, n_max, target_success = 500L) {
  memo <- new.env(parent = emptyenv())
  count_at <- function(i) {
    key <- as.character(i)
    if (!exists(key, envir = memo, inherits = FALSE)) {
      value <- scenario_metrics(cache, pp_grid[i], lambda_prime, interim_counts, n_max)[["success_n"]]
      assign(key, as.integer(value), envir = memo)
    }
    get(key, envir = memo, inherits = FALSE)
  }
  n <- length(pp_grid)
  if (count_at(1L) < target_success || count_at(n) > target_success) return(integer(0L))
  lo <- 1L
  hi <- n
  while (lo < hi) {
    mid <- floor((lo + hi) / 2)
    if (count_at(mid) <= target_success) hi <- mid else lo <- mid + 1L
  }
  first <- lo
  if (count_at(first) != target_success) return(integer(0L))
  lo <- first
  hi <- n
  while (lo < hi) {
    mid <- ceiling((lo + hi) / 2)
    if (count_at(mid) >= target_success) lo <- mid else hi <- mid - 1L
  }
  last <- lo
  if (count_at(last) != target_success) return(integer(0L))
  seq.int(first, last)
}

n_used_vector <- function(cache, pp_threshold, lambda_prime, interim_counts, n_max) {
  boundary <- pp_threshold * (interim_counts / n_max)^lambda_prime
  valid_idx <- which(cache$valid)
  below <- sweep(cache$pp[valid_idx, , drop = FALSE], 2L, boundary, FUN = "<=")
  stopped <- rowSums(below) > 0L
  first_stop <- max.col(below, ties.method = "first")
  n_used <- rep(n_max, cache$attempted)
  n_used[valid_idx[stopped]] <- interim_counts[first_stop[stopped]]
  n_used
}

screen_case <- function(case) {
  design <- designs[[case]]
  h0 <- prepare_cache(cache_path(case, "H0"))
  h1 <- prepare_cache(cache_path(case, "H1"))
  rows <- vector("list", length(lambda_grid))
  for (j in seq_along(lambda_grid)) {
    idx <- find_exact_target_pp(h0, lambda_grid[j], design$interim, design$n_max, 500L)
    if (!length(idx)) next
    rows[[j]] <- do.call(rbind, lapply(idx, function(i) {
      m0 <- scenario_metrics(h0, pp_grid[i], lambda_grid[j], design$interim, design$n_max)
      m1 <- scenario_metrics(h1, pp_grid[i], lambda_grid[j], design$interim, design$n_max)
      data.frame(
        case = case,
        candidate_id = i + (j - 1L) * length(pp_grid),
        c = pp_grid[i],
        gamma = lambda_grid[j],
        H0_T1E = unname(m0["rate"]),
        H0_PET = unname(m0["pet"]),
        H0_ESS = unname(m0["ess"]),
        H1_Power = unname(m1["rate"]),
        H1_PET = unname(m1["pet"]),
        H1_ESS = unname(m1["ess"]),
        stringsAsFactors = FALSE
      )
    }))
  }
  exact <- do.call(rbind, rows[!vapply(rows, is.null, logical(1L))])
  if (is.null(exact) || !nrow(exact)) stop("No exact 5% candidate: ", case)
  exact <- exact[exact$H1_Power >= 0.80 & exact$c >= 0.85 & exact$gamma >= 0.01, , drop = FALSE]
  if (!nrow(exact)) stop("No candidate meets power, c and gamma gates: ", case)
  list(eligible = exact, h0 = h0, design = design)
}

select_case <- function(screened) {
  eligible <- screened$eligible
  design <- screened$design
  h0 <- screened$h0
  strict_id <- eligible$candidate_id[which.min(eligible$H0_ESS)]
  strict <- eligible[eligible$candidate_id == strict_id, , drop = FALSE]
  if (nrow(strict) != 1L) stop("Strict minimum ESS is not unique: ", eligible$case[1])
  n_ref <- n_used_vector(h0, strict$c, strict$gamma, design$interim, design$n_max)
  se <- vapply(seq_len(nrow(eligible)), function(i) {
    n_i <- n_used_vector(h0, eligible$c[i], eligible$gamma[i], design$interim, design$n_max)
    stats::sd(n_i - n_ref) / sqrt(h0$attempted)
  }, numeric(1L))
  near <- eligible[eligible$H0_ESS <= strict$H0_ESS + 2 * se, , drop = FALSE]
  boundary <- near$c <= min(pp_grid) + 1e-9 | near$c >= max(pp_grid) - 1e-9 |
    near$gamma <= 0.01 + 1e-9 | near$gamma >= max(lambda_grid) - 1e-9
  near$boundary <- boundary
  near <- near[order(-near$H1_Power, -near$c, near$H1_PET, near$boundary, near$candidate_id), , drop = FALSE]
  chosen <- near[1L, , drop = FALSE]
  chosen$status <- "NEAR_OPTIMAL_2SE"
  chosen$eligible_n <- nrow(eligible)
  chosen$near_n <- nrow(near)
  chosen
}

