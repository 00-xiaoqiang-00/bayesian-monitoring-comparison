# Runtime checks and cache context for formal RMUT analyses.
source(file.path(.rmut_directory, "..", "..", "common", "rmut_source_paths.R"), local = TRUE)
rmut_stage_version <- "RMUT_ROWPAIR_20261010_V2"

rmut_source_hashes <- function(directory) {
  paths <- rmut_dependency_paths(directory)
  setNames(unname(tools::md5sum(paths)), names(paths))
}

rmut_preflight <- function(directory, case, scope) {
  hashes <- rmut_source_hashes(directory)
  for (path in rmut_dependency_paths(directory)) parse(path, encoding = "UTF-8")
  missing <- c("brms", "rstan", "posterior", "survival", "MASS")[!vapply(c("brms", "rstan", "posterior", "survival", "MASS"), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  cat("PREFLIGHT_OK case=", case, " scope=", scope, " version=", rmut_stage_version, "\n", sep = "")
  invisible(hashes)
}

rmut_resolve_parameters <- function(case, candidate, c, gamma) {
  manifest <- Sys.getenv("RMUT_LOCKED_PARAMETER_MANIFEST", unset = Sys.getenv("RMUT_PARAMETER_MANIFEST", unset = ""))
  status <- "PARAMETER_DEFAULTS"
  if (nzchar(manifest)) {
    tab <- read.delim(manifest, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("case", "candidate_id", "c", "gamma")
    if (!all(required %in% names(tab))) stop("Invalid RMUT parameter manifest")
    if ("method" %in% names(tab)) tab <- tab[tab$method == "RMUT",,drop=FALSE]
    tab$case[tab$case == "IMpower131"] <- "IMpower"
    tab$case[tab$case == "KEYNOTE-598"] <- "KEYNOTE598"
    if (nrow(tab) != 4L || anyDuplicated(tab$case) || !setequal(tab$case,c("IMpower","KEYNOTE598","OAK","survival_colon"))) stop("Manifest must contain exactly four unique RMUT cases")
    ids <- suppressWarnings(as.numeric(tab$candidate_id))
    if (any(!is.finite(ids)) || any(ids<1|ids!=trunc(ids)) || any(!is.finite(tab$c)) || any(tab$c<=0|tab$c>=1) || any(!is.finite(tab$gamma)) || any(tab$gamma<0)) stop("Invalid manifest parameters")
    row <- tab[tab$case == case, , drop = FALSE]
    if (nrow(row) != 1L) stop("Manifest must contain one row for case ", case)
    candidate <- as.character(row$candidate_id); c <- row$c; gamma <- row$gamma
    status <- if ("status" %in% names(row)) as.character(row$status) else "EXTERNAL_PARAMETER_MANIFEST"
  }
  if (length(c) != 1L || !is.finite(c) || c <= 0 || c >= 1 || length(gamma) != 1L || !is.finite(gamma) || gamma < 0 || length(candidate) != 1L || !nzchar(candidate)) stop("Invalid RMUT decision parameters")
  list(candidate_id = candidate, c = c, gamma = gamma, status = status)
}

rmut_cache_context <- function(directory, case, scope, scenario, n, data_seed, model_seed, weight, tau) {
  list(version = rmut_stage_version, case = case, scope = scope, scenario = scenario,
       n = n, data_seed = data_seed, model_seed = model_seed, weight = weight, tau = tau,
       source_hashes = rmut_source_hashes(directory),
       entry_hash = unname(tools::md5sum(.rmut_entry_file)), pair_seed_base = 720260922L,
       normal_kind = "Inversion")
}

rmut_stamp_cache <- function(stats, context) {
  attr(stats, "rmut_corrected_context") <- context
  stats
}

rmut_corrected_cache_is_current <- function(stats, n, n_interims, context) {
  if (!is.list(stats) || length(stats) != n || !identical(attr(stats, "rmut_corrected_context"), context)) return(FALSE)
  ids <- vapply(stats, function(x) {
    id <- if (is.list(x)) x$attempt_id else NULL
    if (!is.numeric(id) || length(id) != 1L || !is.finite(id) || id != trunc(id) || id < 1 || id > n) return(NA_integer_)
    as.integer(id)
  }, integer(1))
  if (!identical(ids,seq_len(n))) return(FALSE)
  all(vapply(stats, function(x) {
    is.list(x) && length(x$PP) == n_interims && length(x$PP_final) == 1L &&
      is.logical(x$fit_ok) && length(x$fit_ok) == 1L && !is.na(x$fit_ok) &&
      is.data.frame(x$pairing_records) && nrow(x$pairing_records) == n_interims + 1L &&
      all(c("stage","pair_seed","n_death","n_composite","n_paired","unequal") %in% names(x$pairing_records)) &&
      identical(x$pairing_records$stage,c(paste0("interim_",seq_len(n_interims)),"final")) &&
      all(is.finite(x$pairing_records$pair_seed)) &&
      (!x$fit_ok || (all(is.finite(c(x$PP,x$PP_final))) && all(c(x$PP,x$PP_final) >= 0 & c(x$PP,x$PP_final) <= 1)))
  }, logical(1)))
}

rmut_draw_pair_record <- function(draws_d, draws_c, seed, stage, fit_d = NULL, fit_c = NULL) {
  da <- attr(fit_d,"rmut_fit_audit"); ca <- attr(fit_c,"rmut_fit_audit")
  nd <- if (!is.null(draws_d)) nrow(draws_d) else if (!is.null(da)) da$posterior_rows else 0L
  nc <- if (!is.null(draws_c)) nrow(draws_c) else if (!is.null(ca)) ca$posterior_rows else 0L
  record <- data.frame(stage = stage, pair_seed = seed, n_death = nd, n_composite = nc,
             n_paired = min(nd,nc), unequal = nd != nc,
             death_attempt = if (is.null(attr(fit_d,"iv_attempt"))) NA_character_ else attr(fit_d,"iv_attempt"),
             composite_attempt = if (is.null(attr(fit_c,"iv_attempt"))) NA_character_ else attr(fit_c,"iv_attempt"), stringsAsFactors = FALSE)
  attr(record,"death_fit_audit") <- da;attr(record,"composite_fit_audit") <- ca
  record
}

rmut_validation_result <- function(j, interim_pp, final_pp, records) {
  pairing <- if (length(records)) do.call(rbind, records) else data.frame()
  ok <- all(is.finite(c(interim_pp,final_pp))) && length(final_pp) == 1L
  list(attempt_id = as.integer(j), fit_ok = ok, failure_reason = if (ok) "" else "one_or_more_AFT_or_PP_failures",
       PP = interim_pp, PP_final = final_pp, pairing_records = pairing,
       fit_pair_audit = lapply(records,function(x) list(death=attr(x,"death_fit_audit"),composite=attr(x,"composite_fit_audit"))))
}

rmut_attach_calibration_pairing <- function(result, case_id, scenario_code) {
  records <- result$fit_records
  stages <- unique(records$stage)
  result$pairing_records <- do.call(rbind, lapply(stages, function(stage) {
    d <- records[records$stage == stage & records$component == "death", , drop = FALSE]
    c <- records[records$stage == stage & records$component == "composite", , drop = FALSE]
    look <- if (stage == "final") length(result$PP) + 1L else as.integer(sub("interim_", "", stage))
    data.frame(stage = stage, pair_seed = rmut_pair_seed(result$attempt_id,look,scenario_code,case_id),
               n_death = d$n_draws, n_composite = c$n_draws, n_paired = min(d$n_draws,c$n_draws),
               unequal = d$n_draws != c$n_draws, death_attempt = d$attempt, composite_attempt = c$attempt,
               stringsAsFactors = FALSE)
  }))
  result
}

rmut_write_validation_audit <- function(stats, output_dir, scenario, cutoffs, final_cutoff, interim_counts, nmax) {
  valid <- vapply(stats, function(x) isTRUE(x$fit_ok), logical(1))
  stop_at <- vapply(stats, function(x) {
    if (!isTRUE(x$fit_ok)) return(0L)
    hits <- which(x$PP <= cutoffs)
    if (length(hits)) hits[1L] else 0L
  }, integer(1))
  go <- vapply(seq_along(stats), function(i) valid[i] && stop_at[i] == 0L && stats[[i]]$PP_final > final_cutoff, logical(1))
  ess <- rep(nmax,length(stats)); ess[stop_at > 0L] <- interim_counts[stop_at[stop_at > 0L]]
  summary <- data.frame(scenario = scenario, total = length(stats), valid = sum(valid), failures = sum(!valid), go = sum(go),
                        go_probability = mean(go), PET = mean(stop_at > 0L), ESS = mean(ess),
                        final_no_go = sum(valid & stop_at == 0L & !go),
                        conditional_valid_go_probability = if (any(valid)) mean(go[valid]) else NA_real_,
                        failure_rule = "no-go; ESS=Nmax; fixed total denominator")
  pp <- do.call(rbind,lapply(seq_along(stats),function(i) data.frame(replicate_id=i,fit_ok=valid[i],
       as.list(setNames(stats[[i]]$PP,paste0("PP_interim_",seq_along(cutoffs)))),PP_final=stats[[i]]$PP_final,
       stop_look=stop_at[i],go=go[i],ESS=ess[i],check.names=FALSE)))
  write.csv(summary,file.path(output_dir,paste0(scenario,"_summary.csv")),row.names=FALSE)
  write.csv(pp,file.path(output_dir,paste0(scenario,"_replicates.csv")),row.names=FALSE)
  pair_rows <- lapply(seq_along(stats),function(i) {
    rows <- stats[[i]]$pairing_records
    if (!nrow(rows)) return(NULL)
    cbind(replicate_id=i,rows)
  })
  pairs <- do.call(rbind,pair_rows)
  if (!is.null(pairs)) write.csv(pairs,file.path(output_dir,paste0(scenario,"_pairing.csv")),row.names=FALSE)
  invisible(summary)
}

# Posterior draws for calibration and formal validation.
extract_valid_aft_draws <- function(fit) {
  if (inherits(fit, "try-error")) return(NULL)

  variables <- c("b_Intercept", "b_group", "sigma")
  draws <- try(posterior::as_draws_array(fit, variable = variables), silent = TRUE)
  if (inherits(draws, "try-error")) return(NULL)

  diagnostics <- posterior::summarise_draws(draws)
  if (any(!is.finite(diagnostics$rhat)) ||
      any(diagnostics$rhat > 1.01) ||
      any(diagnostics$ess_bulk < 400) ||
      any(diagnostics$ess_tail < 400)) {
    return(NULL)
  }

  as.data.frame(fit)[, variables, drop = FALSE]
}
