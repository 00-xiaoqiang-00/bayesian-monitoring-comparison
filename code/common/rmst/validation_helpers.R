# RMST simulation decision and probability summaries.
#
# Converts replicate-level posterior probabilities to decision tables,
# probability estimates, MCSEs, Wilson intervals and final no-go counts.

rmst_iv_scenario_key <- function(label) {
  if (grepl("H0", label, fixed = TRUE)) {
    return("H0")
  }
  if (grepl("H1_observed", label, fixed = TRUE) ||
      grepl("真实", label, fixed = TRUE)) {
    return("H1_observed")
  }
  "H1_standard"
}

rmst_iv_wilson_interval <- function(successes, denominator, conf_level = 0.95) {
  if (!is.finite(denominator) || denominator <= 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  z_value <- qnorm(1 - (1 - conf_level) / 2)
  estimate <- successes / denominator
  adjustment <- z_value^2 / denominator
  center <- (estimate + adjustment / 2) / (1 + adjustment)
  half_width <- (
    z_value *
      sqrt(
        estimate * (1 - estimate) / denominator +
          z_value^2 / (4 * denominator^2)
      )
  ) / (1 + adjustment)

  c(
    lower = max(0, center - half_width),
    upper = min(1, center + half_width)
  )
}

rmst_iv_probability_row <- function(
    scenario,
    analysis_set,
    metric,
    numerator,
    denominator) {
  estimate <- if (denominator > 0L) numerator / denominator else NA_real_
  mcse <- if (denominator > 0L) {
    sqrt(estimate * (1 - estimate) / denominator)
  } else {
    NA_real_
  }
  interval <- rmst_iv_wilson_interval(numerator, denominator)

  data.frame(
    scenario = scenario,
    analysis_set = analysis_set,
    metric = metric,
    numerator = as.integer(numerator),
    denominator = as.integer(denominator),
    estimate = estimate,
    mcse = mcse,
    mc_interval_method = "Wilson 95%",
    mc_ci_lower = unname(interval["lower"]),
    mc_ci_upper = unname(interval["upper"]),
    stringsAsFactors = FALSE
  )
}

rmst_iv_replace_scenario_rows <- function(existing, incoming, key_columns) {
  if (is.null(existing) || nrow(existing) == 0L) {
    return(incoming)
  }

  keep <- rep(TRUE, nrow(existing))
  for (row_index in seq_len(nrow(incoming))) {
    same_row <- rep(TRUE, nrow(existing))
    for (column_name in key_columns) {
      same_row <- same_row &
        as.character(existing[[column_name]]) ==
          as.character(incoming[[column_name]][row_index])
    }
    keep[same_row] <- FALSE
  }

  rbind(existing[keep, , drop = FALSE], incoming)
}

rmst_iv_update_combined_csv <- function(
    incoming,
    combined_file,
    key_columns) {
  existing <- NULL
  if (file.exists(combined_file)) {
    existing <- try(
      read.csv(
        combined_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      silent = TRUE
    )
    if (inherits(existing, "try-error")) {
      existing <- NULL
    }
  }

  combined <- rmst_iv_replace_scenario_rows(
    existing = existing,
    incoming = incoming,
    key_columns = key_columns
  )
  write.csv(combined, combined_file, row.names = FALSE, na = "")
}

rmst_iv_write_a05_a06_outputs <- function(
    stats_list,
    label,
    c_primes,
    N.max,
    interim_fracs,
    PP_final_threshold,
    lambda_prime_bayes,
    output_dir) {
  scenario <- rmst_iv_scenario_key(label)
  n_total <- length(stats_list)
  n_looks <- length(interim_fracs)

  if (n_total == 0L) {
    stop(sprintf("场景 %s 没有重复级结果。", scenario))
  }

  pp_matrix <- matrix(
    NA_real_,
    nrow = n_total,
    ncol = n_looks,
    dimnames = list(NULL, paste0("PP_interim_", seq_len(n_looks)))
  )
  pp_final <- rep(NA_real_, n_total)
  fit_ok <- rep(FALSE, n_total)
  failure_reason <- rep("", n_total)
  early_futility <- rep(NA, n_total)
  stop_stage <- rep(NA_integer_, n_total)
  sample_size_used <- rep(NA_integer_, n_total)
  go <- rep(NA, n_total)
  final_no_go <- rep(NA, n_total)

  for (replicate_index in seq_len(n_total)) {
    stats <- stats_list[[replicate_index]]
    structure_ok <- is.list(stats) &&
      length(stats$PP) == n_looks &&
      length(stats$PP_final) == 1L

    if (structure_ok) {
      pp_matrix[replicate_index, ] <- as.numeric(stats$PP)
      pp_final[replicate_index] <- as.numeric(stats$PP_final)
    }

    if (structure_ok) {
      first_nonfinite <- which(
        !is.finite(pp_matrix[replicate_index, ])
      )
      evaluable_looks <- if (length(first_nonfinite) > 0L) {
        seq_len(first_nonfinite[1] - 1L)
      } else {
        seq_len(n_looks)
      }
      stop_indices <- evaluable_looks[
        pp_matrix[replicate_index, evaluable_looks] <=
          c_primes[evaluable_looks]
      ]
      early_outcome_known <- length(stop_indices) > 0L
      reached_final_with_result <-
        length(first_nonfinite) == 0L &&
        is.finite(pp_final[replicate_index])
      fit_ok[replicate_index] <-
        early_outcome_known || reached_final_with_result
    }

    if (!fit_ok[replicate_index]) {
      reason <- if (is.list(stats) && length(stats$failure_reason) == 1L) {
        as.character(stats$failure_reason)
      } else {
        "missing_or_nonfinite_replicate_output"
      }
      if (!nzchar(reason)) {
        reason <- "missing_or_nonfinite_replicate_output"
      }
      failure_reason[replicate_index] <- reason
      next
    }

    early_futility[replicate_index] <- length(stop_indices) > 0L

    if (early_futility[replicate_index]) {
      first_stop <- stop_indices[1]
      stop_stage[replicate_index] <- first_stop
      sample_size_used[replicate_index] <- floor(
        N.max * interim_fracs[first_stop]
      )
    } else {
      stop_stage[replicate_index] <- n_looks + 1L
      sample_size_used[replicate_index] <- as.integer(N.max)
    }

    go[replicate_index] <- !early_futility[replicate_index] &&
      pp_final[replicate_index] > PP_final_threshold
    final_no_go[replicate_index] <- !early_futility[replicate_index] &&
      pp_final[replicate_index] <= PP_final_threshold
  }

  decision_table <- data.frame(
    replicate_id = seq_len(n_total),
    scenario = scenario,
    fit_ok = fit_ok,
    failure_reason = failure_reason,
    pp_matrix,
    PP_final = pp_final,
    early_futility = early_futility,
    stop_stage = stop_stage,
    sample_size_used = sample_size_used,
    go = go,
    final_no_go = final_no_go,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  valid_index <- which(fit_ok)
  n_valid <- length(valid_index)
  n_failed <- n_total - n_valid

  valid_pet_count <- sum(early_futility[valid_index])
  valid_go_count <- sum(go[valid_index])
  valid_final_no_go_count <- sum(final_no_go[valid_index])

  a05_table <- rbind(
    rmst_iv_probability_row(
      scenario, "valid_only", "PET",
      valid_pet_count, n_valid
    ),
    rmst_iv_probability_row(
      scenario, "valid_only", "go_probability",
      valid_go_count, n_valid
    ),
    rmst_iv_probability_row(
      scenario, "valid_only", "final_no_go",
      valid_final_no_go_count, n_valid
    ),
    rmst_iv_probability_row(
      scenario, "all_replicates_failure_as_no_go", "PET",
      valid_pet_count, n_total
    ),
    rmst_iv_probability_row(
      scenario, "all_replicates_failure_as_no_go", "go_probability",
      valid_go_count, n_total
    ),
    rmst_iv_probability_row(
      scenario, "all_replicates_failure_as_no_go", "final_no_go",
      valid_final_no_go_count + n_failed, n_total
    )
  )

  a06_table <- data.frame(
    scenario = scenario,
    analysis_set = c(
      "valid_only",
      "all_replicates_failure_as_no_go"
    ),
    denominator = c(n_valid, n_total),
    PET_count = c(valid_pet_count, valid_pet_count),
    go_count = c(valid_go_count, valid_go_count),
    final_no_go_count = c(
      valid_final_no_go_count,
      valid_final_no_go_count + n_failed
    ),
    fit_failure_count = c(0L, n_failed),
    partition_count = c(
      valid_pet_count + valid_go_count + valid_final_no_go_count,
      valid_pet_count + valid_go_count +
        valid_final_no_go_count + n_failed
    ),
    partition_probability = c(
      if (n_valid > 0L) {
        (
          valid_pet_count + valid_go_count +
            valid_final_no_go_count
        ) / n_valid
      } else {
        NA_real_
      },
      (
        valid_pet_count + valid_go_count +
          valid_final_no_go_count + n_failed
      ) / n_total
    ),
    identity_ok = c(
      n_valid > 0L &&
        valid_pet_count + valid_go_count +
          valid_final_no_go_count == n_valid,
      valid_pet_count + valid_go_count +
        valid_final_no_go_count + n_failed == n_total
    ),
    PP_final_threshold = PP_final_threshold,
    lambda_prime_bayes = lambda_prime_bayes,
    futility_operator = "PP <= b_k",
    stringsAsFactors = FALSE
  )

  a07_table <- data.frame(
    scenario = scenario,
    analysis_set = "valid_only",
    n_valid = n_valid,
    mean_ESS = if (n_valid > 0L) {
      mean(sample_size_used[valid_index])
    } else {
      NA_real_
    },
    ESS_over_Nmax = if (n_valid > 0L) {
      mean(sample_size_used[valid_index]) / N.max
    } else {
      NA_real_
    },
    stop_stage_1 = sum(stop_stage[valid_index] == 1L),
    stop_stage_2 = sum(stop_stage[valid_index] == 2L),
    stop_stage_3 = sum(stop_stage[valid_index] == 3L),
    reached_final = sum(stop_stage[valid_index] == n_looks + 1L),
    stringsAsFactors = FALSE
  )

  decision_rds <- file.path(
    output_dir,
    sprintf("replicate_level_%s.rds", scenario)
  )
  decision_csv <- file.path(
    output_dir,
    sprintf("replicate_level_%s.csv", scenario)
  )
  a05_file <- file.path(
    output_dir,
    sprintf("A05_probability_metrics_%s.csv", scenario)
  )
  a06_file <- file.path(
    output_dir,
    sprintf("A06_partition_check_%s.csv", scenario)
  )
  a07_file <- file.path(
    output_dir,
    sprintf("A07_ess_stop_stage_%s.csv", scenario)
  )

  saveRDS(decision_table, decision_rds)
  write.csv(decision_table, decision_csv, row.names = FALSE, na = "")
  write.csv(a05_table, a05_file, row.names = FALSE, na = "")
  write.csv(a06_table, a06_file, row.names = FALSE, na = "")
  write.csv(a07_table, a07_file, row.names = FALSE, na = "")

  rmst_iv_update_combined_csv(
    incoming = a05_table,
    combined_file = file.path(
      output_dir,
      "A05_probability_metrics_all_scenarios.csv"
    ),
    key_columns = c("scenario", "analysis_set", "metric")
  )
  rmst_iv_update_combined_csv(
    incoming = a06_table,
    combined_file = file.path(
      output_dir,
      "A06_partition_check_all_scenarios.csv"
    ),
    key_columns = c("scenario", "analysis_set")
  )
  rmst_iv_update_combined_csv(
    incoming = a07_table,
    combined_file = file.path(
      output_dir,
      "A07_ess_stop_stage_all_scenarios.csv"
    ),
    key_columns = c("scenario", "analysis_set")
  )

  list(
    decision_table = decision_table,
    a05_table = a05_table,
    a06_table = a06_table,
    a07_table = a07_table,
    n_total = n_total,
    n_valid = n_valid,
    n_failed = n_failed,
    valid_indices = valid_index
  )
}
