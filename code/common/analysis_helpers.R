# Shared monitoring helpers; source into the analysis entry environment.

monitoring_seed_from_label <- function(label) {
  label_codes <- utf8ToInt(enc2utf8(as.character(label)))
  label_offset <- if (length(label_codes)) {
    sum(seq_along(label_codes) * label_codes)
  } else {
    0
  }
  as.integer((iv_data_seed_base + label_offset) %% .Machine$integer.max)
}

calibration_seed_from_label <- function(label) {
  label_codes <- utf8ToInt(enc2utf8(as.character(label)))
  label_offset <- if (length(label_codes)) {
    sum(seq_along(label_codes) * label_codes)
  } else {
    0
  }
  as.integer((cal_data_seed_base + label_offset) %% .Machine$integer.max)
}

monitoring_rho_seed_from_label <- function(label) {
  scenario_tag <- if (grepl("H0", label, fixed = TRUE)) "H0" else if (
    grepl("H1_observed|Realcase", label)
  ) "H1_observed" else "H1_standard"
  rho_data_seed(iv_master_seed, iv_case, scenario_tag, iv_rho_level)
}

monitoring_rmst_weibull <- function(mu, alpha, tau) {
  scale <- mu / gamma(1 + 1 / alpha)
  x <- (tau / scale)^alpha
  
  res <- mu * pgamma(x, shape = 1 / alpha, rate = 1)
  
  # 检查 RMST 是否在 0 至 tau 范围内。
  stopifnot(all(res >= 0))
  stopifnot(all(res <= tau + 1e-8))
  
  return(res)
}

monitoring_valid_group_draws <- function(fit) {
  if (inherits(fit, "try-error")) return(NULL)

  draws <- try(posterior::as_draws_array(fit, variable = "b_group"), silent = TRUE)
  if (inherits(draws, "try-error")) return(NULL)

  diagnostics <- posterior::summarise_draws(draws)
  if (any(!is.finite(diagnostics$rhat)) ||
      any(diagnostics$rhat > 1.01) ||
      any(diagnostics$ess_bulk < 400) ||
      any(diagnostics$ess_tail < 400)) {
    return(NULL)
  }

  posterior::as_draws_matrix(draws)[, "b_group"]
}

monitoring_extract_valid_pp <- function(fit) {
  if (inherits(fit, "try-error")) return(NA_real_)

  draws <- try(posterior::as_draws_array(fit, variable = "b_group"), silent = TRUE)
  if (inherits(draws, "try-error")) return(NA_real_)

  diagnostics <- posterior::summarise_draws(draws)
  if (any(!is.finite(diagnostics$rhat)) ||
      any(diagnostics$rhat > 1.01) ||
      any(diagnostics$ess_bulk < 400) ||
      any(diagnostics$ess_tail < 400)) {
    return(NA_real_)
  }

  draw_matrix <- posterior::as_draws_matrix(draws)
  mean(draw_matrix[, "b_group"] < 0)
}

monitoring_write_calibration_output <- function() {
  cat("=========================================================================\n")
  cat("                   自适应设计决策参数校准网格报告\n")
  cat("=========================================================================\n")
  cat(sprintf("有效模拟次数：H0 = %d, H1 = %d\n", N.valid_H0, N.valid_H1))
  cat("候选网格汇总；正式决策参数由锁定参数表提供。\n\n")
  
  cat("【所有候选参数表现列表】\n")
  cat(sprintf("%-10s %-8s %-10s %-12s %-10s %-10s %-12s %-10s\n", 
              "PP_Thresh", "Lambda", "H0_T1E", "H0_StopRate", "H0_ESS", "H1_Power", "H1_StopRate", "H1_ESS"))
  for (i in 1:nrow(candidates_df)) {
    row <- candidates_df[i, ]
    cat(sprintf("%-10.2f %-8.2f %-10s %-12s %-10.1f %-10s %-12s %-10.1f\n",
                row$PP_Thresh, row$Lambda, 
                sprintf("%.1f%%", row$H0_T1E * 100), 
                sprintf("%.1f%%", row$H0_StopRate * 100), 
                row$H0_ESS, 
                sprintf("%.1f%%", row$H1_Power * 100), 
                sprintf("%.1f%%", row$H1_StopRate * 100), 
                row$H1_ESS))
  }
  
}

monitoring_rho_write_metadata <- function() {
  writeLines(c(
    paste0("version=", iv_version), paste0("project=", iv_project_id),
    paste0("method=", iv_method), paste0("case=", iv_case),
    paste0("candidate_id=", iv_candidate_id), paste0("run_tag=", iv_run_tag),
    "sensitivity_only=true", "recalibration=false",
    "sensitivity_parameter=KM_target_refitted_latent_DGM_rho",
    "parameter_source=locked_KM_target_parameter_table",
    paste0("parameter_table=", iv_parameter_table_path),
    paste0("parameter_table_sha256=", iv_parameter_table_sha256),
    paste0("sensitivity_rho=", format(iv_sensitivity_rho, digits = 15)),
    paste0("baseline_rho=", format(iv_baseline_rho, digits = 15)),
    paste0("rho_level=", iv_rho_level),
    paste0("fixed_interim_schedule=", iv_fixed_interim_schedule),
    "scenarios=H0,H1_standard,H1_observed",
    "baseline_decision_parameters_locked=true",
    "realcase_or_H1_observed_not_used_for_parameter_selection=true",
    "cache_isolation=KM-target rho-specific output directory",
    paste0("locked_PP_final_threshold=", format(iv_locked_pp, digits = 15)),
    paste0("locked_lambda_prime_bayes=", format(iv_locked_lambda, digits = 15)),
    "futility_rule=PP <= b_k",
    "locked_RMUT_weight=0.5",
    paste0("master_seed=", iv_master_seed),
    paste0("formal_data_seed_base=", iv_formal_data_seed_base),
    paste0("formal_model_seed_base=", iv_formal_model_seed_base),
    paste0("data_seed_base=", iv_data_seed_base),
    paste0("model_seed_base=", iv_model_seed_base),
    "seed_policy=data stream is shared across methods for case/scenario/rho; model stream is method-specific; smoke uses explicit offsets",
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ), file.path(iv_output_dir, "RUN_METADATA.txt"))
}

monitoring_rho_write_completion_marker <- function() {
  if (sink.number(type = "output") != 0L) {
    stop("Refusing to write RUN_COMPLETE.txt while an output sink is active.")
  }
  marker_file <- file.path(iv_output_dir, "RUN_COMPLETE.txt")
  temporary_file <- tempfile(
    pattern = ".RUN_COMPLETE_", tmpdir = iv_output_dir, fileext = ".tmp"
  )
  on.exit({
    if (file.exists(temporary_file)) unlink(temporary_file)
  }, add = TRUE)
  writeLines(
    c(
      paste0("version=", iv_version),
      paste0("project=", iv_project_id),
      paste0("rho_level=", iv_rho_level),
      paste0("sensitivity_rho=", format(iv_sensitivity_rho, digits = 15)),
      paste0("run_tag=", iv_run_tag),
      paste0("completed_at=", format(
        Sys.time(), tz = "Asia/Shanghai", usetz = TRUE
      ))
    ),
    temporary_file
  )
  if (!isTRUE(file.rename(temporary_file, marker_file))) {
    stop("Failed to atomically rename completion marker: ", marker_file)
  }
  invisible(marker_file)
}

cox_IMpower_validation_scenario <- function(label, mu_trt, Sigma_trt, mu_ctrl, Sigma_ctrl) {
  scenario_tag <- ifelse(grepl("H0", label), "H0", ifelse(grepl("H1_standard", label), "H1_standard", "Realcase"))
  cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_tag, cache_tag)
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      mu.trt = mu_trt, Sigma.trt = Sigma_trt, 
      mu.ctrl = mu_ctrl, Sigma.ctrl = Sigma_ctrl, cens_upper = 1e6
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      df_j$time <- pmin(time_to_first, df_j$censor_t)
      df_j$status <- as.integer(time_to_first <= df_j$censor_t)
      
      interim_PP <- numeric(length(interim_fracs))
      original_censor_t <- df_j$censor_t
      simulation_ok <- TRUE
      
      for (k in seq_along(interim_fracs)) {
        n_k <- floor(N.max * interim_fracs[k])
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        valid_rows <- censor_t_interim > 0
        df_interim <- data.frame(
          time = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        fit_bayes_interim <- try(
          update(base_model, newdata = df_interim, chains = 4, cores = 1,
                 iter = bayes_iter, warmup = bayes_warmup,
                 seed = iv_model_seed_base + j, refresh = 0, silent = 2),
          silent = TRUE
        )
        post_samples <- get_valid_group_draws(fit_bayes_interim)
        if (is.null(post_samples)) {
          simulation_ok <- FALSE
          break
        }
        interim_PP[k] <- mean(post_samples < 0)
      }

      if (!simulation_ok) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      
      final_admin_followup <- pmax(41 - entry_times, 0)
      final_censor_t <- pmin(original_censor_t, final_admin_followup)
      final_time <- pmin(time_to_first, final_censor_t)
      final_status <- as.integer(time_to_first <= final_censor_t)
      final_valid_rows <- final_censor_t > 0
      df_final_clean <- data.frame(
        time = final_time[final_valid_rows],
        status = final_status[final_valid_rows],
        group = df_j$group[final_valid_rows]
      )
      fit_bayes_final <- try(
        update(base_model, newdata = df_final_clean, chains = 4, cores = 1,
               iter = bayes_iter, warmup = bayes_warmup,
               seed = iv_model_seed_base + j, refresh = 0, silent = 2),
        silent = TRUE
      )
      post_samples <- get_valid_group_draws(fit_bayes_final)
      if (is.null(post_samples)) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      PP_final <- mean(post_samples < 0)
      
      list(ok = TRUE, PP = interim_PP, PP_final = PP_final)
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  stats_list <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list)
  if (length(stats_list) == 0) stop("没有通过 MCMC 诊断的有效验证模拟。")
  cat(sprintf("  - 通过 MCMC 诊断的有效模拟: %d/%d\n",
              length(stats_list), N.sim_realcase))

  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes
  
  stopped_real <- sapply(stats_list, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(stats_list, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(N.max * interim_fracs[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(stats_list, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

cox_KEYNOTE598_validation_scenario <- function(label, scenario) {
  scenario_tag <- scenario
  cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_tag, cache_tag)
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      cens_upper = trial_settings$cens_upper,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      time_to_first <- pmin(df_j$death_t, df_j$progression_t)
      final_censor_t <- pmin(
        df_j$censor_t,
        keynote598_admin_followup(
          entry_times, trial_settings$final_analysis_time
        )
      )
      df_j$time <- pmin(time_to_first, final_censor_t)
      df_j$status <- as.integer(time_to_first <= final_censor_t)
      
      interim_PP <- numeric(length(interim_fracs))
      original_censor_t <- df_j$censor_t
      simulation_ok <- TRUE
      
      for (k in seq_along(interim_fracs)) {
    n_k <- interim_counts[k]
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        valid_rows <- censor_t_interim > 0
        df_interim <- data.frame(
          time = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        fit_bayes_interim <- try(
          update(base_model, newdata = df_interim, chains = 4, cores = 1,
                 iter = bayes_iter, warmup = bayes_warmup,
                 seed = iv_model_seed_base + j, refresh = 0, silent = 2),
          silent = TRUE
        )
        post_samples <- get_valid_group_draws(fit_bayes_interim)
        if (is.null(post_samples)) {
          simulation_ok <- FALSE
          break
        }
        interim_PP[k] <- mean(post_samples < 0)
      }

      if (!simulation_ok) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      
      df_final_clean <- df_j[, c("time", "status", "group")]
      fit_bayes_final <- try(
        update(base_model, newdata = df_final_clean, chains = 4, cores = 1,
               iter = bayes_iter, warmup = bayes_warmup,
               seed = iv_model_seed_base + j, refresh = 0, silent = 2),
        silent = TRUE
      )
      post_samples <- get_valid_group_draws(fit_bayes_final)
      if (is.null(post_samples)) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      PP_final <- mean(post_samples < 0)
      
      list(ok = TRUE, PP = interim_PP, PP_final = PP_final)
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  stats_list <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list)
  if (length(stats_list) == 0) stop("没有通过 MCMC 诊断的有效验证模拟。")
  cat(sprintf("  - 通过 MCMC 诊断的有效模拟: %d/%d\n",
              length(stats_list), N.sim_realcase))

  c_primes <- PP_final_threshold * (interim_counts / N.max)^lambda_prime_bayes
  
  stopped_real <- sapply(stats_list, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(stats_list, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(interim_counts[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(stats_list, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

cox_OAK_validation_scenario <- function(label, scenario) {
  scenario_tag <- scenario
  cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_tag, cache_tag)
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      cens_upper = 109.3,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      final_censor_t <- pmin(df_j$censor_t, oak_admin_followup(entry_times, 27.9))
      df_j$time <- pmin(time_to_first, final_censor_t)
      df_j$status <- as.integer(time_to_first <= final_censor_t)
      
      interim_PP <- numeric(length(interim_fracs))
      original_censor_t <- df_j$censor_t
      simulation_ok <- TRUE
      
      for (k in seq_along(interim_fracs)) {
        n_k <- floor(N.max * interim_fracs[k])
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        valid_rows <- censor_t_interim > 0
        df_interim <- data.frame(
          time = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        fit_bayes_interim <- try(
          update(base_model, newdata = df_interim, chains = 4, cores = 1,
                 iter = bayes_iter, warmup = bayes_warmup,
                 seed = iv_model_seed_base + j, refresh = 0, silent = 2),
          silent = TRUE
        )
        post_samples <- get_valid_group_draws(fit_bayes_interim)
        if (is.null(post_samples)) {
          simulation_ok <- FALSE
          break
        }
        interim_PP[k] <- mean(post_samples < 0)
      }

      if (!simulation_ok) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      
      df_final_clean <- df_j[, c("time", "status", "group")]
      fit_bayes_final <- try(
        update(base_model, newdata = df_final_clean, chains = 4, cores = 1,
               iter = bayes_iter, warmup = bayes_warmup,
               seed = iv_model_seed_base + j, refresh = 0, silent = 2),
        silent = TRUE
      )
      post_samples <- get_valid_group_draws(fit_bayes_final)
      if (is.null(post_samples)) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      PP_final <- mean(post_samples < 0)
      
      list(ok = TRUE, PP = interim_PP, PP_final = PP_final)
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  stats_list <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list)
  if (length(stats_list) == 0) stop("没有通过 MCMC 诊断的有效验证模拟。")
  cat(sprintf("  - 通过 MCMC 诊断的有效模拟: %d/%d\n",
              length(stats_list), N.sim_realcase))

  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes
  
  stopped_real <- sapply(stats_list, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(stats_list, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(N.max * interim_fracs[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(stats_list, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

cox_survival_colon_validation_scenario <- function(label, scenario) {
  scenario_tag <- scenario
  cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_tag, cache_tag)
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      cens_upper = trial_settings$cens_upper,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      final_censor_t <- pmin(
        df_j$censor_t,
        survival_colon_admin_followup(
          entry_times, trial_settings$final_analysis_time
        )
      )
      df_j$time <- pmin(time_to_first, final_censor_t)
      df_j$status <- as.integer(time_to_first <= final_censor_t)
      
      interim_PP <- numeric(length(interim_fracs))
      original_censor_t <- df_j$censor_t
      simulation_ok <- TRUE
      
      for (k in seq_along(interim_fracs)) {
        n_k <- floor(N.max * interim_fracs[k])
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        valid_rows <- censor_t_interim > 0
        df_interim <- data.frame(
          time = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        fit_bayes_interim <- try(
          update(base_model, newdata = df_interim, chains = 4, cores = 1,
                 iter = bayes_iter, warmup = bayes_warmup,
                 seed = iv_model_seed_base + j, refresh = 0, silent = 2),
          silent = TRUE
        )
        post_samples <- get_valid_group_draws(fit_bayes_interim)
        if (is.null(post_samples)) {
          simulation_ok <- FALSE
          break
        }
        interim_PP[k] <- mean(post_samples < 0)
      }

      if (!simulation_ok) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      
      df_final_clean <- df_j[, c("time", "status", "group")]
      fit_bayes_final <- try(
        update(base_model, newdata = df_final_clean, chains = 4, cores = 1,
               iter = bayes_iter, warmup = bayes_warmup,
               seed = iv_model_seed_base + j, refresh = 0, silent = 2),
        silent = TRUE
      )
      post_samples <- get_valid_group_draws(fit_bayes_final)
      if (is.null(post_samples)) {
        return(list(ok = FALSE, PP = interim_PP, PP_final = NA_real_))
      }
      PP_final <- mean(post_samples < 0)
      
      list(ok = TRUE, PP = interim_PP, PP_final = PP_final)
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  stats_list <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list)
  if (length(stats_list) == 0) stop("没有通过 MCMC 诊断的有效验证模拟。")
  cat(sprintf("  - 通过 MCMC 诊断的有效模拟: %d/%d\n",
              length(stats_list), N.sim_realcase))

  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes
  
  stopped_real <- sapply(stats_list, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(stats_list, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(N.max * interim_fracs[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(stats_list, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

monitoring_rebuild_colon_sigma <- function(Sigma_baseline, rho, label) {
  if (!is.matrix(Sigma_baseline) || !identical(dim(Sigma_baseline), c(2L, 2L)) ||
      any(!is.finite(Sigma_baseline)) || any(diag(Sigma_baseline) <= 0)) {
    stop(label, " baseline covariance is invalid.")
  }
  baseline_reconstructed_rho <- Sigma_baseline[1L, 2L] /
    sqrt(Sigma_baseline[1L, 1L] * Sigma_baseline[2L, 2L])
  if (!isTRUE(all.equal(baseline_reconstructed_rho, iv_baseline_rho,
                        tolerance = 1e-12)) ||
      !isTRUE(all.equal(Sigma_baseline[2L, 1L], Sigma_baseline[1L, 2L],
                        tolerance = 1e-12))) {
    stop(label, " does not reconstruct the declared locked baseline rho.")
  }
  Sigma <- Sigma_baseline
  covariance <- rho * sqrt(Sigma_baseline[1L, 1L] * Sigma_baseline[2L, 2L])
  Sigma[1L, 2L] <- covariance
  Sigma[2L, 1L] <- covariance
  if (!identical(diag(Sigma), diag(Sigma_baseline)) ||
      min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop(label, " candidate covariance is not positive-definite.")
  }
  Sigma
}

monitoring_rebuild_sigma_for_rho <- function(Sigma_baseline, rho, label) {
  if (!is.matrix(Sigma_baseline) || !identical(dim(Sigma_baseline), c(2L, 2L)) ||
      any(!is.finite(Sigma_baseline)) || any(diag(Sigma_baseline) <= 0)) {
    stop(label, " baseline covariance must be a finite 2x2 matrix with positive diagonal.")
  }
  Sigma <- Sigma_baseline
  covariance <- rho * sqrt(Sigma_baseline[1L, 1L] * Sigma_baseline[2L, 2L])
  Sigma[1L, 2L] <- covariance
  Sigma[2L, 1L] <- covariance
  if (!identical(diag(Sigma), diag(Sigma_baseline)) ||
      min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop(label, " rho covariance reconstruction failed.")
  }
  Sigma
}

monitoring_assert_baseline_sigma <- function(Sigma_baseline, label) {
  expected <- iv_baseline_rho * sqrt(prod(diag(Sigma_baseline)))
  if (abs(Sigma_baseline[1L, 2L] - expected) > 5e-5 ||
      abs(Sigma_baseline[2L, 1L] - expected) > 5e-5) {
    stop(label, " does not reconstruct the locked baseline rho.")
  }
  invisible(TRUE)
}

# RMST scenario calculations for the formal and interim schedules.

rmst_IMpower_validation_scenario <- function(label, mu_trt, Sigma_trt, mu_ctrl, Sigma_ctrl) {
  cache_suffix <- if (N.sim_realcase == 10000) "10k_v2_" else "v2_"
  cache_file <- sprintf("validate_stats_%s%s.rds", cache_suffix, ifelse(grepl("H0", label), "H0", ifelse(grepl("H1_sens2", label), "H1_sens2", ifelse(grepl("H1_sens", label), "H1_sens", ifelse(grepl("H1", label), "H1", "Realcase")))))
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      mu.trt = mu_trt, Sigma.trt = Sigma_trt, 
      mu.ctrl = mu_ctrl, Sigma.ctrl = Sigma_ctrl, cens_upper = 1e6
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      
      # 仅使用潜在非致死事件时间和死亡时间均大于零的记录。
      df_j <- df_j[df_j$recurrence_t > 0 & df_j$death_t > 0, ]
      if (nrow(df_j) < 10) {
        return(list(
          ok = FALSE,
          failure_reason = "insufficient_generated_rows",
          PP = rep(NA_real_, length(interim_fracs)),
          PP_final = NA_real_
        ))
      }
      
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      df_j$time <- pmin(time_to_first, df_j$censor_t)
      df_j$status <- as.integer(time_to_first <= df_j$censor_t)
      
      interim_PP <- rep(NA_real_, length(interim_fracs))
      original_censor_t <- df_j$censor_t
      
      # 对每个期中分析做自适应截断和模型拟合
      for (k in seq_along(interim_fracs)) {
        n_k <- floor(nrow(df_j) * interim_fracs[k])
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        
        valid_rows <- censor_t_interim > 0 & time_interim > 0
        
        df_interim <- data.frame(
          time   = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group  = df_j$group[1:n_k][valid_rows]
        )
        
        tau_k <- tau_fixed[k]
        
        if (is.finite(tau_k) && tau_k > 0) {
          fit_bayes_interim <- try(
            update(
              base_model,
              newdata = df_interim,
              recompile = FALSE,
              chains  = bayes_chains,
              cores   = 1,
              iter    = bayes_iter,
              warmup  = bayes_warmup,
              seed    = iv_model_seed_base + 10L * j + k,
              refresh = 0,
              silent  = 2,
              control = list(
                adapt_delta   = bayes_adapt,
                max_treedepth = bayes_treedepth
              )
            ),
            silent = TRUE
          )
          
          if (!inherits(fit_bayes_interim, "try-error")) {
            newdat_sim <- data.frame(group = c(0, 1))
            mu_draws_sim <- try(posterior_linpred(fit_bayes_interim, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
            alpha_draws_sim <- try(posterior_linpred(fit_bayes_interim, dpar = "shape", transform = TRUE), silent = TRUE)
            
            if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
              mu0_draws <- mu_draws_sim[, 1]
              mu1_draws <- mu_draws_sim[, 2]
              alpha_vec <- alpha_draws_sim[, 1]
              
              rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_k)
              rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_k)
              
              interim_PP[k] <- mean(rmst1_draws > rmst0_draws)
            }
          }
        }
      }
      
      # 按预设最终日历截点施加行政删失。
      Time.final <- 41
      obs_window_final <- pmax(Time.final - entry_times, 0)
      censor_t_final <- obs_window_final
      
      sub_df <- df_j[, c("time", "status", "group")]
      sub_df$time <- pmin(df_j$time, censor_t_final)
      sub_df$status <- as.integer(df_j$status == 1 & df_j$time <= censor_t_final)
      tau_final <- tau_fixed[length(tau_fixed)]
      
      PP_final <- NA_real_
      if (is.finite(tau_final) && tau_final > 0) {
        fit_bayes_final <- try(
          update(
            base_model,
            newdata = sub_df,
            recompile = FALSE,
            chains  = bayes_chains,
            cores   = 1,
            iter    = bayes_iter,
            warmup  = bayes_warmup,
            seed    = iv_model_seed_base + 10L * j + 9L,
            refresh = 0,
            silent  = 2,
            control = list(
              adapt_delta   = bayes_adapt,
              max_treedepth = bayes_treedepth
            )
          ),
          silent = TRUE
        )
        if (!inherits(fit_bayes_final, "try-error")) {
          newdat_sim <- data.frame(group = c(0, 1))
          mu_draws_sim <- try(posterior_linpred(fit_bayes_final, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
          alpha_draws_sim <- try(posterior_linpred(fit_bayes_final, dpar = "shape", transform = TRUE), silent = TRUE)
          
          if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
            mu0_draws <- mu_draws_sim[, 1]
            mu1_draws <- mu_draws_sim[, 2]
            alpha_vec <- alpha_draws_sim[, 1]
            
            rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_final)
            rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_final)
            
            PP_final <- mean(rmst1_draws > rmst0_draws)
          }
        }
      }
      
      fit_ok <- length(interim_PP) == length(interim_fracs) &&
        all(is.finite(interim_PP)) &&
        length(PP_final) == 1L &&
        is.finite(PP_final)
      list(
        ok = fit_ok,
        failure_reason = if (fit_ok) {
          ""
        } else {
          "mcmc_fit_or_posterior_extraction_failed"
        },
        PP = interim_PP,
        PP_final = PP_final
      )
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes

  iv_a05_a06 <- rmst_iv_write_a05_a06_outputs(
    stats_list = stats_list,
    label = label,
    c_primes = c_primes,
    N.max = N.max,
    interim_fracs = interim_fracs,
    PP_final_threshold = PP_final_threshold,
    lambda_prime_bayes = lambda_prime_bayes,
    output_dir = iv_output_dir
  )
  valid_stats <- stats_list[iv_a05_a06$valid_indices]
  if (length(valid_stats) == 0L) {
    stop(sprintf("场景 %s 没有有效的RMST拟合结果。", label))
  }
  
  stopped_real <- sapply(valid_stats, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(valid_stats, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(N.max * interim_fracs[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(valid_stats, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf(
    "  - 有效模拟次数             : %d/%d\n",
    iv_a05_a06$n_valid,
    iv_a05_a06$n_total
  ))
  cat(sprintf(
    "  - 拟合失败次数             : %d\n",
    iv_a05_a06$n_failed
  ))
  cat("  - A05输出                  : 概率、MCSE、Wilson 95% Monte Carlo区间\n")
  cat("  - A06输出                  : final no-go直接计数及PET+go+final no-go恒等核对\n")
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

rmst_KEYNOTE598_validation_scenario <- function(label, scenario) {
  cache_suffix <- if (N.sim_realcase == 10000) "10k_" else ""
  cache_file <- sprintf(
    "validate_stats_%s%s_KEYNOTE598_v2_tau4_6_9_12_priorInterceptLog20_weak.rds",
    cache_suffix, scenario
  )
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      cens_upper = trial_settings$cens_upper,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      
      # 仅使用潜在非致死事件时间和死亡时间均大于零的记录。
      df_j <- df_j[df_j$progression_t > 0 & df_j$death_t > 0, ]
      if (nrow(df_j) < 10) {
        return(list(
          ok = FALSE,
          failure_reason = "insufficient_generated_rows",
          PP = rep(NA_real_, length(interim_fracs)),
          PP_final = NA_real_
        ))
      }
      
      time_to_first <- pmin(df_j$death_t, df_j$progression_t)
      df_j$time <- pmin(time_to_first, df_j$censor_t)
      df_j$status <- as.integer(time_to_first <= df_j$censor_t)
      
      interim_PP <- rep(NA_real_, length(interim_fracs))
      original_censor_t <- df_j$censor_t
      
      # 对每个期中分析做自适应截断和模型拟合
      for (k in seq_along(interim_fracs)) {
    n_k <- interim_counts[k]
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        
        valid_rows <- censor_t_interim > 0 & time_interim > 0
        
        df_interim <- data.frame(
          time   = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group  = df_j$group[1:n_k][valid_rows]
        )
        
        tau_k <- tau_fixed[k]
        
        if (is.finite(tau_k) && tau_k > 0) {
          fit_bayes_interim <- try(
            update(
              base_model,
              newdata = df_interim,
              chains  = 4,
              cores   = 1,
              iter    = bayes_iter,
              warmup  = bayes_warmup,
              seed    = iv_model_seed_base + j,
              refresh = 0,
              silent  = 2,
              control = list(
                adapt_delta   = bayes_adapt,
                max_treedepth = bayes_treedepth
              )
            ),
            silent = TRUE
          )
          
          if (!inherits(fit_bayes_interim, "try-error")) {
            newdat_sim <- data.frame(group = c(0, 1))
            mu_draws_sim <- try(posterior_linpred(fit_bayes_interim, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
            alpha_draws_sim <- try(posterior_linpred(fit_bayes_interim, dpar = "shape", transform = TRUE), silent = TRUE)
            
            if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
              mu0_draws <- mu_draws_sim[, 1]
              mu1_draws <- mu_draws_sim[, 2]
              alpha_vec <- alpha_draws_sim[, 1]
              
              rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_k)
              rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_k)
              
              interim_PP[k] <- mean(rmst1_draws > rmst0_draws)
            }
          }
        }
      }
      
      # 按预设最终日历截点施加行政删失。
      Time.final <- trial_settings$final_analysis_time
      obs_window_final <- pmax(Time.final - entry_times, 0)
      censor_t_final <- obs_window_final
      
      sub_df <- df_j[, c("time", "status", "group")]
      sub_df$time <- pmin(df_j$time, censor_t_final)
      sub_df$status <- as.integer(df_j$status == 1 & df_j$time <= censor_t_final)
      tau_final <- tau_fixed[length(tau_fixed)]
      
      PP_final <- NA_real_
      if (is.finite(tau_final) && tau_final > 0) {
        fit_bayes_final <- try(
          update(
            base_model,
            newdata = sub_df,
            chains  = 4,
            cores   = 1,
            iter    = bayes_iter,
            warmup  = bayes_warmup,
            seed    = iv_model_seed_base + j,
            refresh = 0,
            silent  = 2,
            control = list(
              adapt_delta   = bayes_adapt,
              max_treedepth = bayes_treedepth
            )
          ),
          silent = TRUE
        )
        if (!inherits(fit_bayes_final, "try-error")) {
          newdat_sim <- data.frame(group = c(0, 1))
          mu_draws_sim <- try(posterior_linpred(fit_bayes_final, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
          alpha_draws_sim <- try(posterior_linpred(fit_bayes_final, dpar = "shape", transform = TRUE), silent = TRUE)
          
          if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
            mu0_draws <- mu_draws_sim[, 1]
            mu1_draws <- mu_draws_sim[, 2]
            alpha_vec <- alpha_draws_sim[, 1]
            
            rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_final)
            rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_final)
            
            PP_final <- mean(rmst1_draws > rmst0_draws)
          }
        }
      }
      
      fit_ok <- length(interim_PP) == length(interim_fracs) &&
        all(is.finite(interim_PP)) &&
        length(PP_final) == 1L &&
        is.finite(PP_final)
      list(
        ok = fit_ok,
        failure_reason = if (fit_ok) {
          ""
        } else {
          "mcmc_fit_or_posterior_extraction_failed"
        },
        PP = interim_PP,
        PP_final = PP_final
      )
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  c_primes <- PP_final_threshold * (interim_counts / N.max)^lambda_prime_bayes

  iv_a05_a06 <- rmst_iv_write_a05_a06_outputs(
    stats_list = stats_list,
    label = label,
    c_primes = c_primes,
    N.max = N.max,
    interim_fracs = interim_fracs,
    PP_final_threshold = PP_final_threshold,
    lambda_prime_bayes = lambda_prime_bayes,
    output_dir = iv_output_dir
  )
  valid_stats <- stats_list[iv_a05_a06$valid_indices]
  if (length(valid_stats) == 0L) {
    stop(sprintf("场景 %s 没有有效的RMST拟合结果。", label))
  }
  
  stopped_real <- sapply(valid_stats, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(valid_stats, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(interim_counts[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(valid_stats, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf(
    "  - 有效模拟次数             : %d/%d\n",
    iv_a05_a06$n_valid,
    iv_a05_a06$n_total
  ))
  cat(sprintf(
    "  - 拟合失败次数             : %d\n",
    iv_a05_a06$n_failed
  ))
  cat("  - A05输出                  : 概率、MCSE、Wilson 95% Monte Carlo区间\n")
  cat("  - A06输出                  : final no-go直接计数及PET+go+final no-go恒等核对\n")
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

rmst_OAK_validation_scenario <- function(label, scenario) {
  cache_suffix <- if (N.sim_realcase == 10000) "10k_" else ""
  cache_file1 <- sprintf("validate_stats_%s%s_OAK_mix_v4_tau3_4p5_6_18.rds", cache_suffix, scenario)
  cache_file2 <- sprintf("validate_stats_%s%s_OAK_mix_v4_tau3_4p5_6_18_priorN0_1p5_explicit.rds", cache_suffix, scenario)
  cache_file <- file.path(iv_cache_dir, basename(cache_file2))
  use_cache <- FALSE
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      cens_upper = 109.3,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      
      # 仅使用潜在非致死事件时间和死亡时间均大于零的记录。
      df_j <- df_j[df_j$recurrence_t > 0 & df_j$death_t > 0, ]
      if (nrow(df_j) < 10) {
        return(list(
          ok = FALSE,
          failure_reason = "insufficient_generated_rows",
          PP = rep(NA_real_, length(interim_fracs)),
          PP_final = NA_real_
        ))
      }
      
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      df_j$time <- pmin(time_to_first, df_j$censor_t)
      df_j$status <- as.integer(time_to_first <= df_j$censor_t)
      
      interim_PP <- rep(NA_real_, length(interim_fracs))
      original_censor_t <- df_j$censor_t
      
      # 对每个期中分析做自适应截断和模型拟合
      for (k in seq_along(interim_fracs)) {
        n_k <- floor(nrow(df_j) * interim_fracs[k])
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        
        valid_rows <- censor_t_interim > 0 & time_interim > 0
        
        df_interim <- data.frame(
          time   = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group  = df_j$group[1:n_k][valid_rows]
        )
        
        tau_k <- tau_fixed[k]
        
        if (is.finite(tau_k) && tau_k > 0) {
          fit_bayes_interim <- try(
            update(
              base_model,
              newdata = df_interim,
              chains  = 4,
              cores   = 1,
              iter    = bayes_iter,
              warmup  = bayes_warmup,
              seed    = iv_model_seed_base + j,
              refresh = 0,
              silent  = 2,
              control = list(
                adapt_delta   = bayes_adapt,
                max_treedepth = bayes_treedepth
              )
            ),
            silent = TRUE
          )
          
          if (!inherits(fit_bayes_interim, "try-error")) {
            newdat_sim <- data.frame(group = c(0, 1))
            mu_draws_sim <- try(posterior_linpred(fit_bayes_interim, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
            alpha_draws_sim <- try(posterior_linpred(fit_bayes_interim, dpar = "shape", transform = TRUE), silent = TRUE)
            
            if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
              mu0_draws <- mu_draws_sim[, 1]
              mu1_draws <- mu_draws_sim[, 2]
              alpha_vec <- alpha_draws_sim[, 1]
              
              rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_k)
              rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_k)
              
              interim_PP[k] <- mean(rmst1_draws > rmst0_draws)
            }
          }
        }
      }
      
      # 按预设最终日历截点施加行政删失。
      Time.final <- 27.9
      obs_window_final <- pmax(Time.final - entry_times, 0)
      censor_t_final <- obs_window_final
      
      sub_df <- df_j[, c("time", "status", "group")]
      sub_df$time <- pmin(df_j$time, censor_t_final)
      sub_df$status <- as.integer(df_j$status == 1 & df_j$time <= censor_t_final)
      tau_final <- tau_fixed[length(tau_fixed)]
      
      PP_final <- NA_real_
      if (is.finite(tau_final) && tau_final > 0) {
        fit_bayes_final <- try(
          update(
            base_model,
            newdata = sub_df,
            chains  = 4,
            cores   = 1,
            iter    = bayes_iter,
            warmup  = bayes_warmup,
            seed    = iv_model_seed_base + j,
            refresh = 0,
            silent  = 2,
            control = list(
              adapt_delta   = bayes_adapt,
              max_treedepth = bayes_treedepth
            )
          ),
          silent = TRUE
        )
        if (!inherits(fit_bayes_final, "try-error")) {
          newdat_sim <- data.frame(group = c(0, 1))
          mu_draws_sim <- try(posterior_linpred(fit_bayes_final, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
          alpha_draws_sim <- try(posterior_linpred(fit_bayes_final, dpar = "shape", transform = TRUE), silent = TRUE)
          
          if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
            mu0_draws <- mu_draws_sim[, 1]
            mu1_draws <- mu_draws_sim[, 2]
            alpha_vec <- alpha_draws_sim[, 1]
            
            rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_final)
            rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_final)
            
            PP_final <- mean(rmst1_draws > rmst0_draws)
          }
        }
      }
      
      fit_ok <- length(interim_PP) == length(interim_fracs) &&
        all(is.finite(interim_PP)) &&
        length(PP_final) == 1L &&
        is.finite(PP_final)
      list(
        ok = fit_ok,
        failure_reason = if (fit_ok) {
          ""
        } else {
          "mcmc_fit_or_posterior_extraction_failed"
        },
        PP = interim_PP,
        PP_final = PP_final
      )
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes

  iv_a05_a06 <- rmst_iv_write_a05_a06_outputs(
    stats_list = stats_list,
    label = label,
    c_primes = c_primes,
    N.max = N.max,
    interim_fracs = interim_fracs,
    PP_final_threshold = PP_final_threshold,
    lambda_prime_bayes = lambda_prime_bayes,
    output_dir = iv_output_dir
  )
  valid_stats <- stats_list[iv_a05_a06$valid_indices]
  if (length(valid_stats) == 0L) {
    stop(sprintf("场景 %s 没有有效的RMST拟合结果。", label))
  }
  
  stopped_real <- sapply(valid_stats, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(valid_stats, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(N.max * interim_fracs[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(valid_stats, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf(
    "  - 有效模拟次数             : %d/%d\n",
    iv_a05_a06$n_valid,
    iv_a05_a06$n_total
  ))
  cat(sprintf(
    "  - 拟合失败次数             : %d\n",
    iv_a05_a06$n_failed
  ))
  cat("  - A05输出                  : 概率、MCSE、Wilson 95% Monte Carlo区间\n")
  cat("  - A06输出                  : final no-go直接计数及PET+go+final no-go恒等核对\n")
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

rmst_survival_colon_validation_scenario <- function(label, scenario) {
  cache_suffix <- if (N.sim_realcase == 10000) "10k_" else ""
  cache_file <- sprintf(
    "validate_stats_%s%s_survival_colon_v1_tau17p1_25p7_34p3_60_priorN0_1p5_explicit.rds",
    cache_suffix, scenario
  )
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (length(stats_list_temp) == N.sim_realcase) {
      cat(sprintf("\nCompatible cache found for [%s]  cache; loading it...\n", label))
      stats_list <- stats_list_temp
      use_cache <- TRUE
    }
  }
  
  if (!use_cache) {
    cat(sprintf("\n【开始进行 %s 场景仿真，模拟次数 = %d】...\n", label, N.sim_realcase))
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      cens_upper = trial_settings$cens_upper,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      
      # 仅使用潜在非致死事件时间和死亡时间均大于零的记录。
      df_j <- df_j[df_j$recurrence_t > 0 & df_j$death_t > 0, ]
      if (nrow(df_j) < 10) {
        return(list(
          ok = FALSE,
          failure_reason = "insufficient_generated_rows",
          PP = rep(NA_real_, length(interim_fracs)),
          PP_final = NA_real_
        ))
      }
      
      time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
      df_j$time <- pmin(time_to_first, df_j$censor_t)
      df_j$status <- as.integer(time_to_first <= df_j$censor_t)
      
      interim_PP <- rep(NA_real_, length(interim_fracs))
      original_censor_t <- df_j$censor_t
      
      # 对每个期中分析做自适应截断和模型拟合
      for (k in seq_along(interim_fracs)) {
        n_k <- floor(nrow(df_j) * interim_fracs[k])
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        
        time_interim <- pmin(time_to_first[1:n_k], censor_t_interim)
        status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
        
        valid_rows <- censor_t_interim > 0 & time_interim > 0
        
        df_interim <- data.frame(
          time   = time_interim[valid_rows],
          status = status_interim[valid_rows],
          group  = df_j$group[1:n_k][valid_rows]
        )
        
        tau_k <- tau_fixed[k]
        
        if (is.finite(tau_k) && tau_k > 0) {
          fit_bayes_interim <- try(
            update(
              base_model,
              newdata = df_interim,
              chains  = 4,
              cores   = 1,
              iter    = bayes_iter,
              warmup  = bayes_warmup,
              seed    = iv_model_seed_base + j,
              refresh = 0,
              silent  = 2,
              control = list(
                adapt_delta   = bayes_adapt,
                max_treedepth = bayes_treedepth
              )
            ),
            silent = TRUE
          )
          
          if (!inherits(fit_bayes_interim, "try-error")) {
            newdat_sim <- data.frame(group = c(0, 1))
            mu_draws_sim <- try(posterior_linpred(fit_bayes_interim, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
            alpha_draws_sim <- try(posterior_linpred(fit_bayes_interim, dpar = "shape", transform = TRUE), silent = TRUE)
            
            if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
              mu0_draws <- mu_draws_sim[, 1]
              mu1_draws <- mu_draws_sim[, 2]
              alpha_vec <- alpha_draws_sim[, 1]
              
              rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_k)
              rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_k)
              
              interim_PP[k] <- mean(rmst1_draws > rmst0_draws)
            }
          }
        }
      }
      
      # 按预设最终日历截点施加行政删失。
      Time.final <- trial_settings$final_analysis_time
      obs_window_final <- pmax(Time.final - entry_times, 0)
      censor_t_final <- obs_window_final
      
      sub_df <- df_j[, c("time", "status", "group")]
      sub_df$time <- pmin(df_j$time, censor_t_final)
      sub_df$status <- as.integer(df_j$status == 1 & df_j$time <= censor_t_final)
      tau_final <- tau_fixed[length(tau_fixed)]
      
      PP_final <- NA_real_
      if (is.finite(tau_final) && tau_final > 0) {
        fit_bayes_final <- try(
          update(
            base_model,
            newdata = sub_df,
            chains  = 4,
            cores   = 1,
            iter    = bayes_iter,
            warmup  = bayes_warmup,
            seed    = iv_model_seed_base + j,
            refresh = 0,
            silent  = 2,
            control = list(
              adapt_delta   = bayes_adapt,
              max_treedepth = bayes_treedepth
            )
          ),
          silent = TRUE
        )
        if (!inherits(fit_bayes_final, "try-error")) {
          newdat_sim <- data.frame(group = c(0, 1))
          mu_draws_sim <- try(posterior_linpred(fit_bayes_final, newdata = newdat_sim, transform = TRUE, re_formula = NA), silent = TRUE)
          alpha_draws_sim <- try(posterior_linpred(fit_bayes_final, dpar = "shape", transform = TRUE), silent = TRUE)
          
          if (!inherits(mu_draws_sim, "try-error") && !inherits(alpha_draws_sim, "try-error")) {
            mu0_draws <- mu_draws_sim[, 1]
            mu1_draws <- mu_draws_sim[, 2]
            alpha_vec <- alpha_draws_sim[, 1]
            
            rmst0_draws <- rmst_weibull(mu = mu0_draws, alpha = alpha_vec, tau = tau_final)
            rmst1_draws <- rmst_weibull(mu = mu1_draws, alpha = alpha_vec, tau = tau_final)
            
            PP_final <- mean(rmst1_draws > rmst0_draws)
          }
        }
      }
      
      fit_ok <- length(interim_PP) == length(interim_fracs) &&
        all(is.finite(interim_PP)) &&
        length(PP_final) == 1L &&
        is.finite(PP_final)
      list(
        ok = fit_ok,
        failure_reason = if (fit_ok) {
          ""
        } else {
          "mcmc_fit_or_posterior_extraction_failed"
        },
        PP = interim_PP,
        PP_final = PP_final
      )
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    saveRDS(stats_list, cache_file)
  }
  
  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes

  iv_a05_a06 <- rmst_iv_write_a05_a06_outputs(
    stats_list = stats_list,
    label = label,
    c_primes = c_primes,
    N.max = N.max,
    interim_fracs = interim_fracs,
    PP_final_threshold = PP_final_threshold,
    lambda_prime_bayes = lambda_prime_bayes,
    output_dir = iv_output_dir
  )
  valid_stats <- stats_list[iv_a05_a06$valid_indices]
  if (length(valid_stats) == 0L) {
    stop(sprintf("场景 %s 没有有效的RMST拟合结果。", label))
  }
  
  stopped_real <- sapply(valid_stats, function(stats) {
    any(stats$PP <= c_primes, na.rm = TRUE)
  })
  ess_real <- sapply(valid_stats, function(stats) {
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      return(N.max * interim_fracs[stop_idx[1]])
    } else {
      return(N.max)
    }
  })
  power_real <- mean(sapply(valid_stats, function(stats) {
    stop_val <- any(stats$PP <= c_primes, na.rm = TRUE)
    !stop_val && !is.na(stats$PP_final) && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf(
    "  - 有效模拟次数             : %d/%d\n",
    iv_a05_a06$n_valid,
    iv_a05_a06$n_total
  ))
  cat(sprintf(
    "  - 拟合失败次数             : %d\n",
    iv_a05_a06$n_failed
  ))
  cat("  - A05输出                  : 概率、MCSE、Wilson 95% Monte Carlo区间\n")
  cat("  - A06输出                  : final no-go直接计数及PET+go+final no-go恒等核对\n")
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}
