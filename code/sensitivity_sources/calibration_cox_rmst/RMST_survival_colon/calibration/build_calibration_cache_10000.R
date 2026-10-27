source(local({
  .file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .starts <- c(if (length(.file_arg)) dirname(normalizePath(sub("^--file=", "", .file_arg[1]), mustWork = TRUE)) else character(), getwd())
  .shared <- NULL
  for (.start in unique(.starts)) {
    .dir <- normalizePath(.start, mustWork = TRUE)
    repeat {
      .candidate <- file.path(.dir, "code", "common", "source_path.R")
      if (file.exists(.candidate)) { .shared <- .candidate; break }
      .parent <- dirname(.dir)
      if (identical(.parent, .dir)) break
      .dir <- .parent
    }
    if (!is.null(.shared)) break
  }
  if (is.null(.shared)) stop("Cannot locate code/common/source_path.R")
  .shared
}), local = TRUE)
source(project_source_path("analysis_helpers.R"), local = TRUE)

# CALIBRATION CACHE BUILD ONLY
# Method: RMST; case: survival_colon
# Formal run: 10,000 newly generated attempts per scenario.
# This script generates posterior-probability caches for parameter selection.
# It does not select, lock, or independently validate decision parameters.
# Base MCMC: 4 chains, 2,000 iterations, 1,000 warmup,
#            adapt_delta=0.95, max_treedepth=12.
# Retry: same simulated dataset; 4 chains, 4,000 iterations, 2,000 warmup,
#        adapt_delta=0.99, max_treedepth=15.
# Failed attempts remain explicit and are never encoded as PP=0.

cal_version <- "CAL10000_V2_STABLE_MCMC_20260730"
cal_project_id <- "RMST_survival_colon"
cal_method <- "RMST"
cal_case <- "survival_colon"
cal_data_seed_base <- 824001L
cal_model_seed_base <- 924001L

cal_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(cal_file_arg) != 1L) {
  stop("Run this calibration-cache file with Rscript.")
}
cal_script_file <- normalizePath(sub("^--file=", "", cal_file_arg), mustWork = TRUE)
cal_version_dir <- dirname(cal_script_file)
cal_project_dir <- normalizePath(file.path(cal_version_dir, ".."), mustWork = TRUE)
setwd(cal_project_dir)

cal_smoke_test <- identical(Sys.getenv("CAL_SMOKE_TEST", unset = "0"), "1")
if (cal_smoke_test) {
  cal_smoke_n <- suppressWarnings(as.integer(Sys.getenv("N_SIM_REALCASE", unset = "")))
  if (!is.finite(cal_smoke_n) || cal_smoke_n < 1L || cal_smoke_n > 20L) {
    stop("Smoke mode requires N_SIM_REALCASE between 1 and 20.")
  }
  cal_run_tag <- sprintf("smoke_n%d", cal_smoke_n)
} else {
  Sys.setenv(N_SIM_REALCASE = "10000")
  cal_run_tag <- "formal_10000"
}

cal_output_dir <- file.path(cal_version_dir, "outputs", cal_run_tag)
cal_cache_dir <- file.path(cal_output_dir, "cache")
dir.create(cal_cache_dir, recursive = TRUE, showWarnings = FALSE)

source(
  project_source_path("rmst", "stable_fit.R"),
  local = TRUE
)

cal_seed_from_label <- calibration_seed_from_label

writeLines(
  c(
    paste0("version=", cal_version),
    paste0("project=", cal_project_id),
    paste0("method=", cal_method),
    paste0("case=", cal_case),
    paste0("run_tag=", cal_run_tag),
    "futility_rule=PP <= b_k",
    paste0("data_seed_base=", cal_data_seed_base),
    paste0("model_seed_base=", cal_model_seed_base),
    "replicate_output=attempt_id,fit_ok,failure_reason,PP_interims,PP_final,fit_attempts",
    "cache_purpose=parameter_calibration_only",
    "failure_rule=failed_attempts_are_explicit_and_excluded_from_valid_denominator",
    "mcmc_retry=same_simulated_dataset_with_stricter_sampler_settings",
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  file.path(cal_output_dir, "RUN_METADATA.txt")
)


# PP/lambda must be supplied from this case's calibration output.

# 1. 加载依赖并配置编译环境
# Optional compiler settings may be supplied through R_MAKEVARS_USER.
Sys.setenv(PKG_CXXFLAGS = "-DTBB_INTERFACE_NEW")
old_path <- Sys.getenv("PATH")
r_env_bin <- file.path(dirname(dirname(R.home())), "bin")
Sys.setenv(PATH = paste(r_env_bin, old_path, sep = .Platform$path.sep))

options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

library(survival)
library(brms)
library(posterior)

rmst_weibull <- monitoring_rmst_weibull

# 2. 重复次数与结果输出设置
n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!cal_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal calibration-cache generation must use exactly 10,000 replicates.")
}
results_file <- if (N.sim_realcase == 10000) "calibration_cache_build_results_10k_priorN0_1p5_explicit.txt" else "calibration_cache_build_results_priorN0_1p5_explicit.txt"
results_file <- file.path(cal_output_dir, basename(results_file))
sink(results_file, split = TRUE)

source(project_source_path("dgm", "survival_colon_shared.R"), local = TRUE)
trial_settings <- survival_colon_trial_settings()

# 4. 决策参数
# Locked before calibration-cache generation; no environment override is permitted.
N.max <- trial_settings$N.max
interim_fracs <- trial_settings$interim_fracs
recruit.int <- trial_settings$recruit.int
tau_fixed <- trial_settings$tau_fixed
bayes_chains <- 4L
bayes_iter <- 2000L
bayes_warmup <- 1000L
bayes_adapt <- 0.95
bayes_treedepth <- 12L
bayes_retry_iter <- 4000L
bayes_retry_warmup <- 2000L
bayes_retry_adapt <- 0.99
bayes_retry_treedepth <- 15L

# 打乱受试者以模拟真实随机入组顺序
# Planned entry times used by all newly generated validation replicates.
entry_times <- recruit.int * (0:(N.max - 1))

# 预编译基础贝叶斯模型 (Weibull 家族)
cat("\n正在预编译基础贝叶斯 Weibull 模型 (以极大提速后续分析)...\n")
df_dummy <- data.frame(
  time   = runif(100, 0.5, 40),
  status = sample(c(0, 1), 100, replace = TRUE),
  group  = sample(c(0, 1), 100, replace = TRUE)
)
base_model <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = weibull(),
  prior   = c(
    prior(normal(0, 1.5), class = "b", coef = "group"),
    prior(normal(log(50), 1.5), class = "Intercept"),
    prior(lognormal(0, 1), class = "shape")
  ),
  chains  = 1,
  cores   = 1,
  iter    = 10,
  warmup  = 5,
  seed    = 12345,
  refresh = 0,
  silent  = 2
)
cat("Base model precompilation completed.\n\n")

# Generate the scenario replicates.
# Simulate the configured scenarios.
cat("\n=========================================================================\n")
cat(sprintf("正在对三场景进行 %d 次自适应设计模拟特征验证 (RMST 法)...\n", N.sim_realcase))
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
source(project_source_path("dgm", "survival_colon_shared.R"), local = TRUE)


cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30



run_validation_scenario <- function(label, scenario) {
  scenario_key <- if (grepl("H0", label, fixed = TRUE)) {
    "H0"
  } else if (grepl("H1_observed", label, fixed = TRUE)) {
    "H1_observed"
  } else if (grepl("H1", label, fixed = TRUE)) {
    "H1_standard"
  } else {
    "H1_observed"
  }
  cache_file <- file.path(
    cal_cache_dir,
    sprintf(
      "calibration_stats_%s_%s.rds",
      scenario_key,
      cal_run_tag
    )
  )
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
    
    set.seed(cal_seed_from_label(label))
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
      interim_attempt <- rep(NA_character_, length(interim_fracs))
      interim_failure <- rep("", length(interim_fracs))
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
          fit_result <- rmst_fit_with_stable_retry(
            base_model = base_model,
            newdata = df_interim,
            seed = cal_model_seed_base + j + k * 10000L
          )
          interim_attempt[k] <- fit_result$attempt
          interim_failure[k] <- fit_result$failure_reason
          if (isTRUE(fit_result$ok)) {
            rmst0_draws <- rmst_weibull(
              mu = fit_result$draws$mu0,
              alpha = fit_result$draws$alpha,
              tau = tau_k
            )
            rmst1_draws <- rmst_weibull(
              mu = fit_result$draws$mu1,
              alpha = fit_result$draws$alpha,
              tau = tau_k
            )
            interim_PP[k] <- mean(rmst1_draws > rmst0_draws)
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
      tau_final <- tau_fixed[4]
      
      PP_final <- NA_real_
      final_attempt <- NA_character_
      final_failure <- ""
      if (is.finite(tau_final) && tau_final > 0) {
        fit_result <- rmst_fit_with_stable_retry(
          base_model = base_model,
          newdata = sub_df,
          seed = cal_model_seed_base + j + 900000L
        )
        final_attempt <- fit_result$attempt
        final_failure <- fit_result$failure_reason
        if (isTRUE(fit_result$ok)) {
          rmst0_draws <- rmst_weibull(
            mu = fit_result$draws$mu0,
            alpha = fit_result$draws$alpha,
            tau = tau_final
          )
          rmst1_draws <- rmst_weibull(
            mu = fit_result$draws$mu1,
            alpha = fit_result$draws$alpha,
            tau = tau_final
          )
          PP_final <- mean(rmst1_draws > rmst0_draws)
        }
      }
      
      fit_ok <- length(interim_PP) == length(interim_fracs) &&
        all(is.finite(interim_PP)) &&
        length(PP_final) == 1L &&
        is.finite(PP_final)
      failure_parts <- c(
        sprintf(
          "interim_%d:%s",
          which(nzchar(interim_failure)),
          interim_failure[nzchar(interim_failure)]
        ),
        if (nzchar(final_failure)) {
          paste0("final:", final_failure)
        } else {
          character(0)
        }
      )
      list(
        ok = fit_ok,
        failure_reason = if (fit_ok) {
          ""
        } else if (length(failure_parts) > 0L) {
          paste(failure_parts, collapse = ";")
        } else {
          "missing_or_nonfinite_posterior_probability"
        },
        PP = interim_PP,
        PP_final = PP_final,
        fit_attempts = c(interim_attempt, final_attempt)
      )
    }, mc.cores = num_cores_validate)
    
    saveRDS(stats_list, cache_file)
  }
  
  valid_index <- which(vapply(
    stats_list,
    function(x) {
      is.list(x) &&
        isTRUE(x$ok) &&
        length(x$PP) == length(interim_fracs) &&
        all(is.finite(x$PP)) &&
        length(x$PP_final) == 1L &&
        is.finite(x$PP_final)
    },
    logical(1)
  ))
  n_total <- length(stats_list)
  n_valid <- length(valid_index)
  n_failed <- n_total - n_valid
  attempt_values <- unlist(
    lapply(stats_list[valid_index], function(x) x$fit_attempts),
    use.names = FALSE
  )
  summary_row <- data.frame(
    scenario = scenario_key,
    n_attempted = n_total,
    n_valid = n_valid,
    n_failed = n_failed,
    valid_fraction = n_valid / n_total,
    base_fit_count = sum(attempt_values == "base", na.rm = TRUE),
    retry_fit_count = sum(attempt_values == "retry", na.rm = TRUE),
    created_at = format(
      Sys.time(),
      tz = "Asia/Shanghai",
      usetz = TRUE
    ),
    stringsAsFactors = FALSE
  )
  summary_file <- file.path(
    cal_output_dir,
    sprintf("calibration_cache_summary_%s.csv", scenario_key)
  )
  write.csv(summary_row, summary_file, row.names = FALSE)
  failure_file <- file.path(
    cal_output_dir,
    sprintf("calibration_failures_%s.csv", scenario_key)
  )
  failed_index <- which(!(seq_len(n_total) %in% valid_index))
  failure_rows <- data.frame(
    attempt_id = failed_index,
    scenario = rep(scenario_key, length(failed_index)),
    failure_reason = vapply(
      stats_list[failed_index],
      function(x) {
        if (is.list(x) && length(x$failure_reason) == 1L) {
          as.character(x$failure_reason)
        } else {
          "invalid_replicate_structure"
        }
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(failure_rows, failure_file, row.names = FALSE)
  cat(sprintf(
    "\n[%s] attempted=%d, valid=%d, failed=%d, valid_fraction=%.6f\n",
    scenario_key,
    n_total,
    n_valid,
    n_failed,
    n_valid / n_total
  ))

}

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 [survival_colon]", "H0")
run_validation_scenario("H1_calibration [standard effect]", "H1_calibration")
run_validation_scenario("H1_observed [Lev+5FU observed effect]", "H1_observed")

sink()
cat(sprintf(
  "\nRMST calibration-cache generation completed: %s\n",
  cal_output_dir
))
