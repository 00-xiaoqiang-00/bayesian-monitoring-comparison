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

# Formal RMUT analysis with a locked configuration signature.
.rmut_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.rmut_file_arg) != 1L) stop("Run with Rscript")
.rmut_entry_file <- normalizePath(sub("^--file=", "", .rmut_file_arg), mustWork = TRUE)
.rmut_directory <- dirname(.rmut_entry_file)
source(file.path(.rmut_directory, "..", "..", "common", "production_runtime.R"), local = TRUE)
source(file.path(.rmut_directory, "..", "..", "common", "rmut_pairing.R"), local = TRUE)
if ("--preflight-only" %in% commandArgs(trailingOnly = TRUE)) {
  rmut_preflight(.rmut_directory, "survival_colon", "calibration")
  quit(status = 0L)
}
# Formal data generation uses Inversion; pairing restores this RNG kind and state.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
# CALIBRATION CACHE GENERATION ONLY
# Method: RMUT; case: survival_colon
# Reference parameters for cache-build diagnostics only: PP_final_threshold = 0.885;
#                             lambda_prime_bayes = 0.00
# Formal run: 10,000 newly generated calibration replicates.
# Futility rule: stop at interim look k when PP <= b_k.
# This file generates calibration PP caches and never selects decision parameters.
# Runtime outputs use the configured output directory.

cal_version <- "CAL10000_ROWPAIR_20261010_V2"
cal_project_id <- "weight_survival_colon_RMUT"
cal_method <- "RMUT"
cal_case <- "survival_colon"
cal_reference_pp <- 0.885
cal_reference_lambda <- 0.00
cal_data_seed_base <- 634001L
cal_model_seed_base <- 734001L

if ("--syntax-check" %in% commandArgs(trailingOnly = TRUE)) {
  cat("SYNTAX_OK\n")
  quit(status = 0L)
}

cal_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(cal_file_arg) != 1L) {
  stop("Run this calibration-cache file with Rscript.")
}
cal_script_file <- normalizePath(sub("^--file=", "", cal_file_arg), mustWork = TRUE)
cal_version_dir <- dirname(cal_script_file)
cal_project_dir <- normalizePath(file.path(cal_version_dir, ".."), mustWork = TRUE)
setwd(cal_version_dir)

cal_smoke_test <- identical(Sys.getenv("CAL_SMOKE_TEST", unset = "0"), "1")
if (cal_smoke_test) {
  cal_smoke_n <- suppressWarnings(as.integer(Sys.getenv("N_SIM_REALCASE", unset = "")))
  if (!is.finite(cal_smoke_n) || cal_smoke_n < 1L || cal_smoke_n > 20L) {
    stop("Smoke mode requires N_SIM_REALCASE between 1 and 20.")
  }
  cal_smoke_suffix <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    Sys.getenv("CAL_SMOKE_TAG_SUFFIX", unset = "")
  )
  cal_run_tag <- sprintf(
    "smoke_n%d%s",
    cal_smoke_n,
    if (nzchar(cal_smoke_suffix)) paste0("_", cal_smoke_suffix) else ""
  )
} else {
  Sys.setenv(N_SIM_REALCASE = "10000")
  cal_run_tag <- "calibration_10000"
}

cal_output_dir <- file.path(cal_version_dir, "outputs", cal_run_tag)
cal_cache_dir <- file.path(cal_output_dir, "cache")
dir.create(cal_cache_dir, recursive = TRUE, showWarnings = FALSE)

cal_seed_from_label <- calibration_seed_from_label

writeLines(
  c(
    paste0("version=", cal_version),
    paste0("project=", cal_project_id),
    paste0("method=", cal_method),
    paste0("case=", cal_case),
    paste0("run_tag=", cal_run_tag),
    "purpose=parameter_calibration_cache_generation",
    "decision_parameters_selected=false",
    "reference_metrics_only=true",
    paste0("reference_PP_final_threshold=", format(cal_reference_pp, digits = 15)),
    paste0("reference_lambda_prime_bayes=", format(cal_reference_lambda, digits = 15)),
    "futility_rule=PP <= b_k",
    "mcmc_base=4 chains; iter 2000; warmup 1000; adapt_delta 0.95; max_treedepth 12",
    "mcmc_retry=4 chains; iter 4000; warmup 2000; adapt_delta 0.99; max_treedepth 15",
    "retry_policy=retry only the failed AFT component on the same simulated dataset",
    "replicate_output=attempt_id,fit_ok,failure_reason,fit_attempts,PP_interims,PP_final",
    "scenario_outputs=replicate_csv,failure_csv,summary_csv,cache_rds",
    "failure_rule=failed fits remain explicit and are excluded from the valid denominator",
    paste0("data_seed_base=", cal_data_seed_base),
    paste0("model_seed_base=", cal_model_seed_base),
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  file.path(cal_output_dir, "RUN_METADATA.txt")
)

# RMUT posterior-decision calculation.

# Reference thresholds are used only for diagnostic summaries of calibration paths.

# 1. 加载依赖并配置编译环境
.rmut_makevars <- Sys.getenv("RMUT_MAKEVARS", unset = "")
if (file.exists(.rmut_makevars)) Sys.setenv(R_MAKEVARS_USER = .rmut_makevars)
Sys.setenv(PKG_CXXFLAGS = "-DTBB_INTERFACE_NEW")
old_path <- Sys.getenv("PATH")
r_env_bin <- file.path(dirname(dirname(R.home())), "bin")
Sys.setenv(PATH = paste(r_env_bin, old_path, sep = .Platform$path.sep))

options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

library(survival)
library(brms)
library(posterior)

# 2. 分析配置
sink(file.path(cal_output_dir, "calibration_cache_build_results.txt"), split = TRUE, append = TRUE)

source(project_source_path("dgm", "survival_colon_shared.R"), local = TRUE)
trial_settings <- survival_colon_trial_settings()

# Calibration reference thresholds and analysis settings.
# Reference thresholds for calibration-path diagnostic summaries.
# Reference values are used only to print diagnostics; the saved PP cache is decision-parameter independent.
PP_final_threshold <- cal_reference_pp
lambda_prime_bayes <- cal_reference_lambda
N.max <- trial_settings$N.max
interim_fracs <- trial_settings$interim_fracs
recruit.int <- trial_settings$recruit.int
bayes_iter <- 2000
bayes_warmup <- 1000
bayes_chains <- 4L
bayes_adapt <- 0.95
bayes_treedepth <- 12L
bayes_retry_iter <- 4000L
bayes_retry_warmup <- 2000L
bayes_retry_adapt <- 0.99
bayes_retry_treedepth <- 15L
w_R          <- as.numeric(Sys.getenv("RMUT_WEIGHT", unset = "0.5"))
tau_fixed    <- trial_settings$tau_fixed
stopifnot(length(tau_fixed) == length(interim_fracs) + 1L)
stopifnot(w_R >= 0, w_R <= 1, all(diff(tau_fixed) > 0))

prior_death <- c(
  prior(normal(0, 1.5), class = "b", coef = "group"),
  prior(normal(log(60), 1.5), class = "Intercept"),
  prior(lognormal(log(0.8), 0.75), class = "sigma")
)
prior_composite <- c(
  prior(normal(0, 1.5), class = "b", coef = "group"),
  prior(normal(log(30), 1.5), class = "Intercept"),
  prior(lognormal(log(0.8), 0.75), class = "sigma")
)


# Planned entry times used by all newly generated calibration replicates.
entry_times <- recruit.int * (0:(N.max - 1))

# RMUT 后验计算函数
source(file.path(cal_version_dir, "..", "..", "common", "rmut_pairing.R"), local = TRUE)

# 仅使用通过收敛和有效样本量筛查的后验抽样。


source(
  file.path(cal_version_dir, "..", "..", "common", "rmut_calibration_fit_helper.R"),
  local = TRUE
)

# 分别预编译两个组成模型，以使用终点特异的截距先验。
cat("\n正在预编译死亡和首次复发或死亡 AFT 模型...\n")
df_dummy <- data.frame(
  time   = runif(100, 0.5, 40),
  status = sample(c(0, 1), 100, replace = TRUE),
  group  = sample(c(0, 1), 100, replace = TRUE)
)
base_model_death <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = lognormal(),
  prior   = prior_death,
  chains  = 1,
  cores   = 1,
  iter    = 10,
  warmup  = 5,
  seed    = 12345,
  refresh = 0,
  silent  = 2
)
base_model_composite <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = lognormal(),
  prior   = prior_composite,
  chains  = 1,
  cores   = 1,
  iter    = 10,
  warmup  = 5,
  seed    = 22345,
  refresh = 0,
  silent  = 2
)
cat("Component-model precompilation completed.\n\n")

# Generate the scenario replicates.
# Posterior-probability paths for the three configured scenarios.
cat("\n=========================================================================\n")
cat("准备对三场景进行自适应设计模拟特征验证 (三状态RMUT法)...\n")
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
source(project_source_path("dgm", "survival_colon_shared.R"), local = TRUE)

cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30
if (length(num_cores_validate) != 1L || !is.finite(num_cores_validate) || num_cores_validate < 1L) stop("NUM_CORES must be a positive integer")

n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!cal_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal calibration must use exactly 10,000 replicates.")
}
analysis_version <- "survival_colon_v1_RMUT3STATE_v3_weakpriors_twoAFT_ROWPAIR_20261010_V2"
cache_signature <- sprintf(
  "%s_w%.2f_tau%s",
  analysis_version, w_R, paste(tau_fixed, collapse = "-")
)

run_calibration_scenario <- function(label, scenario) {
  scenario_id <- scenario
  rmut_scenario_code <- if (grepl("H0", scenario_id)) 1L else if (grepl("observed|Realcase", scenario_id)) 3L else 2L
  rmut_context <- rmut_cache_context(cal_version_dir, "survival_colon", "calibration", label, N.sim_realcase, cal_seed_from_label(label), cal_model_seed_base, w_R, tau_fixed)
  cache_file <- sprintf("calib_stats_%s_%s.rds", scenario_id, cache_signature)
  use_cache <- FALSE
  cache_file <- file.path(cal_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (rmut_corrected_cache_is_current(stats_list_temp, N.sim_realcase, length(interim_fracs), rmut_context)) {
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
      
      interim_PP <- rep(NA_real_, length(interim_fracs))
      fit_records <- list()
      original_censor_t <- df_j$censor_t
      
      for (k in seq_along(interim_fracs)) {
        n_k <- floor(N.max * interim_fracs[k])
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        
        death_event_t <- df_j$death_t[1:n_k]
        death_time    <- pmin(death_event_t, censor_t_interim)
        death_status  <- as.integer(death_event_t <= censor_t_interim)
        
        composite_event_t <- pmin(
          df_j$recurrence_t[1:n_k],
          df_j$death_t[1:n_k]
        )
        composite_time   <- pmin(composite_event_t, censor_t_interim)
        composite_status <- as.integer(composite_event_t <= censor_t_interim)
        
        valid_rows <- censor_t_interim > 0
        
        df_death_interim <- data.frame(
          time = death_time[valid_rows],
          status = death_status[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        
        df_composite_interim <- data.frame(
          time = composite_time[valid_rows],
          status = composite_status[valid_rows],
          group = df_j$group[1:n_k][valid_rows]
        )
        
        tau_k <- tau_fixed[k]
        
        if (is.finite(tau_k) && tau_k > 0) {
          death_result <- rmut_fit_with_stable_retry(
            base_model_death,
            df_death_interim,
            cal_model_seed_base + j + 1000L * k
          )
          composite_result <- rmut_fit_with_stable_retry(
            base_model_composite,
            df_composite_interim,
            cal_model_seed_base + j + 100000L + 1000L * k
          )
          fit_records[[length(fit_records) + 1L]] <- rmut_fit_record(
            death_result,
            sprintf("interim_%d", k),
            "death"
          )
          fit_records[[length(fit_records) + 1L]] <- rmut_fit_record(
            composite_result,
            sprintf("interim_%d", k),
            "composite"
          )
          interim_PP[k] <- rmut_pp_from_fit_results(
            death_result,
            composite_result,
            tau_k,
            w_R,
            pair_seed = rmut_pair_seed(j, k, rmut_scenario_code, 4L)
          )
        }
      }
      
      final_censor_t <- pmin(
        original_censor_t,
        survival_colon_admin_followup(
          entry_times, trial_settings$final_analysis_time
        )
      )
      death_event_t <- df_j$death_t
      death_time    <- pmin(death_event_t, final_censor_t)
      death_status  <- as.integer(death_event_t <= final_censor_t)
      
      composite_event_t <- pmin(df_j$recurrence_t, df_j$death_t)
      composite_time   <- pmin(composite_event_t, final_censor_t)
      composite_status <- as.integer(composite_event_t <= final_censor_t)
      
      df_death_final <- data.frame(
        time = death_time,
        status = death_status,
        group = df_j$group
      )
      
      df_composite_final <- data.frame(
        time = composite_time,
        status = composite_status,
        group = df_j$group
      )
      
      tau_final <- tau_fixed[4]
      
      PP_final <- NA_real_
      if (is.finite(tau_final) && tau_final > 0) {
        death_result <- rmut_fit_with_stable_retry(
          base_model_death,
          df_death_final,
          cal_model_seed_base + j + 9000L
        )
        composite_result <- rmut_fit_with_stable_retry(
          base_model_composite,
          df_composite_final,
          cal_model_seed_base + j + 109000L
        )
        fit_records[[length(fit_records) + 1L]] <- rmut_fit_record(
          death_result,
          "final",
          "death"
        )
        fit_records[[length(fit_records) + 1L]] <- rmut_fit_record(
          composite_result,
          "final",
          "composite"
        )
        PP_final <- rmut_pp_from_fit_results(
          death_result,
          composite_result,
          tau_final,
          w_R,
            pair_seed = rmut_pair_seed(j, length(tau_fixed), rmut_scenario_code, 4L)
        )
      }
      
      rmut_attach_calibration_pairing(rmut_make_replicate_result(
        attempt_id = j,
        interim_pp = interim_PP,
        pp_final = PP_final,
        fit_records = fit_records
      ), 4L, rmut_scenario_code)
    }, mc.cores = num_cores_validate)
    
    stats_list <- rmut_stamp_cache(stats_list, rmut_context)
    saveRDS(stats_list, cache_file)
  }
  
  valid_stats <- Filter(function(stats) {
    is.list(stats) && isTRUE(stats$fit_ok)
  }, stats_list)
  if (length(valid_stats) == 0L) {
    rmut_write_scenario_audit(
      stats_list,
      cal_output_dir,
      scenario_id,
      label,
      list(PET = NA_real_, ESS = NA_real_, success_rate = NA_real_)
    )
    stop(sprintf("场景 %s 没有有效的MCMC模拟结果。", label))
  }
  if (length(valid_stats) < length(stats_list)) {
    cat(sprintf("警告：%s 排除了 %d 次MCMC失败的模拟。\n",
                label, length(stats_list) - length(valid_stats)))
  }

  c_primes <- PP_final_threshold * (floor(N.max * interim_fracs) / N.max)^lambda_prime_bayes
  
  stopped_real <- sapply(valid_stats, function(stats) {
    any(stats$PP <= c_primes)
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
    stop_val <- any(stats$PP <= c_primes)
    !stop_val && stats$PP_final > PP_final_threshold
  }))

  rmut_write_scenario_audit(
    stats_list = stats_list,
    output_dir = cal_output_dir,
    scenario_id = scenario_id,
    scenario_label = label,
    metrics = list(
      PET = mean(stopped_real),
      ESS = mean(ess_real),
      success_rate = power_real
    )
  )
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf("  - 有效模拟次数             : %d/%d\n", length(valid_stats), length(stats_list)))
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 检验功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

# Generate H0, H1-standard and H1-observed PP paths; only H0/H1-standard enter parameter selection.
run_calibration_scenario("H0 [survival_colon]", "H0")
run_calibration_scenario("H1_calibration [standard effect]", "H1_calibration")
run_calibration_scenario("H1_observed [Lev+5FU observed effect]", "H1_observed")

sink()
