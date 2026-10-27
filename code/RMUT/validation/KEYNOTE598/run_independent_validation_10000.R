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
  rmut_preflight(.rmut_directory, "KEYNOTE598", "validation")
  quit(status = 0L)
}
# Formal data generation uses Inversion; pairing restores this RNG kind and state.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
# Independent Monte Carlo validation.
# Method: RMUT; case: KEYNOTE-598
# Default thresholds: PP_final_threshold = 0.971; lambda_prime_bayes = 0.405.
# A verified frozen parameter manifest may override these defaults for a formal run.
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_ROWPAIR_20261010_V2"
iv_project_id <- "weight_KEYNOTE598_RMUT"
iv_method <- "RMUT"
iv_case <- "KEYNOTE-598"
iv_candidate_id <- "31204"
iv_locked_pp <- 0.971
iv_locked_lambda <- 0.405
.rmut_parameters <- rmut_resolve_parameters("KEYNOTE598", iv_candidate_id, iv_locked_pp, iv_locked_lambda)
iv_candidate_id <- .rmut_parameters$candidate_id
iv_locked_pp <- .rmut_parameters$c
iv_locked_lambda <- .rmut_parameters$gamma
iv_data_seed_base <- 1832001L
iv_model_seed_base <- 1932001L

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this independent-validation file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_project_dir <- normalizePath(file.path(iv_version_dir, ".."), mustWork = TRUE)
setwd(iv_version_dir)

iv_smoke_test <- identical(Sys.getenv("IV_SMOKE_TEST", unset = "0"), "1")
if (iv_smoke_test) {
  iv_smoke_n <- suppressWarnings(as.integer(Sys.getenv("N_SIM_REALCASE", unset = "")))
  if (!is.finite(iv_smoke_n) || iv_smoke_n < 1L || iv_smoke_n > 20L) {
    stop("Smoke mode requires N_SIM_REALCASE between 1 and 20.")
  }
  iv_run_tag <- sprintf("smoke_n%d", iv_smoke_n)
} else {
  Sys.setenv(N_SIM_REALCASE = "10000")
  iv_run_tag <- "formal_10000"
}

iv_formal_data_seed_base <- iv_data_seed_base
iv_formal_model_seed_base <- iv_model_seed_base
if (iv_smoke_test) {
  iv_data_seed_base <- iv_data_seed_base + 500000000L
  iv_model_seed_base <- iv_model_seed_base + 500000000L
}

iv_output_dir <- file.path(iv_version_dir, "outputs", iv_run_tag)
iv_cache_dir <- file.path(iv_output_dir, "cache")
dir.create(iv_cache_dir, recursive = TRUE, showWarnings = FALSE)

iv_seed_from_label <- monitoring_seed_from_label

writeLines(
  c(
    paste0("version=", iv_version),
    paste0("project=", iv_project_id),
    paste0("method=", iv_method),
    paste0("case=", iv_case),
    paste0("candidate_id=", iv_candidate_id),
    paste0("parameter_status=", .rmut_parameters$status),
    paste0("pairing_version=", rmut_stage_version),
    paste0("run_tag=", iv_run_tag),
    paste0("locked_PP_final_threshold=", format(iv_locked_pp, digits = 15)),
    paste0("locked_lambda_prime_bayes=", format(iv_locked_lambda, digits = 15)),
    "futility_rule=PP <= b_k",
    "mcmc_base=4 chains; iter 2000; warmup 1000; adapt_delta 0.95; max_treedepth 12",
    "mcmc_retry=4 chains; iter 4000; warmup 2000; adapt_delta 0.99; max_treedepth 15",
    "retry_policy=retry only the failed AFT component on the same simulated dataset",
    paste0("data_seed_base=", iv_data_seed_base),
    paste0("model_seed_base=", iv_model_seed_base),
    paste0("formal_data_seed_base=", iv_formal_data_seed_base),
    paste0("formal_model_seed_base=", iv_formal_model_seed_base),
    "seed_policy=smoke test uses +500000000 offset; formal streams remain untouched",
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  file.path(iv_output_dir, "RUN_METADATA.txt")
)

# RMUT posterior-decision calculation.
# KEYNOTE598 案例的独立模拟验证。

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
sink(file.path(iv_output_dir, "KEYNOTE598_RMUT_validation_summary.txt"),
  split = TRUE,
  append = TRUE
)

source(project_source_path("dgm", "KEYNOTE598_shared.R"), local = TRUE)
trial_settings <- keynote598_trial_settings(tau_fixed = c(4, 6, 9, 18))
# Generate and analyse the scenario replicates.
# 4. 配置自适应设计参数；使用本文件锁定值
# Locked decision thresholds are specified below.
# Use the thresholds resolved from the defaults or the validated parameter manifest.
PP_final_threshold <- iv_locked_pp
lambda_prime_bayes <- iv_locked_lambda
N.max <- trial_settings$N.max
interim_fracs <- trial_settings$interim_fracs
interim_counts <- c(228L, 340L, 454L)
stopifnot(
  length(interim_counts) == length(interim_fracs),
  all(interim_counts > 0L),
  all(interim_counts < N.max)
)
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
  prior(normal(log(25), 1.5), class = "Intercept"),
  prior(lognormal(log(0.8), 0.75), class = "sigma")
)
prior_composite <- c(
  prior(normal(0, 1.5), class = "b", coef = "group"),
  prior(normal(log(10), 1.5), class = "Intercept"),
  prior(lognormal(log(0.8), 0.75), class = "sigma")
)


# Planned entry times used by all newly generated validation replicates.
entry_times <- recruit.int * (0:(N.max - 1))

# RMUT 后验计算函数
source(file.path(iv_version_dir, "..", "..", "common", "rmut_pairing.R"), local = TRUE)

# 仅使用通过收敛和有效样本量筛查的后验抽样。


source(
  file.path(iv_version_dir, "..", "..", "common", "rmut_stable_fit_helper.R"),
  local = TRUE
)

# 分别预编译两个组成模型，以使用终点特异的截距先验。
cat("\n正在预编译死亡和首次进展或死亡 AFT 模型...\n")
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
# H0, H1-standard and H1-observed scenario analyses.
cat("\n=========================================================================\n")
cat("准备对三场景进行自适应设计模拟特征验证 (三状态RMUT法)...\n")
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
source(project_source_path("dgm", "KEYNOTE598_shared.R"), local = TRUE)

cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30
if (length(num_cores_validate) != 1L || !is.finite(num_cores_validate) || num_cores_validate < 1L) stop("NUM_CORES must be a positive integer")

n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal independent validation must use exactly 10,000 replicates.")
}
analysis_version <- "KEYNOTE598_v2_RMUT3STATE_v4_centered_weakpriors_twoAFT_ROWPAIR_20261010_V2"
cache_signature <- sprintf(
  "%s_w%.2f_tau%s",
  analysis_version, w_R, paste(tau_fixed, collapse = "-")
)

run_validation_scenario <- function(label, scenario) {
  scenario_id <- scenario
  rmut_scenario_code <- if (grepl("H0", scenario_id)) 1L else if (grepl("observed|Realcase", scenario_id)) 3L else 2L
  rmut_context <- rmut_cache_context(iv_version_dir, "KEYNOTE598", "validation", label, N.sim_realcase, iv_seed_from_label(label), iv_model_seed_base, w_R, tau_fixed)
  cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_id, cache_signature)
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
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
    
    set.seed(iv_seed_from_label(label))
    sim_array <- data.simulation(
      N.sim = N.sim_realcase, N.max = N.max, 
      cens_upper = trial_settings$cens_upper,
      scenario = scenario
    )
    
    stats_list <- mclapply(seq_len(N.sim_realcase), function(j) {
      df_j <- as.data.frame(sim_array[, , j])
      
      interim_PP <- rep(NA_real_, length(interim_fracs))
      pairing_records <- list()
      original_censor_t <- df_j$censor_t
      
      for (k in seq_along(interim_fracs)) {
    n_k <- interim_counts[k]
        Time.current <- recruit.int * (n_k - 1)
        obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
        censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
        
        death_event_t <- df_j$death_t[1:n_k]
        death_time    <- pmin(death_event_t, censor_t_interim)
        death_status  <- as.integer(death_event_t <= censor_t_interim)
        
        composite_event_t <- pmin(
          df_j$progression_t[1:n_k],
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
          fit_death_interim <- rmut_fit_with_stable_retry(
            base_model_death,
            df_death_interim,
            iv_model_seed_base + j
          )
          fit_composite_interim <- rmut_fit_with_stable_retry(
            base_model_composite,
            df_composite_interim,
            iv_model_seed_base + j + 100000L
          )
          
          pairing_records[[k]] <- rmut_draw_pair_record(NULL, NULL, rmut_pair_seed(j, k, rmut_scenario_code, 2L), sprintf("interim_%d", k), fit_death_interim, fit_composite_interim)
          if (!inherits(fit_death_interim, "try-error") && !inherits(fit_composite_interim, "try-error")) {
            draws_d <- extract_valid_aft_draws(fit_death_interim)
            draws_c <- extract_valid_aft_draws(fit_composite_interim)
            
            if (!is.null(draws_d) && !is.null(draws_c) &&
                all(c("b_Intercept", "b_group", "sigma") %in% colnames(draws_d)) &&
                all(c("b_Intercept", "b_group", "sigma") %in% colnames(draws_c))) {
              
              pairing_records[[k]] <- rmut_draw_pair_record(draws_d, draws_c, rmut_pair_seed(j, k, rmut_scenario_code, 2L), sprintf("interim_%d", k), fit_death_interim, fit_composite_interim)
              rmut_ctrl <- compute_rmut_samples(
                mu_D    = draws_d$b_Intercept,
                sigma_D = draws_d$sigma,
                mu_C    = draws_c$b_Intercept,
                sigma_C = draws_c$sigma,
                tau     = tau_k,
                w_R     = w_R,
                pair_seed = rmut_pair_seed(j, k, rmut_scenario_code, 2L)
              )
              
              rmut_trt <- compute_rmut_samples(
                mu_D    = draws_d$b_Intercept + draws_d$b_group,
                sigma_D = draws_d$sigma,
                mu_C    = draws_c$b_Intercept + draws_c$b_group,
                sigma_C = draws_c$sigma,
                tau     = tau_k,
                w_R     = w_R,
                pair_seed = rmut_pair_seed(j, k, rmut_scenario_code, 2L)
              )
              
              interim_PP[k] <- mean((rmut_trt - rmut_ctrl) > 0)
            }
          }
        }
      }
      
      final_censor_t <- pmin(
        original_censor_t,
        keynote598_admin_followup(
          entry_times, trial_settings$final_analysis_time
        )
      )
      death_event_t <- df_j$death_t
      death_time    <- pmin(death_event_t, final_censor_t)
      death_status  <- as.integer(death_event_t <= final_censor_t)
      
      composite_event_t <- pmin(df_j$progression_t, df_j$death_t)
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
        fit_death_final <- rmut_fit_with_stable_retry(
          base_model_death,
          df_death_final,
          iv_model_seed_base + j
        )
        fit_composite_final <- rmut_fit_with_stable_retry(
          base_model_composite,
          df_composite_final,
          iv_model_seed_base + j + 100000L
        )
        
        pairing_records[[length(tau_fixed)]] <- rmut_draw_pair_record(NULL, NULL, rmut_pair_seed(j, length(tau_fixed), rmut_scenario_code, 2L), "final", fit_death_final, fit_composite_final)
        if (!inherits(fit_death_final, "try-error") && !inherits(fit_composite_final, "try-error")) {
          draws_d <- extract_valid_aft_draws(fit_death_final)
          draws_c <- extract_valid_aft_draws(fit_composite_final)
          
          if (!is.null(draws_d) && !is.null(draws_c) &&
              all(c("b_Intercept", "b_group", "sigma") %in% colnames(draws_d)) &&
              all(c("b_Intercept", "b_group", "sigma") %in% colnames(draws_c))) {
            
            pairing_records[[length(tau_fixed)]] <- rmut_draw_pair_record(draws_d, draws_c, rmut_pair_seed(j, length(tau_fixed), rmut_scenario_code, 2L), "final", fit_death_final, fit_composite_final)
              rmut_ctrl <- compute_rmut_samples(
              mu_D    = draws_d$b_Intercept,
              sigma_D = draws_d$sigma,
              mu_C    = draws_c$b_Intercept,
              sigma_C = draws_c$sigma,
              tau     = tau_final,
              w_R     = w_R,
                pair_seed = rmut_pair_seed(j, length(tau_fixed), rmut_scenario_code, 2L)
            )
            
            rmut_trt <- compute_rmut_samples(
              mu_D    = draws_d$b_Intercept + draws_d$b_group,
              sigma_D = draws_d$sigma,
              mu_C    = draws_c$b_Intercept + draws_c$b_group,
              sigma_C = draws_c$sigma,
              tau     = tau_final,
              w_R     = w_R,
                pair_seed = rmut_pair_seed(j, length(tau_fixed), rmut_scenario_code, 2L)
            )
            
            PP_final <- mean((rmut_trt - rmut_ctrl) > 0)
          }
        }
      }
      
      rmut_validation_result(j, interim_PP, PP_final, pairing_records)
    }, mc.cores = num_cores_validate, mc.preschedule = FALSE)
    
    stats_list <- rmut_stamp_cache(stats_list, rmut_context)
    saveRDS(stats_list, cache_file)
  }
  
  rmut_write_validation_audit(stats_list, iv_output_dir, scenario_id, PP_final_threshold * (interim_counts/N.max)^lambda_prime_bayes, PP_final_threshold, interim_counts, N.max)

  valid_stats <- Filter(function(stats) {
    is.list(stats) && length(stats$PP) == length(interim_fracs) &&
      all(is.finite(stats$PP)) && is.finite(stats$PP_final)
  }, stats_list)
  if (length(valid_stats) == 0L) {
    stop(sprintf("场景 %s 没有有效的MCMC模拟结果。", label))
  }
  if (length(valid_stats) < length(stats_list)) {
    cat(sprintf("警告：%s 排除了 %d 次MCMC失败的模拟。\n",
                label, length(stats_list) - length(valid_stats)))
  }

  c_primes <- PP_final_threshold * (interim_counts / N.max)^lambda_prime_bayes
  
  stopped_real <- sapply(valid_stats, function(stats) {
    any(stats$PP <= c_primes)
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
    stop_val <- any(stats$PP <= c_primes)
    !stop_val && stats$PP_final > PP_final_threshold
  }))
  
  cat(sprintf("\nScenario: %s:\n", label))
  cat(sprintf("  - 有效模拟次数             : %d/%d\n", length(valid_stats), length(stats_list)))
  cat(sprintf("  - 预期早停率 (PET)     : %.1f%%\n", mean(stopped_real) * 100))
  cat(sprintf("  - 平均期望样本量 (ESS) : %.1f 人 (占最大样本量的 %.1f%%)\n", mean(ess_real), (mean(ess_real) / N.max) * 100))
  cat(sprintf("  - 条件于有效拟合的功效 / 一类错误率 : %.1f%%\n", power_real * 100))
  cat("=========================================================================\n")
}

# Evaluate H0, H1-standard and H1-observed using the configured thresholds.
run_validation_scenario("H0 [KEYNOTE598]", "H0")
run_validation_scenario("H1_calibration [standard effect]", "H1_calibration")
run_validation_scenario("H1_observed [pembrolizumab + ipilimumab observed effect]", "H1_observed")

sink()
