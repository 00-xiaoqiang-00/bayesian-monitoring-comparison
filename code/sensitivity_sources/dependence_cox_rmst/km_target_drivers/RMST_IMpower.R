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

# Independent Monte Carlo validation.
# Method: RMST; case: IMpower131
# Decision parameters are hard-coded from locked candidate 3502.
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# KM-target rho sensitivity validation; no recalibration or decision-parameter selection.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_SENSITIVITY_RHO_KM_TARGET"
iv_dgm_version <- "IMpower131_bivariate_lognormal_v6_E05_locked_20260801"
iv_project_id <- "RMST_IMpower"
iv_method <- "RMST"
iv_case <- "IMpower131"
iv_candidate_id <- "3502"
iv_locked_pp <- 0.863
iv_locked_lambda <- 0.180
iv_baseline_data_seed_base <- 1821001L
iv_baseline_model_seed_base <- 1921001L

if (!is.finite(iv_locked_pp) || iv_locked_pp <= 0 || iv_locked_pp >= 1) {
  stop("Set LOCKED_PP_FINAL_THRESHOLD to the parameter locked from the E05 v6 calibration.")
}
if (!is.finite(iv_locked_lambda) || iv_locked_lambda < 0) {
  stop("Set LOCKED_LAMBDA_PRIME to the parameter locked from the E05 v6 calibration.")
}

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this independent-validation file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_project_dir <- normalizePath(file.path(iv_version_dir, ".."), mustWork = TRUE)
setwd(iv_project_dir)

iv_baseline_rho <- 0.5
iv_allowed_rho <- c(low = 0.3, high = 0.7)
iv_rho_level <- Sys.getenv("SENSITIVITY_RHO_LEVEL", unset = "")
iv_rho_raw <- Sys.getenv("SENSITIVITY_RHO", unset = "")
if (!nzchar(iv_rho_level) || !iv_rho_level %in% names(iv_allowed_rho)) {
  stop("SENSITIVITY_RHO_LEVEL is required and must be exactly 'low' or 'high'.")
}
if (!nzchar(iv_rho_raw)) {
  stop("SENSITIVITY_RHO is required.")
}
iv_sensitivity_rho <- suppressWarnings(as.numeric(iv_rho_raw))
if (length(iv_sensitivity_rho) != 1L || !is.finite(iv_sensitivity_rho) ||
    abs(iv_sensitivity_rho) >= 1) {
  stop("SENSITIVITY_RHO must be one finite numeric value with abs(rho) < 1.")
}
iv_expected_rho <- unname(iv_allowed_rho[[iv_rho_level]])
if (!isTRUE(all.equal(iv_sensitivity_rho, iv_expected_rho, tolerance = 1e-12))) {
  stop(sprintf(
    "SENSITIVITY_RHO mismatch: level '%s' requires rho=%s for %s.",
    iv_rho_level, format(iv_expected_rho, digits = 15), iv_case
  ))
}
iv_rho_safe_tag <- gsub(".", "p", format(
  iv_sensitivity_rho, scientific = FALSE, trim = TRUE, digits = 15
), fixed = TRUE)
iv_rho_safe_tag <- gsub("-", "m", iv_rho_safe_tag, fixed = TRUE)
iv_fixed_interim_schedule <- "272,408,544 (40%,60%,80%)"

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

iv_seed_offset <- 0L
# Seed bases are set by the KM-target helper after project data-generation code is sourced.

iv_output_root <- Sys.getenv("RHO_SENSITIVITY_OUTPUT_ROOT", unset = iv_version_dir)
iv_output_dir <- file.path(
  iv_output_root, "outputs", iv_method, iv_project_id, paste0("rho_", iv_rho_safe_tag), iv_run_tag
)
iv_cache_dir <- file.path(iv_output_dir, "cache")
dir.create(iv_cache_dir, recursive = TRUE, showWarnings = FALSE)

iv_write_completion_marker <- monitoring_rho_write_completion_marker

source(
  project_source_path("rmst", "validation_helpers.R"),
  local = TRUE
)



# RUN_METADATA.txt is written after KM-target parameters and seeds are locked.


# 锁定参数: PP_final_threshold = 0.863, lambda_prime_bayes = 0.180

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
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal independent validation must use exactly 10,000 replicates.")
}
results_file <- if (N.sim_realcase == 10000) "independent_validation_results_10k.txt" else "independent_validation_results.txt"
results_file <- file.path(iv_output_dir, basename(results_file))
sink(results_file, split = TRUE)

# Simulated data are generated by the case DGM.
# 4. 配置自适应设计参数
PP_final_threshold <- iv_locked_pp
lambda_prime_bayes <- iv_locked_lambda
N.max <- 680
interim_fracs <- c(0.40, 0.60, 0.80)
recruit.int <- 22 / 680
tau_fixed <- c(6, 9, 12, 18)
bayes_chains <- 4L
bayes_iter <- 2000L
bayes_warmup <- 1000L
bayes_adapt <- 0.95
bayes_treedepth <- 12L

# 打乱受试者以模拟真实随机入组顺序
# Planned entry times used by all newly generated validation replicates.
entry_times <- recruit.int * (0:(N.max - 1))

# 预编译基础贝叶斯模型（Weibull家族）
cat("\n正在预编译基础贝叶斯 Weibull 模型...\n")
df_dummy <- data.frame(
  time   = seq(0.5, 40, length.out = 100),
  status = rep(c(0L, 1L), 50),
  group  = rep(c(0L, 1L), each = 50)
)

rmst_priors <- c(
  prior(normal(0, 1.5), class = "b", coef = "group"),
  prior(normal(0, 2.5), class = "Intercept"),
  prior(lognormal(0, 1), class = "shape")
)

base_model <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = weibull(link = "log"),
  prior   = rmst_priors,
  chains  = 1,
  cores   = 1,
  iter    = 100,
  warmup  = 50,
  seed    = 12345,
  refresh = 0,
  silent  = 2
)
cat("Base model precompilation completed.\n\n")


# Simulate the configured scenarios.
cat("\n=========================================================================\n")
cat(sprintf("正在对H0/H1-standard/H1-observed三种场景进行 %d 次自适应设计模拟特征验证 (RMST 法)...\n", N.sim_realcase))
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
iv_data_generator_file <- Sys.getenv("RHO_SENSITIVITY_DATA_GENERATOR", unset = "")
if (!nzchar(iv_data_generator_file)) stop("RHO_SENSITIVITY_DATA_GENERATOR is required.")
source(normalizePath(iv_data_generator_file, mustWork = TRUE), local = TRUE)


cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 60



# 删失感知参数（通过 survreg 对数正态 AFT 模型估计，正确处理右删失）


iv_rebuild_sigma_for_rho <- monitoring_rebuild_sigma_for_rho

iv_assert_baseline_sigma <- monitoring_assert_baseline_sigma

mu.ctrl.baseline <- c(log(6.1202), log(13.4685))   # Latent progression / death log-location
mu.ctrl <- mu.ctrl.baseline
stopifnot(identical(mu.ctrl, mu.ctrl.baseline))
Sigma.ctrl.baseline <- matrix(c(0.4209, 0.3287, 0.3287, 1.0268), nrow = 2)
iv_assert_baseline_sigma(Sigma.ctrl.baseline, "IMpower control")
Sigma.ctrl <- iv_rebuild_sigma_for_rho(
  Sigma.ctrl.baseline, iv_sensitivity_rho, "IMpower control"
)

mu.trt.real <- c(log(7.5610), log(14.3431))   # H1-observed treatment scenario
Sigma.trt.real <- matrix(c(0.8207, 0.5298, 0.5298, 1.3679), nrow = 2)

# H1-standard: recurrence and death medians of 8 and 18 months.
mu.trt.H1_sens2.baseline <- c(log(8.00), log(18.00))
mu.trt.H1_sens2 <- mu.trt.H1_sens2.baseline
stopifnot(identical(mu.trt.H1_sens2, mu.trt.H1_sens2.baseline))
Sigma.trt.H1_sens2.baseline <- Sigma.ctrl.baseline
iv_assert_baseline_sigma(Sigma.trt.H1_sens2.baseline, "IMpower H1-standard")
Sigma.trt.H1_sens2 <- iv_rebuild_sigma_for_rho(
  Sigma.trt.H1_sens2.baseline, iv_sensitivity_rho, "IMpower H1-standard"
)


# Load case- and scenario-specific KM-target parameters and seed rules.
iv_rho_shared_dir <- normalizePath(file.path(iv_version_dir, "..", "..", "..", "common", "rho"), mustWork = TRUE)
Sys.setenv(RHO_SENSITIVITY_VERSION_DIR = iv_rho_shared_dir)
iv_rho_helper <- Sys.getenv("RHO_SENSITIVITY_HELPER", unset = file.path(iv_rho_shared_dir, "rho_sensitivity_parameter_helper.R"))
if (!nzchar(iv_rho_helper)) stop("RHO_SENSITIVITY_HELPER is required.")
source(iv_rho_helper, local = TRUE)
iv_master_seed <- as.integer(Sys.getenv("RHO_SENSITIVITY_MASTER_SEED", unset = "8252026"))
if (!is.finite(iv_master_seed)) stop("Invalid RHO_SENSITIVITY_MASTER_SEED.")
iv_data_seed_base <- rho_data_seed_base(iv_master_seed, iv_case, iv_rho_level)
iv_model_seed_base <- rho_model_seed_base(iv_master_seed, iv_case, iv_method, iv_rho_level)
if (iv_smoke_test) {
  iv_data_seed_base <- as.integer((iv_data_seed_base + 500000000L) %% (.Machine$integer.max - 1))
  iv_model_seed_base <- as.integer((iv_model_seed_base + 700000000L) %% (.Machine$integer.max - 1))
}
iv_formal_data_seed_base <- rho_data_seed_base(iv_master_seed, iv_case, iv_rho_level)
iv_formal_model_seed_base <- rho_model_seed_base(iv_master_seed, iv_case, iv_method, iv_rho_level)
iv_sensitivity_formal_data_seed_base <- iv_formal_data_seed_base
iv_sensitivity_formal_model_seed_base <- iv_formal_model_seed_base
iv_seed_from_label <- monitoring_rho_seed_from_label
iv_parameter_table_sha256 <- Sys.getenv("RHO_SENSITIVITY_PARAMETER_SHA256", unset = "not_supplied")
iv_parameter_table_path <- rho_parameter_file
iv_write_km_metadata <- monitoring_rho_write_metadata

iv_ctrl_km <- rho_make_impower_parameters("H0", "control", iv_rho_level)
iv_ctrl_h1_standard_km <- rho_make_impower_parameters("H1-standard", "control", iv_rho_level)
iv_ctrl_h1_observed_km <- rho_make_impower_parameters("H1-observed", "control", iv_rho_level)
iv_h1_standard_km <- rho_make_impower_parameters("H1-standard", "treatment", iv_rho_level)
iv_h1_observed_km <- rho_make_impower_parameters("H1-observed", "treatment", iv_rho_level)
mu.ctrl <- iv_ctrl_km$mu; Sigma.ctrl <- iv_ctrl_km$Sigma
mu.ctrl.H1_standard <- iv_ctrl_h1_standard_km$mu; Sigma.ctrl.H1_standard <- iv_ctrl_h1_standard_km$Sigma
mu.ctrl.H1_observed <- iv_ctrl_h1_observed_km$mu; Sigma.ctrl.H1_observed <- iv_ctrl_h1_observed_km$Sigma
mu.trt.H1_standard <- iv_h1_standard_km$mu; Sigma.trt.H1_standard <- iv_h1_standard_km$Sigma
mu.trt.H1 <- iv_h1_standard_km$mu; Sigma.trt.H1 <- iv_h1_standard_km$Sigma
mu.trt.real <- iv_h1_observed_km$mu; Sigma.trt.real <- iv_h1_observed_km$Sigma
mu.trt.H1_sens2 <- mu.trt.H1_standard; Sigma.trt.H1_sens2 <- Sigma.trt.H1_standard
iv_write_km_metadata()

run_validation_scenario <- function(label, mu_trt, Sigma_trt, mu_ctrl, Sigma_ctrl) {
  expected_cache_context <- rho_cache_context(
    c(entry = iv_script_file, dgm = iv_data_generator_file,
      helper = iv_rho_helper, parameters = iv_parameter_table_path,
      analysis_helpers = project_source_path("analysis_helpers.R"),
      source_paths = project_source_path("source_path.R")),
    case = iv_case, method = iv_method, scenario = label, n = N.sim_realcase,
    data_seed = iv_seed_from_label(label), model_seed_base = iv_model_seed_base,
    rho_level = iv_rho_level)
  cache_suffix <- if (N.sim_realcase == 10000) "10k_v2_" else "v2_"
  cache_file <- sprintf("validate_stats_%s%s.rds", cache_suffix, ifelse(grepl("H0", label), "H0", ifelse(grepl("H0", label), "H0", ifelse(grepl("H1_observed", label), "H1_observed", "H1_standard"))))
  use_cache <- FALSE
  cache_file <- file.path(iv_cache_dir, basename(cache_file))
  if (file.exists(cache_file)) {
    stats_list_temp <- readRDS(cache_file)
    if (rho_cache_is_current(stats_list_temp, expected_cache_context, N.sim_realcase)) {
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
      tau_final <- tau_fixed[4]
      
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
    
    saveRDS(rho_cache_stamp(stats_list, expected_cache_context), cache_file)
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

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 (原假设/药物无效) [校准用场景]", mu.ctrl, Sigma.ctrl, mu.ctrl, Sigma.ctrl)
run_validation_scenario("H1_standard [KM-target standard effect]", mu.trt.H1_standard, Sigma.trt.H1_standard, mu.ctrl.H1_standard, Sigma.ctrl.H1_standard)
run_validation_scenario("H1_observed [KM-target observed effect]", mu.trt.real, Sigma.trt.real, mu.ctrl.H1_observed, Sigma.ctrl.H1_observed)

sink()
cat("\nSensitivity validation completed; A05/A06 replicate and summary outputs were written.\n")

iv_write_completion_marker()
