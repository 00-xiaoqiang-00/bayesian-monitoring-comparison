source(local({
  .file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  .script <- if (length(.file_arg)) sub("^--file=", "", .file_arg[1]) else ""
  .starts <- c(if (nzchar(.script) && file.exists(.script)) dirname(normalizePath(.script, mustWork = TRUE)) else character(), getwd())
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
# Method: Cox; case: KEYNOTE-598
# Locked decision parameters: PP_final_threshold = 0.851;
#                             lambda_prime_bayes = 0.175
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# KM-target rho sensitivity validation; no recalibration or decision-parameter selection.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_SENSITIVITY_RHO_KM_TARGET"
iv_project_id <- "cox_KEYNOTE598"
iv_method <- "Cox"
iv_case <- "KEYNOTE-598"
iv_candidate_id <- "6737"
iv_locked_pp <- 0.851
iv_locked_lambda <- 0.175
iv_baseline_data_seed_base <- 1812001L
iv_baseline_model_seed_base <- 1912001L

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
iv_fixed_interim_schedule <- "228,340,454 (locked exact KEYNOTE counts for 40/60/80 looks)"

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



# RUN_METADATA.txt is written after KM-target parameters and seeds are locked.


# Locked candidate 6737: PP = 0.851, lambda = 0.175.

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

# 2. 分析配置
sink(file.path(iv_output_dir, "independent_validation_results.txt"), split = TRUE)

iv_data_generator_file <- Sys.getenv("RHO_SENSITIVITY_DATA_GENERATOR", unset = "")
if (!nzchar(iv_data_generator_file)) stop("RHO_SENSITIVITY_DATA_GENERATOR is required.")
source(normalizePath(iv_data_generator_file, mustWork = TRUE), local = TRUE)
# Generate and analyse the scenario replicates.
trial_settings <- keynote598_trial_settings(tau_fixed = c(6, 9, 12, 18))
# 4. 决策参数
# Locked before independent validation; no environment override is permitted.
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
bayes_warmup <- 1000 # 每链保留 1000；4 链共 4000 个 post-warmup draws
cache_tag <- "km_target_rho_sensitivity_commonCutoff41_iter2000_warmup1000"

get_valid_group_draws <- monitoring_valid_group_draws

# 按照随机化临床试验的实际情况，打乱受试者入组顺序
# 随机打乱分组后用于中期分析。
# Planned entry times used by all newly generated validation replicates.
entry_times <- recruit.int * (0:(N.max - 1))

# 预编译基础贝叶斯模型
cat("\n正在预编译基础贝叶斯模型 (以极大提速后续分析)...\n")
df_dummy <- data.frame(
  time   = runif(100, 0.5, 40),
  status = sample(c(0, 1), 100, replace = TRUE),
  group  = sample(c(0, 1), 100, replace = TRUE)
)
base_model <- brm(
  formula = time | cens(1 - status) ~ group,
  data    = df_dummy,
  family  = brmsfamily("cox"),
  prior   = prior(normal(0, 1.5), class = "b", coef = "group"),
  chains  = 1,
  cores   = 1,
  iter    = 10,
  warmup  = 5,
  seed    = 12345,
  refresh = 0,
  silent  = 2
)
cat("Base model precompilation completed.\n\n")

# Run independent validation with the case data generator.
# Simulate the configured scenarios.
cat("\n=========================================================================\n")
cat("正在对H0/H1-standard/H1-observed三种场景进行自适应设计模拟验证 (Cox 法)...\n")
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
iv_data_generator_file <- Sys.getenv("RHO_SENSITIVITY_DATA_GENERATOR", unset = "")
if (!nzchar(iv_data_generator_file)) stop("RHO_SENSITIVITY_DATA_GENERATOR is required.")
source(normalizePath(iv_data_generator_file, mustWork = TRUE), local = TRUE)

cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30

n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal independent validation must use exactly 10,000 replicates.")
}
# 删失感知参数。第一维必须是进展时间，不是 PFS。


iv_baseline_parameters <- keynote598_generator_parameters()
iv_sensitivity_parameters <- iv_baseline_parameters
iv_keynote_marginals <- c("pi", "mu_R", "sigma_R", "mu_D", "sigma_D")
for (iv_arm_name in c("control", "H1_calibration")) {
  iv_baseline_arm <- iv_baseline_parameters[[iv_arm_name]]
  if (!isTRUE(all.equal(iv_baseline_arm$rho, iv_baseline_rho, tolerance = 1e-12))) {
    stop(iv_arm_name, " does not use the declared locked baseline rho.")
  }
  for (iv_k in seq_along(iv_baseline_arm$Sigma)) {
    iv_expected_baseline <- keynote598_component_sigma(
      iv_baseline_arm$sigma_R[iv_k], iv_baseline_arm$sigma_D[iv_k],
      iv_baseline_rho
    )
    if (!isTRUE(all.equal(iv_baseline_arm$Sigma[[iv_k]],
                          iv_expected_baseline, tolerance = 1e-12))) {
      stop(iv_arm_name, " component ", iv_k,
           " does not reconstruct the locked baseline covariance.")
    }
  }
  iv_sensitivity_parameters[[iv_arm_name]]$rho <- iv_sensitivity_rho
  iv_sensitivity_parameters[[iv_arm_name]]$Sigma <- lapply(
    seq_along(iv_baseline_arm$Sigma),
    function(iv_k) keynote598_component_sigma(
      iv_baseline_arm$sigma_R[iv_k], iv_baseline_arm$sigma_D[iv_k],
      iv_sensitivity_rho
    )
  )
  for (iv_name in iv_keynote_marginals) {
    if (!isTRUE(all.equal(iv_sensitivity_parameters[[iv_arm_name]][[iv_name]],
                          iv_baseline_arm[[iv_name]], tolerance = 0))) {
      stop(iv_arm_name, " marginal parameter changed: ", iv_name)
    }
  }
  for (iv_k in seq_along(iv_baseline_arm$Sigma)) {
    if (!identical(diag(iv_sensitivity_parameters[[iv_arm_name]]$Sigma[[iv_k]]),
                   diag(iv_baseline_arm$Sigma[[iv_k]]))) {
      stop(iv_arm_name, " component ", iv_k, " covariance diagonal changed.")
    }
  }
  keynote598_validate_arm(iv_sensitivity_parameters[[iv_arm_name]], iv_arm_name)
}
if (!identical(iv_sensitivity_parameters$H1_observed,
               iv_baseline_parameters$H1_observed) ||
    !identical(iv_sensitivity_parameters$targets,
               iv_baseline_parameters$targets)) {
  stop("Non-target KEYNOTE parameters changed.")
}


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

iv_sensitivity_parameters_by_scenario <- rho_make_generator_parameters_by_scenario(iv_case, iv_rho_level)
iv_sensitivity_parameters <- iv_sensitivity_parameters_by_scenario$H0
if (iv_method == "RMUT" && nzchar(Sys.getenv("RMUT_WEIGHT", unset = ""))) {
  stop("RMUT_WEIGHT override is forbidden in KM-target rho sensitivity validation.")
}
iv_write_km_metadata()

run_validation_scenario <- function(label, scenario) {
  expected_cache_context <- rho_cache_context(
    c(entry = iv_script_file, dgm = iv_data_generator_file,
      helper = iv_rho_helper, parameters = iv_parameter_table_path,
      analysis_helpers = project_source_path("analysis_helpers.R"),
      source_paths = project_source_path("source_path.R")),
    case = iv_case, method = iv_method, scenario = label, n = N.sim_realcase,
    data_seed = iv_seed_from_label(label), model_seed_base = iv_model_seed_base,
    rho_level = iv_rho_level)
  scenario_tag <- if (identical(scenario, "H1_calibration")) "H1_standard" else scenario
  cache_file <- sprintf("validate_stats_%s_%s.rds", scenario_tag, cache_tag)
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
      cens_upper = trial_settings$cens_upper,
      scenario = scenario,
      parameters = iv_sensitivity_parameters_by_scenario[[scenario]]
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
    
    saveRDS(rho_cache_stamp(stats_list, expected_cache_context), cache_file)
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

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 [KEYNOTE598]", "H0")
run_validation_scenario("H1_standard [KM-target standard effect]", "H1_calibration")
run_validation_scenario("H1_observed [KM-target observed effect]", "H1_observed")

sink()

iv_write_completion_marker()
