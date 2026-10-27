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

# Calibration uses H0 and H1-standard posterior-probability paths.
# H1-observed outputs are descriptive and do not enter parameter selection.
# Formal independent evaluation uses the separate indexed validation entry.

# Method: Cox; case: IMpower131
# Reference decision parameters are used only for diagnostic summaries:

# Formal run: 10,000 newly generated calibration replicates.
# Futility rule: stop at interim look k when PP <= b_k.
# This file generates calibration PP caches and descriptive simulation summaries.
# Runtime outputs use the configured output directory.
iv_version <- "CAL10000_V2_COMMON_CUTOFF_41_20260801"
iv_project_id <- "cox_IMpower"
iv_method <- "Cox"
iv_case <- "IMpower131"
iv_locked_pp <- 0.90
iv_locked_lambda <- 0.45
iv_data_seed_base <- 641001L
iv_model_seed_base <- 741001L

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this calibration-cache file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_project_dir <- normalizePath(file.path(iv_version_dir, ".."), mustWork = TRUE)
setwd(iv_project_dir)

iv_smoke_test <- identical(Sys.getenv("CAL_SMOKE_TEST", unset = "0"), "1")
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
    paste0("run_tag=", iv_run_tag),
    "purpose=parameter_calibration_cache_generation",
    "decision_parameters_selected=false",
    "reference_metrics_only=true",
    paste0("reference_PP_final_threshold=", format(iv_locked_pp, digits = 15)),
    paste0("reference_lambda_prime_bayes=", format(iv_locked_lambda, digits = 15)),
    "futility_rule=PP <= b_k",
    "censoring=random censoring disabled with Uniform(0,1e6); common calendar cutoff applied",
    "final_analysis_time=41 months from first enrollment",
    "final_followup=pmin(generated_censoring, 41-entry_time)",
    paste0("data_seed_base=", iv_data_seed_base),
    paste0("model_seed_base=", iv_model_seed_base),
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  file.path(iv_output_dir, "RUN_METADATA.txt")
)




# 1. 加载依赖并配置编译环境
# Optional compiler settings may be supplied through R_MAKEVARS_USER.
Sys.setenv(PKG_CXXFLAGS = "-DTBB_INTERFACE_NEW")
old_path <- Sys.getenv("PATH")
r_env_bin <- Sys.getenv("R_ENV_BIN", unset = "")
if (nzchar(r_env_bin) && dir.exists(r_env_bin)) {
  Sys.setenv(PATH = paste(r_env_bin, old_path, sep = .Platform$path.sep))
}

options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

library(survival)
library(brms)
library(posterior)

# 2. 分析配置
sink(file.path(iv_output_dir, "calibration_cache_build_results.txt"), split = TRUE)


# 4. 决策参数
PP_final_threshold <- iv_locked_pp
lambda_prime_bayes <- iv_locked_lambda
N.max <- 680
interim_fracs <- c(0.40, 0.60, 0.80)
recruit.int <- 22 / 680 # 入组间隔（22月入组680人）
bayes_iter <- 2000
bayes_warmup <- 1000 # 每链保留 1000；4 链共 4000 个 post-warmup draws
cache_tag <- "ctrl_6p1202_13p4685_H1_8_18_real_7p5610_14p3431_cox_commonCutoff41_normal_0_1p5_iter2000_warmup1000_v7"

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


# Simulate the configured scenarios.
cat("\n=========================================================================\n")
cat("正在对三场景进行自适应设计模拟验证 (Cox 法)...\n")
cat("=========================================================================\n")

library(parallel)
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = 1)
}
source(project_source_path("dgm", "IMpower_shared.R"), local = TRUE)

cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30

n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal calibration must use exactly 10,000 replicates.")
}
# 删失感知参数。第一维必须是复发时间，不是 PFS。
mu.ctrl <- c(log(6.1202), log(13.4685))   # 对照组：复发时间 / 死亡时间
Sigma.ctrl <- matrix(c(0.4209, 0.3287, 0.3287, 1.0268), nrow = 2)

mu.trt.real <- c(log(7.5610), log(14.3431))   # 试验组：复发时间 / 死亡时间
Sigma.trt.real <- matrix(c(0.8207, 0.5298, 0.5298, 1.3679), nrow = 2)

# 标准疗效 H1 场景（复发时间=8.00月，死亡时间=18.00月）
mu.trt.H1_standard <- c(log(8.00), log(18.00))
Sigma.trt.H1_standard <- matrix(c(0.4209, 0.3287, 0.3287, 1.0268), nrow = 2)

run_validation_scenario <- function(label, mu_trt, Sigma_trt, mu_ctrl, Sigma_ctrl) {
  scenario_tag <- ifelse(grepl("H0", label), "H0", ifelse(grepl("H1_standard", label), "H1_standard", "Realcase"))
  cache_file <- sprintf("calibration_stats_%s_%s.rds", scenario_tag, cache_tag)
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
    }, mc.cores = num_cores_validate)
    
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

# H0/H1-standard support calibration; H1-observed is descriptive only.
run_validation_scenario("H0 (原假设/药物无效) [校准用场景]", mu.ctrl, Sigma.ctrl, mu.ctrl, Sigma.ctrl)
run_validation_scenario("H1_standard (复发时间=8.0, 死亡时间=18.0) [标准疗效场景]", mu.trt.H1_standard, Sigma.trt.H1_standard, mu.ctrl, Sigma.ctrl)
run_validation_scenario("真实 IMpower131 生存参数 [模拟验证用场景]", mu.trt.real, Sigma.trt.real, mu.ctrl, Sigma.ctrl)

sink()
