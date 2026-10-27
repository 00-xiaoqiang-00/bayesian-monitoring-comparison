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

# Method: Cox; case: OAK
# Purpose: generate fresh H0 and H1-standard PP caches for parameter calibration.
# H1-standard treatment component-level log-scale dispersions equal control
# in the control arm.
# Futility rule used for diagnostic replay: PP <= b_k.
# Generates H0/H1-standard calibration probabilities and grid summaries.
# Final parameter selection is performed by the documented selection entry.

# 1. 加载依赖并配置编译环境
# Optional compiler settings may be supplied through R_MAKEVARS_USER.
Sys.setenv(PKG_CXXFLAGS = "-DTBB_INTERFACE_NEW")
old_path <- Sys.getenv("PATH")
r_env_bin <- Sys.getenv("R_ENV_BIN", unset = "")
if (nzchar(r_env_bin) && dir.exists(r_env_bin)) {
  Sys.setenv(PATH = paste(r_env_bin, old_path, sep = .Platform$path.sep))
}

options(repos = c(CRAN = "https://mirrors.ustc.edu.cn/CRAN/"))

library(MASS)
library(survival)
library(brms)
library(parallel)
library(posterior)

cal_version <- "CAL10000_V2_OAKV5_EQUAL_H1_VARIANCE_20260802"
cal_data_seed_h0 <- 1813001L
cal_data_seed_h1 <- 1813002L
cal_model_seed_base <- 1913001L
cal_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(cal_file_arg) != 1L) stop("Run this file with Rscript.")
cal_script_file <- normalizePath(sub("^--file=", "", cal_file_arg), mustWork = TRUE)
cal_version_dir <- dirname(cal_script_file)
cal_project_dir <- normalizePath(file.path(cal_version_dir, ".."), mustWork = TRUE)
setwd(cal_project_dir)

cal_smoke_test <- identical(Sys.getenv("CAL_SMOKE_TEST", unset = "0"), "1")
if (cal_smoke_test) {
  cal_smoke_n <- suppressWarnings(as.integer(Sys.getenv("N_SIM", unset = "")))
  if (!is.finite(cal_smoke_n) || cal_smoke_n < 1L || cal_smoke_n > 20L) {
    stop("Smoke mode requires N_SIM between 1 and 20.")
  }
  cal_run_tag <- sprintf("smoke_n%d", cal_smoke_n)
} else {
  Sys.setenv(N_SIM = "10000", N_SIM_H0 = "10000", N_SIM_H1 = "10000")
  cal_run_tag <- "formal_10000"
}

cal_output_dir <- file.path(cal_version_dir, "outputs", cal_run_tag)
cal_cache_dir <- file.path(cal_output_dir, "cache")
dir.create(cal_cache_dir, recursive = TRUE, showWarnings = FALSE)

recruit.int <- 8.6 / 850

# 载入数据生成函数
source(project_source_path("dgm", "OAK_shared.R"), local = TRUE)

# --- 校准与模拟参数设置 ---
n_sim_env    <- Sys.getenv("N_SIM")
n_sim_h0_env <- Sys.getenv("N_SIM_H0")
n_sim_h1_env <- Sys.getenv("N_SIM_H1")

N.sim_h0 <- if (n_sim_h0_env != "") as.integer(n_sim_h0_env) else (if (n_sim_env != "") as.integer(n_sim_env) else 10000)
N.sim_h1 <- if (n_sim_h1_env != "") as.integer(n_sim_h1_env) else (if (n_sim_env != "") as.integer(n_sim_env) else 10000)
if (!cal_smoke_test && (N.sim_h0 != 10000L || N.sim_h1 != 10000L)) {
  stop("Formal calibration requires exactly 10,000 H0 and 10,000 H1-standard replicates.")
}

cores_env <- Sys.getenv("NUM_CORES")
num_cores <- if (cores_env != "") as.integer(cores_env) else 30     # 并行核数

N.max        <- 850
cens_upper   <- 109.3
interim_fracs<- c(0.40, 0.60, 0.80)
alpha_final  <- 0.05      # 频率派初始显著性水平
bayes_iter   <- 2000      # 每链总迭代数
bayes_warmup <- 1000      # 每链保留 1000；4 链共 4000 个 post-warmup draws
cache_tag <- "OAK_mix_v5_equal_H1_variance_cox_normal_0_1p5_iter2000_warmup1000"
cache_h0_path <- file.path(cal_cache_dir, sprintf("calib_stats_H0_%s.rds", cache_tag))
cache_h1_path <- file.path(cal_cache_dir, sprintf("calib_stats_H1_%s.rds", cache_tag))

writeLines(
  c(
    paste0("version=", cal_version),
    "project=cox_OAK",
    "method=Cox",
    "case=OAK",
    paste0("run_tag=", cal_run_tag),
    "purpose=parameter_calibration_cache_generation",
    "decision_parameters_selected=false",
    "dgm_version=OAK_v5_equal_H1_component_variance",
    "dgm_source=code/common/dgm/OAK_shared.R",
    "dgm_source_sha256=5d436743409a08e875e4eb05019dea97af8a5841ea0fd7aa51847e9d64f82c83",
    "scenarios=H0,H1-standard",
    "futility_rule=PP <= b_k",
    "N_max=850",
    "interim_fractions=0.40,0.60,0.80",
    "final_followup=27.9 months from first enrollment",
    "mcmc_chains=4",
    paste0("mcmc_iter=", bayes_iter),
    paste0("mcmc_warmup=", bayes_warmup),
    paste0("data_seed_H0=", cal_data_seed_h0),
    paste0("data_seed_H1_standard=", cal_data_seed_h1),
    paste0("model_seed_base=", cal_model_seed_base),
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  file.path(cal_output_dir, "RUN_METADATA.txt")
)

# H0 and H1_calibration arm parameters are selected inside the OAK generator.

# --- MCMC 诊断与后验概率 ---
extract_valid_pp <- monitoring_extract_valid_pp

# --- 模拟与提取统计量函数 ---
run_simulation_and_extract_stats <- function(sim_array, base_model, N.sim, label) {
  cat(sprintf("\n【MCMC】开始对 [%s] 模拟数据进行采样并提取指标 (总次数: %d, 核数: %d)...\n", label, N.sim, num_cores))
  entry_times <- recruit.int * (0:(N.max - 1))
  
  stats_list <- mclapply(seq_len(N.sim), function(j) {
    df_j <- as.data.frame(sim_array[, , j])
    
    # 复合终点
    time_to_first <- pmin(df_j$death_t, df_j$recurrence_t)
    final_censor_t <- pmin(df_j$censor_t, oak_admin_followup(entry_times, 27.9))
    df_j$time   <- pmin(time_to_first, final_censor_t)
    df_j$status <- as.integer(time_to_first <= final_censor_t)

    interim_CP <- numeric(length(interim_fracs))
    interim_PP <- numeric(length(interim_fracs))
    simulation_ok <- TRUE
    
    original_censor_t <- df_j$censor_t
    
    for (k in seq_along(interim_fracs)) {
      n_k <- floor(N.max * interim_fracs[k])

      # 行政随访截断
      Time.current <- recruit.int * (n_k - 1)
      obs_window <- pmax(Time.current - entry_times[1:n_k], 0)
      censor_t_interim <- pmin(original_censor_t[1:n_k], obs_window)
      
      time_interim   <- pmin(time_to_first[1:n_k], censor_t_interim)
      status_interim <- as.integer(time_to_first[1:n_k] <= censor_t_interim)
      
      valid_rows <- censor_t_interim > 0
      df_interim <- data.frame(
        time   = time_interim[valid_rows],
        status = status_interim[valid_rows],
        group  = df_j$group[1:n_k][valid_rows]
      )
      
      # 1. 频率派 CP 计算
      fit_interim <- try(coxph(Surv(time, status) ~ group, data = df_interim), silent = TRUE)
      if (!inherits(fit_interim, "try-error") && fit_interim$nevent > 0) {
        s <- summary(fit_interim)
        beta_hat <- s$coefficients[1, "coef"]
        se       <- s$coefficients[1, "se(coef)"]
        
        t_k <- n_k / N.max
        Z_k <- -beta_hat / se
        CP_k <- 1 - pnorm((qnorm(1 - alpha_final) - sqrt(t_k) * Z_k) / sqrt(1 - t_k))
        interim_CP[k] <- CP_k
      } else {
        interim_CP[k] <- 0
      }
      
      # 2. 计算 MCMC 后验概率 PP_k。
      fit_bayes_interim <- try(
        update(
          base_model,
          newdata = df_interim,
          chains  = 4,
          cores   = 1,
          iter    = bayes_iter,
          warmup  = bayes_warmup,
          seed    = cal_model_seed_base + j + k * 10000L,
          refresh = 0,
          silent  = 2
        ),
        silent = TRUE
      )
      
      interim_PP[k] <- extract_valid_pp(fit_bayes_interim)
      if (is.na(interim_PP[k])) {
        simulation_ok <- FALSE
        break
      }
    }

    if (!simulation_ok) {
      return(list(ok = FALSE, CP = interim_CP, PP = interim_PP,
                  pval = NA_real_, PP_final = NA_real_))
    }
    
    # 最终分析 - 频率派
    fit_final <- try(coxph(Surv(time, status) ~ group, data = df_j), silent = TRUE)
    if (!inherits(fit_final, "try-error")) {
      s <- summary(fit_final)
      pval_final <- s$coefficients[1, "Pr(>|z|)"]
    } else {
      pval_final <- 1.0
    }
    
    # 最终分析的贝叶斯 MCMC。
    df_final_clean <- df_j[, c("time", "status", "group")]
    fit_bayes_final <- try(
      update(
        base_model,
        newdata = df_final_clean,
        chains  = 4,
        cores   = 1,
        iter    = bayes_iter,
        warmup  = bayes_warmup,
        seed    = cal_model_seed_base + j + 90000L,
        refresh = 0,
        silent  = 2
      ),
      silent = TRUE
    )
    PP_final <- extract_valid_pp(fit_bayes_final)
    if (is.na(PP_final)) {
      return(list(ok = FALSE, CP = interim_CP, PP = interim_PP,
                  pval = pval_final, PP_final = NA_real_))
    }
    
    list(ok = TRUE, CP = interim_CP, PP = interim_PP,
         pval = pval_final, PP_final = PP_final)
  }, mc.cores = num_cores)
  
  return(stats_list)
}

# --- 缓存加载或重新计算 ---
use_cache <- FALSE
force_recalculate <- Sys.getenv("FORCE_RECALC") == "TRUE"

if (!force_recalculate && file.exists(cache_h0_path) && file.exists(cache_h1_path)) {
  stats_list_H0_temp <- readRDS(cache_h0_path)
  stats_list_H1_temp <- readRDS(cache_h1_path)
  if (length(stats_list_H0_temp) == N.sim_h0 && length(stats_list_H1_temp) == N.sim_h1) {
    cat("\nCompatible MCMC cache found for the current N_SIM; loading it...\n")
    stats_list_H0 <- stats_list_H0_temp
    stats_list_H1 <- stats_list_H1_temp
    use_cache <- TRUE
  }
}

if (!use_cache) {
  # 生成 H0 模拟数据
  cat(sprintf("\n【1】在 H0（原假设）下生成 %d 次模拟数据集...\n", N.sim_h0))
  set.seed(cal_data_seed_h0)
  sim_array_H0 <- data.simulation(
    N.sim = N.sim_h0, N.max = N.max, 
    cens_upper = cens_upper,
    scenario = "H0"
  )

  # 生成 H1 模拟数据
  cat(sprintf("\n【2】在 H1（标准疗效）下生成 %d 次模拟数据集...\n", N.sim_h1))
  set.seed(cal_data_seed_h1)
  sim_array_H1 <- data.simulation(
    N.sim = N.sim_h1, N.max = N.max, 
    cens_upper = cens_upper,
    scenario = "H1_calibration"
  )

  # 预编译基础贝叶斯模型
  cat("\n正在为所有分析预编译通用的基础贝叶斯 Cox 模型...\n")
  df_dummy <- data.frame(
    time   = runif(N.max, 0.5, cens_upper),
    status = sample(c(0, 1), N.max, replace = TRUE),
    group  = sample(c(0, 1), N.max, replace = TRUE)
  )
  base_model <- brm(
    formula = time | cens(1 - status) ~ group,
    data    = df_dummy,
    family  = brmsfamily("cox"),
    prior   = prior(normal(0, 1.5), class = "b", coef = "group"),
    chains  = 4,
    cores   = 1,
    iter    = bayes_iter,
    warmup  = bayes_warmup,
    seed    = 12345,
    refresh = 0,
    silent  = 2
  )
  cat("Base model precompilation completed.即将开始并行采样分析。\n\n")

  # 运行 MCMC 并提取结果
  stats_list_H0 <- run_simulation_and_extract_stats(sim_array_H0, base_model, N.sim_h0, "H0")
  saveRDS(stats_list_H0, cache_h0_path)

  stats_list_H1 <- run_simulation_and_extract_stats(sim_array_H1, base_model, N.sim_h1, "H1")
  saveRDS(stats_list_H1, cache_h1_path)
}

# 过滤并提取结果
successful_runs_H0 <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list_H0)
N.valid_H0 <- length(successful_runs_H0)

successful_runs_H1 <- Filter(function(x) is.list(x) && isTRUE(x$ok), stats_list_H1)
N.valid_H1 <- length(successful_runs_H1)

if (N.valid_H0 < N.sim_h0) {
  cat(sprintf("警告：H0 仿真有 %d 次报错，实际成功 %d 次。\n", N.sim_h0 - N.valid_H0, N.valid_H0))
}
if (N.valid_H1 < N.sim_h1) {
  cat(sprintf("警告：H1 仿真有 %d 次报错，实际成功 %d 次。\n", N.sim_h1 - N.valid_H1, N.valid_H1))
}

# --- 辅助评价函数 ---
get_h0_metrics <- function(lam, pp_thresh, stats_list, N.valid) {
  c_primes <- pp_thresh * (floor(N.max * interim_fracs) / N.max)^lam
  h0_stopped_count <- 0
  h0_type_I_count  <- 0
  ess_sum <- 0
  for (j in seq_len(N.valid)) {
    stats <- stats_list[[j]]
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      h0_stopped_count <- h0_stopped_count + 1
      ess_sum <- ess_sum + N.max * interim_fracs[stop_idx[1]]
    } else {
      ess_sum <- ess_sum + N.max
      if (!is.na(stats$PP_final) && stats$PP_final > pp_thresh) {
        h0_type_I_count <- h0_type_I_count + 1
      }
    }
  }
  return(list(
    t1e = h0_type_I_count / N.valid, 
    stop_rate = h0_stopped_count / N.valid,
    ess = ess_sum / N.valid
  ))
}

evaluate_h1 <- function(lam, pp_thresh, stats_list_H1, N.valid_H1) {
  c_primes <- pp_thresh * (floor(N.max * interim_fracs) / N.max)^lam
  h1_stopped_count <- 0
  h1_success_count <- 0
  ess_sum <- 0
  for (j in seq_len(N.valid_H1)) {
    stats <- stats_list_H1[[j]]
    stop_idx <- which(stats$PP <= c_primes)
    if (length(stop_idx) > 0) {
      h1_stopped_count <- h1_stopped_count + 1
      ess_sum <- ess_sum + N.max * interim_fracs[stop_idx[1]]
    } else {
      ess_sum <- ess_sum + N.max
      if (!is.na(stats$PP_final) && stats$PP_final > pp_thresh) {
        h1_success_count <- h1_success_count + 1
      }
    }
  }
  return(list(
    Power = h1_success_count / N.valid_H1,
    H1_stop_rate = h1_stopped_count / N.valid_H1,
    ess = ess_sum / N.valid_H1
  ))
}

# --- 网格扫描评估与输出 ---
cat("\n【3】开始执行决策参数网格扫描并保存校准结果...\n")

PP_thresh_grid <- seq(0.90, 0.95, by = 0.01)
lambda_grid <- seq(0.20, 2.00, by = 0.05)

candidates_df <- data.frame(
  PP_Thresh = numeric(),
  Lambda = numeric(),
  H0_T1E = numeric(),
  H0_StopRate = numeric(),
  H0_ESS = numeric(),
  H1_Power = numeric(),
  H1_StopRate = numeric(),
  H1_ESS = numeric(),
  stringsAsFactors = FALSE
)

for (pp_thresh in PP_thresh_grid) {
  for (lam in lambda_grid) {
    h0_res <- get_h0_metrics(lam, pp_thresh, successful_runs_H0, N.valid_H0)
    h1_res <- evaluate_h1(lam, pp_thresh, successful_runs_H1, N.valid_H1)
    
    candidates_df <- rbind(candidates_df, data.frame(
      PP_Thresh = pp_thresh,
      Lambda = lam,
      H0_T1E = h0_res$t1e,
      H0_StopRate = h0_res$stop_rate,
      H0_ESS = h0_res$ess,
      H1_Power = h1_res$Power,
      H1_StopRate = h1_res$H1_stop_rate,
      H1_ESS = h1_res$ess
    ))
  }
}

# Candidate-grid summaries do not select or lock decision parameters.

# 写入 txt 结果文件
write.csv(
  candidates_df,
  file.path(cal_output_dir, "calibration_candidate_grid_reference.csv"),
  row.names = FALSE
)
output_path <- file.path(cal_output_dir, "calibration_candidate_grid_reference.txt")

# 构造表格文本
write_output <- monitoring_write_calibration_output

# 同时输出到控制台和文件
sink(output_path)
write_output()
sink()

# 输出候选网格到终端
write_output()
cat(sprintf("\n完整的校准评估列表已保存至：%s\n", output_path))
