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
# Method: RMST; case: KEYNOTE-598
# Locked decision parameters: PP_final_threshold = 0.864;
#                             lambda_prime_bayes = 0.210
# Formal run: 10,000 newly generated replicates with independent RNG streams.
# Futility rule: stop at interim look k when PP <= b_k.
# This file must not be used for calibration or parameter selection.
# Runtime outputs use the configured output directory.
iv_version <- "IV10000_STANDARD"
iv_project_id <- "RMST_KEYNOTE598"
iv_method <- "RMST"
iv_case <- "KEYNOTE-598"
iv_candidate_id <- "4076"
iv_locked_pp <- 0.864
iv_locked_lambda <- 0.210
iv_data_seed_base <- 1822001L
iv_model_seed_base <- 1922001L

iv_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(iv_file_arg) != 1L) {
  stop("Run this independent-validation file with Rscript.")
}
iv_script_file <- normalizePath(sub("^--file=", "", iv_file_arg), mustWork = TRUE)
iv_version_dir <- dirname(iv_script_file)
iv_project_dir <- normalizePath(file.path(iv_version_dir, ".."), mustWork = TRUE)
setwd(iv_project_dir)

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

source(
  project_source_path("rmst", "validation_helpers.R"),
  local = TRUE
)

iv_seed_from_label <- monitoring_seed_from_label

writeLines(
  c(
    paste0("version=", iv_version),
    paste0("project=", iv_project_id),
    paste0("method=", iv_method),
    paste0("case=", iv_case),
    paste0("candidate_id=", iv_candidate_id),
    paste0("run_tag=", iv_run_tag),
    paste0("locked_PP_final_threshold=", format(iv_locked_pp, digits = 15)),
    paste0("locked_lambda_prime_bayes=", format(iv_locked_lambda, digits = 15)),
    "futility_rule=PP <= b_k",
    paste0("data_seed_base=", iv_data_seed_base),
    paste0("model_seed_base=", iv_model_seed_base),
    paste0("formal_data_seed_base=", iv_formal_data_seed_base),
    paste0("formal_model_seed_base=", iv_formal_model_seed_base),
    "seed_policy=smoke test uses +500000000 offset; formal streams remain untouched",
    "replicate_output=replicate_id,fit_ok,failure_reason,PP_interims,PP_final,early_futility,stop_stage,sample_size_used,go,final_no_go",
    "A05_output=probability,MCSE,Wilson_95pct_Monte_Carlo_interval",
    "A06_output=direct_final_no_go_and_partition_identity_check",
    "failure_handling=reported_both_valid_only_and_all_replicates_failure_as_no_go",
    paste0("started_at=", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
  ),
  file.path(iv_output_dir, "RUN_METADATA.txt")
)


# Locked candidate 4076: PP = 0.864, lambda = 0.210.

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

# Simulation configuration and result logging.
n_sim_realcase_env <- Sys.getenv("N_SIM_REALCASE")
N.sim_realcase <- if (n_sim_realcase_env != "") as.integer(n_sim_realcase_env) else 10000
if (!iv_smoke_test && N.sim_realcase != 10000L) {
  stop("Formal independent validation must use exactly 10,000 replicates.")
}
results_file <- if (N.sim_realcase == 10000) {
  "independent_validation_results_10k_KEYNOTE598_v2_tau4_6_9_12_priorInterceptLog20_weak.txt"
} else {
  "independent_validation_results_KEYNOTE598_v2_tau4_6_9_12_priorInterceptLog20_weak.txt"
}
results_file <- file.path(iv_output_dir, basename(results_file))
sink(results_file, split = TRUE)

source(project_source_path("dgm", "KEYNOTE598_shared.R"), local = TRUE)
trial_settings <- keynote598_trial_settings(tau_fixed = c(4, 6, 9, 12))

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
tau_fixed <- trial_settings$tau_fixed
bayes_iter <- 8000L
bayes_warmup <- 4000L
bayes_adapt <- 0.95
bayes_treedepth <- 12L

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
    prior(normal(log(20), 1.5), class = "Intercept"),
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
source(project_source_path("dgm", "KEYNOTE598_shared.R"), local = TRUE)


cores_env <- Sys.getenv("NUM_CORES")
num_cores_validate <- if (cores_env != "") as.integer(cores_env) else 30



run_validation_scenario <- rmst_KEYNOTE598_validation_scenario

# 依次执行三个场景的模拟验证 (校准与验证场景分离)
run_validation_scenario("H0 [KEYNOTE598]", "H0")
run_validation_scenario("H1_calibration [standard effect]", "H1_calibration")
run_validation_scenario("H1_observed [pembrolizumab + ipilimumab observed effect]", "H1_observed")

sink()
cat(sprintf(
  "\nH0、H1-standard 和 H1-observed 三种模拟场景均已完成，结果保存在 %s\n",
  results_file
))
