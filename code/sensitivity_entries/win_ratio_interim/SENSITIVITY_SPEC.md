# 锁定设计下的中期时点敏感性

评估锁定的 Win Ratio 决策参数在 **50%/80%** 中期时点上的操作特性，并与主分析 **40%/60%/80%** 配对比较。`c_F` 和 `gamma_E` 使用正式参数。

## 分析设定

- 主设计时点：40%、60%、80%、最终。
- 敏感性时点：50%、80%、最终。
- 决策参数、DGM、先验、MCMC、入组、删失和最终分析时间采用正式设定。
- 无效边界：`b(n) = c_F * (n / N_max)^gamma_E`。
- 最终 go：`PP_final > c_F`。
- 中期分析行规则：`interim n-1`；最终：全部 `N.max`。
- 场景：`H0_null`、`H1_calibration`、`H1_observed`。
- `H1_observed` 用于锁参后的描述性评价，不参与选参。
- 中期时点与依赖参数敏感性分别评价。

## 种子与 look_code

- 数据种子：`data_seed_base + replicate_id`，与正式独立评价配对。
- 模型种子：`model_seed_base + replicate_id * 10 + look_code`。
- `look_code` 按分析人数定义：50%→5，80%→3，最终→9。
- 正式40%/60%分析的 `look_code` 为1和2。

## PP路径

- 从 `code/vendor/bayesian_win_ratio/output/validation/<case_key>/formal_10000/cache/` 读取正式80%和最终PP列。
- 计算50% PP并写入本目录 `outputs/`。
- 50% PP按每140个replicate保存checkpoint。

## 运行位置

从仓库根运行 `Rscript code/sensitivity_entries/win_ratio_interim/run_interim_time_sensitivity.R --case=IMpower131`，其他案例键为 `KEYNOTE598`、`OAK` 和 `survival_colon`。入口从自身位置查找敏感性输入，并从 `code/vendor/bayesian_win_ratio/engine/` 加载模型，使用 `cases/` 配置与 `validation/` 锁参表。可用 `PROJECT_ROOT` 指定仓库根。

正式计算需要先完成对应案例的独立验证。`--smoke` 进行烟雾计算；`--dry-run` 检查路径；`--replay-only` 读取已有PP；`--rebuild-look50` 重算50% PP。批量运行可调用本目录的 `run_all_interim_time_sensitivity.sh`；`CONTROL_ROOT` 默认为该脚本目录，`RUN_R` 指定 Rscript。
