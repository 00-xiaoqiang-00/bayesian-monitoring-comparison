from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[2]
DEFAULT_ARCHIVE = PACKAGE / "data/estimands/H1_standard_true_estimands.csv"


@dataclass(frozen=True)
class Component:
    weight: float
    mu_r: float
    sigma_r: float
    mu_d: float
    sigma_d: float
    rho: float


@dataclass(frozen=True)
class Arm:
    components: tuple[Component, ...]


@dataclass(frozen=True)
class Trial:
    nmax: int
    recruitment_duration: float
    final_analysis_time: float
    cens_upper: float
    tau_rmst: float
    tau_rmut: float
    w_r: float = 0.5

    @property
    def recruit_interval(self) -> float:
        return self.recruitment_duration / self.nmax


@dataclass(frozen=True)
class CaseSpec:
    case: str
    version: str
    control: Arm
    standard: Arm
    observed: Arm
    trial: Trial
    generator: str
    status: str


def bivar_arm(mu_r: float, mu_d: float, var_r: float, cov_rd: float, var_d: float) -> Arm:
    sigma_r = math.sqrt(var_r)
    sigma_d = math.sqrt(var_d)
    rho = cov_rd / (sigma_r * sigma_d)
    return Arm((Component(1.0, mu_r, sigma_r, mu_d, sigma_d, rho),))


def mixture_arm(rows: list[tuple[float, float, float, float, float, float]]) -> Arm:
    return Arm(tuple(Component(*row) for row in rows))


def build_specs() -> list[CaseSpec]:
    colon_ctrl = bivar_arm(
        math.log(38.1268), math.log(64.5907),
        3.29767175, 2.19701400, 1.77436755,
    )
    colon_standard = bivar_arm(
        math.log(92.6930), math.log(116.2198),
        3.29767175, 2.19701400, 1.77436755,
    )
    colon_observed = bivar_arm(
        math.log(87.1562), math.log(111.4653),
        4.46376802, 3.41770993, 2.92813972,
    )

    oak_ctrl = mixture_arm([
        (0.13123216948698738, 0.2666273699841142, 0.3000000000000000,
         2.3088355183235585, 1.0177, 0.7),
        (1 - 0.13123216948698738, 1.6274656470707884, 0.8489991311979127,
         2.3088355183235585, 1.0177, 0.7),
    ])
    # OAK H1-standard changes the treatment locations and mixture weight.
    # Component log-scale dispersions and correlations match control.
    oak_standard = mixture_arm([
        (0.10183833564436579, 0.510611455725049, 0.3000000000000000,
         2.636754680552583, 1.0177, 0.7),
        (1 - 0.10183833564436579, 1.9615174106704887, 0.8489991311979127,
         2.636754680552583, 1.0177, 0.7),
    ])
    oak_observed = mixture_arm([
        (0.2623242204350347, 0.3645597370078664, 0.0806101279283579,
         2.636754680552583, 1.1709, 0.7),
        (1 - 0.2623242204350347, 1.677883933719392, 1.2693769377027608,
         2.636754680552583, 1.1709, 0.7),
    ])

    keynote_ctrl = mixture_arm([
        (0.2070545742, math.log(2.4837219024), 0.4340057961,
         math.log(2.7290953670), 1.5000000000, 0.5),
        (0.7929454258, math.log(16.4600817241), 1.5317036442,
         math.log(28.6890202294), 0.9670082477, 0.5),
    ])
    keynote_standard = mixture_arm([
        (0.2070545742, math.log(2.4837219024) + 0.5946115916, 0.4340057961,
         math.log(2.7290953670) + 0.3773863423, 1.5000000000, 0.5),
        (0.7929454258, math.log(16.4600817241) + 0.5946115916, 1.5317036442,
         math.log(28.6890202294) + 0.3773863423, 0.9670082477, 0.5),
    ])
    keynote_observed = mixture_arm([
        (0.1936579096, math.log(2.0212342835), 0.1500000000,
         math.log(2.5120209874), 1.5000000000, 0.5),
        (0.8063420904, math.log(15.2003022720), 1.3810016053,
         math.log(29.2329569146), 1.2402934461, 0.5),
    ])

    imp_formal_ctrl = bivar_arm(
        math.log(6.1202), math.log(13.4685), 0.4209, 0.3287, 1.0268
    )
    imp_formal_standard = bivar_arm(
        math.log(8.0), math.log(18.0), 0.4209, 0.3287, 1.0268
    )
    imp_formal_observed = bivar_arm(
        math.log(7.5610), math.log(14.3431), 0.8207, 0.5298, 1.3679
    )

    return [
        CaseSpec(
            "survival_colon", "formal_v1", colon_ctrl, colon_standard, colon_observed,
            Trial(620, 43, 109, 1900, 60, 60),
            "survival_colon_bivariate_lognormal_v1", "CONSISTENT_ACROSS_METHODS",
        ),
        CaseSpec(
            "OAK", "formal_v5_equal_H1_variance", oak_ctrl, oak_standard, oak_observed,
            Trial(850, 8.6, 27.9, 109.3, 18, 18),
            "OAK_shared_death_progression_mixture_v5_equal_H1_variance",
            "CONSISTENT_ACROSS_METHODS_EQUAL_COMPONENT_VARIANCE",
        ),
        CaseSpec(
            "KEYNOTE-598", "formal_v1", keynote_ctrl, keynote_standard, keynote_observed,
            Trial(568, 19.3, 31.7, 141.8559446075, 12, 18),
            "keynote598_shared_component_lognormal_v1",
            "DGM_CONSISTENT_TAU_DIFFERS_RMST_VS_RMUT",
        ),
        CaseSpec(
            "IMpower131", "formal_E05_v6_common_cutoff_41", imp_formal_ctrl,
            imp_formal_standard, imp_formal_observed,
            Trial(680, 22, 41, 1e6, 18, 18),
            "IMpower131_bivariate_lognormal_v6_E05_locked_20260801",
            "CONSISTENT_ACROSS_METHODS_LOCKED",
        ),
    ]


def arm_survival(arm: Arm, time: np.ndarray, endpoint: str) -> np.ndarray:
    time = np.asarray(time, dtype=float)
    result = np.zeros_like(time)
    positive = time > 0
    result[~positive] = 1.0
    for component in arm.components:
        if endpoint == "D":
            z_d = (np.log(time[positive]) - component.mu_d) / component.sigma_d
            survival = norm.sf(z_d)
        elif endpoint == "C":
            z_r = (np.log(time[positive]) - component.mu_r) / component.sigma_r
            z_d = (np.log(time[positive]) - component.mu_d) / component.sigma_d
            joint_cdf = multivariate_normal.cdf(
                np.column_stack([z_r, z_d]),
                mean=np.zeros(2),
                cov=np.array([[1.0, component.rho], [component.rho, 1.0]]),
            )
            survival = 1.0 - norm.cdf(z_r) - norm.cdf(z_d) + joint_cdf
        else:
            raise ValueError(endpoint)
        result[positive] += component.weight * survival
    return np.clip(result, 0.0, 1.0)


def integrate_curve_difference(control: Arm, treatment: Arm, tau: float, endpoint: str) -> float:
    grid = np.linspace(0.0, tau, 20001)
    difference = arm_survival(treatment, grid, endpoint) - arm_survival(control, grid, endpoint)
    return float(np.trapezoid(difference, grid))


def draw_arm(arm: Arm, n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    weights = np.array([c.weight for c in arm.components], dtype=float)
    component = rng.choice(len(arm.components), size=n, p=weights / weights.sum())
    z_r = rng.standard_normal(n)
    z_independent = rng.standard_normal(n)
    mu_r = np.array([c.mu_r for c in arm.components])[component]
    sigma_r = np.array([c.sigma_r for c in arm.components])[component]
    mu_d = np.array([c.mu_d for c in arm.components])[component]
    sigma_d = np.array([c.sigma_d for c in arm.components])[component]
    rho = np.array([c.rho for c in arm.components])[component]
    z_d = rho * z_r + np.sqrt(1.0 - rho * rho) * z_independent
    return np.exp(mu_r + sigma_r * z_r), np.exp(mu_d + sigma_d * z_d)


def draw_censor(trial: Trial, n: int, rng: np.random.Generator) -> np.ndarray:
    entry_index = rng.integers(0, trial.nmax, size=n)
    administrative = trial.final_analysis_time - entry_index * trial.recruit_interval
    random_censor = rng.uniform(0.0, trial.cens_upper, size=n)
    return np.minimum(administrative, random_censor)


def cox_hr_one_batch(spec: CaseSpec, n_arm: int, rng: np.random.Generator) -> float:
    r0, d0 = draw_arm(spec.control, n_arm, rng)
    r1, d1 = draw_arm(spec.standard, n_arm, rng)
    c0 = draw_censor(spec.trial, n_arm, rng)
    c1 = draw_censor(spec.trial, n_arm, rng)
    tc0 = np.minimum(r0, d0)
    tc1 = np.minimum(r1, d1)
    time0 = np.minimum(tc0, c0)
    time1 = np.minimum(tc1, c1)
    event0 = tc0 <= c0
    event1 = tc1 <= c1
    sorted0 = np.sort(time0)
    sorted1 = np.sort(time1)
    event_times0 = time0[event0]
    event_times1 = time1[event1]
    all_event_times = np.concatenate([event_times0, event_times1])
    all_event_groups = np.concatenate([
        np.zeros(event_times0.size, dtype=float),
        np.ones(event_times1.size, dtype=float),
    ])
    risk0 = n_arm - np.searchsorted(sorted0, all_event_times, side="left")
    risk1 = n_arm - np.searchsorted(sorted1, all_event_times, side="left")

    def score(beta: float) -> float:
        eb = math.exp(beta)
        expected = eb * risk1 / (risk0 + eb * risk1)
        return float(np.sum(all_event_groups - expected))

    beta = brentq(score, -5.0, 5.0)
    return math.exp(beta)


def censor_survival(time: np.ndarray, trial: Trial) -> np.ndarray:
    time = np.asarray(time, dtype=float)
    random_survival = np.clip(1.0 - time / trial.cens_upper, 0.0, 1.0)
    max_index = np.floor((trial.final_analysis_time - time) / trial.recruit_interval).astype(int)
    admin_count = np.clip(max_index + 1, 0, trial.nmax)
    admin_survival = admin_count / trial.nmax
    return random_survival * admin_survival


def wr_one_batch(spec: CaseSpec, n_pairs: int, rng: np.random.Generator) -> dict[str, float]:
    rt, dt = draw_arm(spec.standard, n_pairs, rng)
    rc, dc = draw_arm(spec.control, n_pairs, rng)
    ct = draw_censor(spec.trial, n_pairs, rng)
    cc = draw_censor(spec.trial, n_pairs, rng)

    dd_t = dt <= ct
    dd_c = dc <= cc
    dr_t = rt <= np.minimum(dt, ct)
    dr_c = rc <= np.minimum(dc, cc)
    stored_ct = np.minimum(dt, ct)
    stored_cc = np.minimum(dc, cc)

    death_win = (dd_t & dd_c & (dc < dt)) | (~dd_t & dd_c & (dc < stored_ct))
    death_loss = (dd_t & dd_c & (dt < dc)) | (dd_t & ~dd_c & (dt < stored_cc))
    unresolved = ~(death_win | death_loss)
    recurrence_win = unresolved & ((dr_t & dr_c & (rc < rt)) | (~dr_t & dr_c & (rc < stored_ct)))
    recurrence_loss = unresolved & ((dr_t & dr_c & (rt < rc)) | (dr_t & ~dr_c & (rt < stored_cc)))

    win_time = np.where(death_win, dc, np.where(recurrence_win, rc, np.nan))
    loss_time = np.where(death_loss, dt, np.where(recurrence_loss, rt, np.nan))
    win_g = censor_survival(win_time[np.isfinite(win_time)], spec.trial)
    loss_g = censor_survival(loss_time[np.isfinite(loss_time)], spec.trial)
    win_weights = 1.0 / np.square(win_g)
    loss_weights = 1.0 / np.square(loss_g)
    return {
        "weighted_wins": float(np.sum(win_weights)),
        "weighted_losses": float(np.sum(loss_weights)),
        "death_win_n": int(np.sum(death_win)),
        "death_loss_n": int(np.sum(death_loss)),
        "recurrence_win_n": int(np.sum(recurrence_win)),
        "recurrence_loss_n": int(np.sum(recurrence_loss)),
        "unresolved_n": int(n_pairs - np.sum(death_win | death_loss | recurrence_win | recurrence_loss)),
        "pairs": n_pairs,
    }


def compute_estimands(spec: CaseSpec, seed: int) -> dict[str, float | str]:
    delta_rmst = integrate_curve_difference(
        spec.control, spec.standard, spec.trial.tau_rmst, "C"
    )
    delta_sd = integrate_curve_difference(
        spec.control, spec.standard, spec.trial.tau_rmut, "D"
    )
    delta_sc = integrate_curve_difference(
        spec.control, spec.standard, spec.trial.tau_rmut, "C"
    )
    delta_rmut = (1.0 - spec.trial.w_r) * delta_sd + spec.trial.w_r * delta_sc

    rng_hr = np.random.default_rng(seed)
    hr_batches = np.array([
        cox_hr_one_batch(spec, 100_000, rng_hr) for _ in range(12)
    ])
    hr = float(np.exp(np.mean(np.log(hr_batches))))
    hr_mcse = float(np.std(np.log(hr_batches), ddof=1) / math.sqrt(hr_batches.size) * hr)

    rng_wr = np.random.default_rng(seed + 1_000_000)
    wr_batches = [wr_one_batch(spec, 200_000, rng_wr) for _ in range(20)]
    total_wins = sum(x["weighted_wins"] for x in wr_batches)
    total_losses = sum(x["weighted_losses"] for x in wr_batches)
    wr = total_wins / total_losses
    log_wr_batches = np.log([
        x["weighted_wins"] / x["weighted_losses"] for x in wr_batches
    ])
    wr_mcse = float(np.std(log_wr_batches, ddof=1) / math.sqrt(len(log_wr_batches)) * wr)
    total_pairs = sum(x["pairs"] for x in wr_batches)

    return {
        "case": spec.case,
        "version": spec.version,
        "generator": spec.generator,
        "status": spec.status,
        "wr": wr,
        "wr_mcse": wr_mcse,
        "wr_log": math.log(wr),
        "wr_death_win_fraction": sum(x["death_win_n"] for x in wr_batches) / total_pairs,
        "wr_death_loss_fraction": sum(x["death_loss_n"] for x in wr_batches) / total_pairs,
        "wr_recurrence_win_fraction": sum(x["recurrence_win_n"] for x in wr_batches) / total_pairs,
        "wr_recurrence_loss_fraction": sum(x["recurrence_loss_n"] for x in wr_batches) / total_pairs,
        "wr_unresolved_fraction": sum(x["unresolved_n"] for x in wr_batches) / total_pairs,
        "cox_hr": hr,
        "cox_hr_mcse": hr_mcse,
        "delta_rmst_months": delta_rmst,
        "tau_rmst_months": spec.trial.tau_rmst,
        "delta_rmut_months": delta_rmut,
        "tau_rmut_months": spec.trial.tau_rmut,
        "rmut_weight_recurrence_state": spec.trial.w_r,
        "delta_integral_sd_months": delta_sd,
        "delta_integral_sc_months": delta_sc,
        "cox_and_wr_final_trial_process": (
            f"Nmax={spec.trial.nmax}; recruitment={spec.trial.recruitment_duration}; "
            f"final={spec.trial.final_analysis_time}; random_censor=Uniform(0,{spec.trial.cens_upper:g})"
        ),
        "hr_batches": len(hr_batches),
        "hr_n_per_arm_per_batch": 100_000,
        "wr_batches": len(wr_batches),
        "wr_pairs_per_batch": 200_000,
        "seed": seed,
    }


def parameter_rows(specs: list[CaseSpec], scenario: str) -> list[dict[str, object]]:
    rows = []
    for spec in specs:
        arm_map = {"control": spec.control, "treatment": getattr(spec, scenario)}
        for arm_name, arm in arm_map.items():
            for component_index, component in enumerate(arm.components, start=1):
                rows.append({
                    "case": spec.case,
                    "version": spec.version,
                    "status": spec.status,
                    "generator": spec.generator,
                    "scenario": "H1-standard" if scenario == "standard" else "H1-observed",
                    "arm": arm_name,
                    "component": component_index,
                    "weight": component.weight,
                    "mu_log_R": component.mu_r,
                    "median_R_months": math.exp(component.mu_r),
                    "sigma_log_R": component.sigma_r,
                    "mu_log_D": component.mu_d,
                    "median_D_months": math.exp(component.mu_d),
                    "sigma_log_D": component.sigma_d,
                    "rho_log_R_D": component.rho,
                    "cov_log_R_D": component.rho * component.sigma_r * component.sigma_d,
                })
    return rows


def trial_rows(specs: list[CaseSpec]) -> list[dict[str, object]]:
    rows = []
    for spec in specs:
        t = spec.trial
        rows.append({
            "case": spec.case,
            "version": spec.version,
            "status": spec.status,
            "generator": spec.generator,
            "Nmax": t.nmax,
            "interim_information_fractions": "0.40;0.60;0.80",
            "recruitment_duration_months": t.recruitment_duration,
            "recruitment_interval_months": t.recruit_interval,
            "final_analysis_calendar_month": t.final_analysis_time,
            "minimum_administrative_followup_months": (
                t.final_analysis_time - (t.nmax - 1) * t.recruit_interval
            ),
            "independent_censoring": f"Uniform(0,{t.cens_upper:g}) months",
            "tau_RMST_final_months": t.tau_rmst,
            "tau_RMUT_final_months": t.tau_rmut,
            "RMUT_w_R": t.w_r,
        })
    return rows


def verify_archive(path: Path) -> dict:
    import csv
    config_path = PACKAGE / "data/estimands/computation_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    archive_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    assert archive_hash == config["saved_result"]["sha256"], "Saved-result SHA-256 does not match the archived configuration."
    assert config["saved_result"]["filename"] == DEFAULT_ARCHIVE.name
    script_path = Path(__file__).resolve()
    assert config["portable_script"]["path"] == script_path.relative_to(PACKAGE).as_posix()
    assert hashlib.sha256(script_path.read_bytes()).hexdigest() == config["portable_script"]["sha256"], "Portable-script SHA-256 does not match the archived configuration."
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    specs = build_specs()
    assert len(rows) == 4 and {r["case"] for r in rows} == {s.case for s in specs}
    configurations = {r["case"]: r for r in config["case_configurations"]}
    assert len(configurations) == len(config["case_configurations"]) == len(specs)
    assert config["scenario"] == "H1-standard" and config["integration_grid_points"] == 20001
    for index, (row, spec) in enumerate(zip(rows, specs)):
        assert row["case"] == spec.case and int(row["seed"]) == 2026080100 + index * 100
        assert row["version"] == spec.version and row["generator"] == spec.generator and row["status"] == spec.status
        for key in ("wr", "wr_mcse", "cox_hr", "cox_hr_mcse", "delta_rmst_months", "delta_rmut_months"):
            assert math.isfinite(float(row[key])), (spec.case, key)
        assert int(row["hr_batches"]) == 12 and int(row["hr_n_per_arm_per_batch"]) == 100000
        assert int(row["wr_batches"]) == 20 and int(row["wr_pairs_per_batch"]) == 200000
        for key, expected_value in (("tau_rmst_months", spec.trial.tau_rmst), ("tau_rmut_months", spec.trial.tau_rmut), ("rmut_weight_recurrence_state", spec.trial.w_r)):
            assert math.isclose(float(row[key]), expected_value, rel_tol=1e-12, abs_tol=1e-12), (spec.case, key)
        for key, configured_value in configurations[spec.case].items():
            assert row[key] == str(configured_value), (spec.case, key, "configuration")
        expected = (1-spec.trial.w_r)*float(row["delta_integral_sd_months"]) + spec.trial.w_r*float(row["delta_integral_sc_months"])
        assert math.isclose(float(row["delta_rmut_months"]), expected, rel_tol=1e-12, abs_tol=1e-12)
    with (PACKAGE/"data/configuration/DGM_all_scenarios.csv").open(encoding="utf-8-sig", newline="") as stream:
        current = {(r["case"],r["scenario"],r["arm"],int(r["component"])):r for r in csv.DictReader(stream)}
    count = 0
    fields = ("weight","mu_log_R","median_R_months","sigma_log_R","mu_log_D","median_D_months","sigma_log_D","rho_log_R_D","cov_log_R_D")
    for scenario in ("standard", "observed"):
        for row in parameter_rows(specs, scenario):
            target = current[(row["case"],row["scenario"],row["arm"],row["component"])]
            for field in fields:
                assert math.isclose(float(row[field]),float(target[field]),rel_tol=1e-9,abs_tol=1e-9),(row["case"],field)
                count += 1
    with (PACKAGE/"data/configuration/trial_process.csv").open(encoding="utf-8-sig", newline="") as stream:
        saved_trials = list(csv.DictReader(stream))
    trials = {r["case"]: r for r in saved_trials}
    assert len(trials) == len(saved_trials) == len(specs)
    trial_fields = 0
    for row in trial_rows(specs):
        target = trials[row["case"]]
        for field, expected_value in row.items():
            if isinstance(expected_value, (int, float)):
                assert math.isclose(float(target[field]),expected_value,rel_tol=1e-12,abs_tol=1e-12),(row["case"],field,"trial process")
            else:
                assert target[field] == str(expected_value),(row["case"],field,"trial process")
            trial_fields += 1
    return {"status":"passed","cases":4,"DGM_fields_matched":count,"trial_fields_matched":trial_fields,"archive_sha256":archive_hash,"Monte_Carlo_executed":False}


def main() -> None:
    global np, brentq, multivariate_normal, norm
    import argparse, csv
    parser = argparse.ArgumentParser(description="Verify saved population estimands, or explicitly perform their DGM integration.")
    parser.add_argument("--verify-archive",action="store_true",help="Check the saved results and DGM; this is the default action.")
    parser.add_argument("--archive",type=Path,default=DEFAULT_ARCHIVE)
    parser.add_argument("--recompute",action="store_true",help="Run Monte Carlo and numerical integration.")
    parser.add_argument("--output-dir",type=Path,help="New output directory outside the submission package, required for --recompute.")
    parser.add_argument("--case",choices=[s.case for s in build_specs()],help="Limit optional integration to one case; use that case's specified seed.")
    args = parser.parse_args()
    if not args.recompute:
        print(json.dumps(verify_archive(args.archive.resolve()),ensure_ascii=False))
        return
    if args.verify_archive or args.output_dir is None:
        parser.error("--recompute requires --output-dir and cannot be combined with --verify-archive.")
    output = args.output_dir.resolve()
    if output.exists() or output == PACKAGE or PACKAGE in output.parents:
        parser.error("Choose a new output directory outside the submission package.")
    try:
        import numpy as np
        from scipy.optimize import brentq
        from scipy.stats import multivariate_normal, norm
    except ImportError as error:
        parser.error("Install the integration dependencies in requirements.txt: " + str(error))
    output.mkdir(parents=True)
    results=[]
    for index,spec in enumerate(build_specs()):
        if args.case and spec.case != args.case: continue
        print("Integrating " + spec.case,flush=True)
        results.append(compute_estimands(spec,2026080100+index*100))
    with (output/"H1_standard_true_estimands.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(results[0]))
        writer.writeheader();writer.writerows(results)
    print("Results saved in " + str(output))


if __name__ == "__main__":
    main()
