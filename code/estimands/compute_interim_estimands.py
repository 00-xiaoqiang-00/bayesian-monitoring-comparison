"""Integrate the locked RMST and RMUT contrasts at the scheduled analysis windows."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy

import compute_h1_observed_estimands as integration


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def component_rows(rows, case, scenario, arm):
    selected = [r for r in rows if r["case"] == case and r["scenario"] == scenario and r["arm"] == arm]
    selected.sort(key=lambda r: int(r["component"]))
    if not selected:
        raise ValueError(f"Missing DGM rows for {case}, {scenario}, {arm}.")
    if not math.isclose(sum(float(r["weight"]) for r in selected), 1.0, abs_tol=1e-9):
        raise ValueError(f"Mixture weights do not sum to one: {case}, {scenario}, {arm}.")
    for row in selected:
        if float(row["weight"]) < 0 or min(float(row["sigma_log_R"]), float(row["sigma_log_D"])) <= 0:
            raise ValueError("DGM weights and log-scale standard deviations are invalid.")
        if not -1 < float(row["rho_log_R_D"]) < 1:
            raise ValueError("A nonsingular log-scale correlation is required.")
    return selected


def make_arm(rows):
    return integration.dgm.Arm(tuple(integration.dgm.Component(
        float(row["weight"]), float(row["mu_log_R"]), float(row["sigma_log_R"]),
        float(row["mu_log_D"]), float(row["sigma_log_D"]), float(row["rho_log_R_D"]),
    ) for row in rows))


def integrate_contrasts(control, treatment, tau_rmst, tau_rmut, weight):
    control_arm, treatment_arm = make_arm(control), make_arm(treatment)
    means = {}
    for arm_name, arm in (("control", control_arm), ("treatment", treatment_arm)):
        for horizon in set((tau_rmst, tau_rmut)):
            if horizon <= 0:
                raise ValueError("Restriction windows must be positive.")
            means[arm_name, horizon] = integration.arm_composite_mean(arm, horizon)
    rmst_control, rmst_control_detail = means["control", tau_rmst]
    rmst_treatment, rmst_treatment_detail = means["treatment", tau_rmst]
    composite_control, composite_control_detail = means["control", tau_rmut]
    composite_treatment, composite_treatment_detail = means["treatment", tau_rmut]
    delta_death = integration.arm_death_mean(treatment_arm, tau_rmut) - integration.arm_death_mean(control_arm, tau_rmut)
    return {
        "dRMST": rmst_treatment - rmst_control,
        "dRMUT": (1.0 - weight) * delta_death + weight * (composite_treatment - composite_control),
        "dRMST_quadrature_error_estimate_months": sum(d["estimated_quadrature_absolute_error_months"] for d in (rmst_control_detail, rmst_treatment_detail)),
        "dRMUT_quadrature_error_estimate_months": weight * sum(d["estimated_quadrature_absolute_error_months"] for d in (composite_control_detail, composite_treatment_detail)),
        "dRMST_normal_tail_bound_months": sum(d["omitted_tail_upper_bound_months"] for d in (rmst_control_detail, rmst_treatment_detail)),
        "dRMUT_normal_tail_bound_months": weight * sum(d["omitted_tail_upper_bound_months"] for d in (composite_control_detail, composite_treatment_detail)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, required=True, help="A new directory outside the repository.")
    args = parser.parse_args()
    repository = args.repository_root.resolve()
    output = args.output_dir.resolve()
    if output.exists() or output == repository or repository in output.parents:
        parser.error("Choose a new output directory outside the repository.")
    configuration = repository / "data/configuration"
    dgm_path = configuration / "DGM_all_scenarios.csv"
    schedule_path = configuration / "四案例逐次分析安排.csv"
    trial_path = configuration / "trial_process.csv"
    dgm = read_csv(dgm_path)
    schedule = read_csv(schedule_path)
    trials = {row["case"]: row for row in read_csv(trial_path)}
    cases = list(dict.fromkeys(row["case"] for row in schedule))
    results = []
    for case in cases:
        if case not in trials:
            raise ValueError(f"Missing trial-process specification for {case}.")
        weight = float(trials[case]["RMUT_w_R"])
        if not 0 <= weight <= 1:
            raise ValueError("RMUT utility weight must lie in [0,1].")
        analyses = [r for r in schedule if r["case"] == case and (r["analysis"].startswith("主分析") or r["analysis"] == "最终分析")]
        if len(analyses) != 4:
            raise ValueError(f"Expected three primary interim analyses and a final analysis for {case}.")
        for analysis in analyses:
            tau_rmst = float(analysis["RMST_tau_months"])
            tau_rmut = float(analysis["RMUT_tau_months"])
            for scenario in ("H1-standard", "H1-observed"):
                control = component_rows(dgm, case, scenario, "control")
                treatment = component_rows(dgm, case, scenario, "treatment")
                contrast = integrate_contrasts(control, treatment, tau_rmst, tau_rmut, weight)
                results.append({
                    "case": case, "analysis": analysis["analysis"],
                    "calendar_m": float(analysis["calendar_months"]),
                    "tau_rmst_months": tau_rmst, "tau_rmut_months": tau_rmut,
                    "scenario": scenario, **contrast,
                    "rmut_weight_recurrence_state": weight,
                })
        print("Integrated " + case, flush=True)
    output.mkdir(parents=True)
    result_path = output / "interim_true_estimands.csv"
    with result_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    metadata = {
        "method": "Deterministic conditional-normal restricted-mean quadrature, using the population-effect integration helpers in the same directory. The OS restricted mean has a closed form.",
        "target": "Treatment-minus-control latent-event restricted means at the scheduled method-specific windows; censoring and observed risk sets are not part of these targets.",
        "normal_integration_z_bound": integration.Z_BOUND,
        "quadrature_epsabs": integration.QUAD_ABS,
        "quadrature_epsrel": integration.QUAD_REL,
        "random_numbers_used": False,
        "runtime": {"numpy": np.__version__, "scipy": scipy.__version__},
        "input_sha256": {p.name: sha256(p) for p in (dgm_path, schedule_path, trial_path)},
        "source_sha256": {p.name: sha256(p) for p in (Path(__file__), Path(integration.__file__), Path(integration.dgm.__file__))},
        "result_sha256": {result_path.name: sha256(result_path)},
        "numerical_note": "Restricted means are evaluated by deterministic conditional-normal quadrature. Quadrature errors are estimates; omitted-normal-tail bounds are analytic. These numerical errors do not include floating-point rounding, original-data sampling, or DGM-fit uncertainty.",
    }
    (output / "computation_config.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("Results saved to " + str(output), flush=True)


if __name__ == "__main__":
    main()
