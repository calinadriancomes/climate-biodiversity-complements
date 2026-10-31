"""Rebuild reader tables and check the frozen numerical evidence; no optimization.

Reader calculations are separate from the original-coefficient replay in verify_matrices.py.
"""

from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
import shutil

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "certificates"
SENSITIVITY = ROOT / "certificates" / "sensitivity"
DATA = ROOT / "outputs" / "data"
RESULTS = ROOT / "outputs" / "checks"


def read(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def number(row, name):
    return float(row[name])


def true(v):
    return str(v).lower() in ("true", "1")


def total(rows, name):
    return math.fsum(number(r, name) for r in rows)


def equal(a, b, tol=0.03):
    if abs(a - b) > max(tol, 1e-10 * max(1.0, abs(a), abs(b))):
        raise AssertionError((a, b))


def save_csv(name, rows):
    if not rows:
        raise ValueError("No rows for " + name)
    with (DATA / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        idmap = {
            r["certificate_assignment"]: r["policy_id"]
            for r in read(ROOT / "inputs/policy_identifiers.csv")
        }
        for row in rows:
            row = {k: idmap.get(str(v), v) for k, v in row.items()}
            w.writerow(row)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run(output=ROOT / "outputs"):
    global DATA, RESULTS
    DATA = Path(output) / "data"
    RESULTS = Path(output) / "checks"
    DATA.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    checks = []
    for r in read(CERT / "MANIFEST_SHA256.csv"):
        p = CERT / r["path"]
        assert p.stat().st_size == int(r["bytes"]) and sha(p) == r["sha256"], r["path"]
        checks.append("source_hash:" + r["path"])
    base = read(CERT / "Frozen_Baselines.csv")
    bmap = {r["country"]: r for r in base}
    assert len(base) == 23 and sum(int(r["n_sdu"]) for r in base) == 87644
    assert sum(true(r["cee11"]) for r in base) == 11
    policies = read(CERT / "final_EU23.csv")
    pcandidates = read(CERT / "policies.csv")
    bands = read(CERT / "verified_policy_bands.csv")
    assert len(policies) == 46 and len(pcandidates) == 382 and len(bands) == 1146
    assignments = {r["assignment"]: r for r in pcandidates}
    bcheck = {(r["assignment"], float(r["band"])): r for r in bands}
    policy_names = {
        "national_final": "National climate coefficients",
        "common_final": "Common climate coefficient",
    }
    policy_reader = []
    ecology = []
    for r in policies:
        b = bmap[r["country"]]
        ref = assignments[r["assignment"]]
        for field in ["E", "B", "C20", "net_transfer", "resource_opportunity_cost"]:
            equal(float(r[field]), float(ref[field]))
        equal(float(b["Emax"]) - float(r["E"]), float(r["resource_opportunity_cost"]))
        equal(
            float(r["gross_component_disbursements"])
            - float(r["gross_component_receipts"]),
            float(r["net_transfer"]),
        )
        assert float(r["gross_component_receipts"]) == 0
        for tag, eps in [("exact", 0.0), ("stress6", 1e-6), ("stress4", 1e-4)]:
            assert true(r[tag + "_robust"]) and true(
                bcheck[(r["assignment"], eps)]["robust"]
            )
            for g, target in [
                ("E", float(b["E0"])),
                ("B", float(b["B0"]) + 0.05 * (float(b["Bmax"]) - float(b["B0"]))),
                ("C20", float(b["C0"])),
            ]:
                assert float(r[tag + "_min_" + g]) >= target - 1e-8 * max(
                    1, abs(target)
                )
        label = policy_names[r["label"]]
        policy_reader.append(
            dict(
                country=r["country"],
                group="CEE-11" if true(r["CEE11"]) else "Other 12",
                policy=label,
                climate_coefficient=float(r["lambda_climate"]),
                biodiversity_coefficient=float(r["beta_B"]),
                biodiversity_complement=int(float(r["beta_B"]) > 0),
                E_USD2015_per_year=float(r["E"]),
                C20_index_tCO2e=float(r["C20"]),
                B_national_score=float(r["B"]),
                gross_payments_USD2015_per_year=float(
                    r["gross_component_disbursements"]
                ),
                gross_receipts_USD2015_per_year=float(r["gross_component_receipts"]),
                net_transfer_USD2015_per_year=float(r["net_transfer"]),
                resource_opportunity_cost_USD2015_per_year=float(
                    r["resource_opportunity_cost"]
                ),
                assignment=r["assignment"],
            )
        )
        er = dict(
            country=r["country"],
            group="CEE-11" if true(r["CEE11"]) else "Other 12",
            policy=label,
            E_target=float(b["E0"]),
            B_target=float(b["B0"]) + 0.05 * (float(b["Bmax"]) - float(b["B0"])),
            C20_target=float(b["C0"]),
        )
        for tag in ["exact", "stress6", "stress4"]:
            for g in ["E", "B", "C20"]:
                er[tag + "_minimum_" + g] = float(r[tag + "_min_" + g])
            er[tag + "_robust"] = true(r[tag + "_robust"])
        ecology.append(er)
        checks.append("policy_and_three_bands:" + r["label"] + ":" + r["country"])
    save_csv("policies.csv", policy_reader)
    save_csv("ecological_guarantees.csv", ecology)
    for name, k in [("national_final", 4), ("common_final", 21)]:
        rr = [r for r in policies if r["label"] == name]
        assert len(rr) == 23 and sum(float(r["beta_B"]) > 0 for r in rr) == k
        if name == "common_final":
            assert all(float(r["lambda_climate"]) == 0.1 for r in rr)
        else:
            assert {r["country"] for r in rr if float(r["beta_B"]) > 0} == {
                "Belgium",
                "Finland",
                "Germany",
                "Sweden",
            }
    sens = read(SENSITIVITY / "Independent_Breadth_Sensitivity.csv")
    breadth = []
    stress_bounds = []
    for r in sens:
        alpha = float(r["alpha"])
        nat = int(r["minimum_national_B_complements"])
        low = int(r["minimum_common_B_complements_lower_bound"])
        assert nat == 23 - int(r["climate_only_national_max"])
        assert low == 23 - int(r["climate_only_common_max"])
        assert low - nat == int(r["minimum_extra_B_complements_lower_bound"])
        breadth.append(
            dict(
                alpha=alpha,
                climate_only_national_countries=23 - nat,
                maximum_common_climate_only_countries=23 - low,
                national_B_minimum=nat,
                common_B_minimum_lower=low,
                common_B_minimum_upper=21,
                extra_B_minimum_lower=low - nat,
                CEE11_national_B_minimum=0,
                CEE11_common_B_minimum_lower=int(
                    r["CEE11_minimum_common_B_complements_lower_bound"]
                ),
                CEE11_common_B_minimum_upper=11,
                choice_convention="All exact best responses",
                source="Coverage; target nesting",
            )
        )
        stress_bounds.append(
            dict(
                alpha=alpha,
                epsilon=1e-4,
                national_B_minimum_lower=nat,
                national_B_minimum_upper=4,
                common_B_minimum_lower=low,
                common_B_minimum_upper=21,
                extra_B_minimum_lower=low - 4,
                status="Logical bounds from nested response sets and the existing 5% constructions; no new price search",
            )
        )
        checks.append("breadth_accounting:" + str(alpha))
    save_csv("breadth.csv", breadth)
    save_csv("payoff_band_logical_bounds.csv", stress_bounds)
    witnesses = read(SENSITIVITY / "ClimateOnly_Alpha_Witnesses.csv")
    unimplemented = {}
    for alpha in [0.0, 0.025, 0.05]:
        wr = [
            r
            for r in witnesses
            if float(r["alpha"]) == alpha
            and r["test_kind"] == "country_individual_witness"
        ]
        assert len(wr) == next(
            r["climate_only_national_countries"] for r in breadth if r["alpha"] == alpha
        )
        assert all(true(r["robust_all_exact_best_responses"]) for r in wr)
        unimplemented[str(alpha)] = sorted(set(bmap) - {r["country"] for r in wr})
    save_csv("climate_only_witnesses.csv", witnesses)
    scopes = []
    for r in base:
        scopes.append(
            dict(
                country=r["country"],
                included=True,
                group="CEE-11" if true(r["cee11"]) else "Other 12",
                n_sdu=int(r["n_sdu"]),
                E0_USD2015_per_year=float(r["E0"]),
                C20_0_index_tCO2e=float(r["C0"]),
                B0_national_score=float(r["B0"]),
                Bmax_national_score=float(r["Bmax"]),
                Emax_USD2015_per_year=float(r["Emax"]),
                coverage_reason="All available EU members; fixed before outcome comparisons",
            )
        )
    for c in ["Cyprus", "Ireland", "Luxembourg", "Malta"]:
        scopes.append(
            dict(
                country=c,
                included=False,
                group="Not represented",
                n_sdu="",
                E0_USD2015_per_year="",
                C20_0_index_tCO2e="",
                B0_national_score="",
                Bmax_national_score="",
                Emax_USD2015_per_year="",
                coverage_reason="Absent from EU source coverage",
            )
        )
    save_csv("scope.csv", sorted(scopes, key=lambda r: r["country"]))
    accounts = read(CERT / "fixed_allocation_accounts.csv")
    part = read(CERT / "final_recentered_participation.csv")
    am = {(r["policy"], r["country"]): r for r in accounts}
    pm = {(r["policy"], r["country"]): r for r in part}
    cash = []
    cash_country = []
    for policy, label in policy_names.items():
        rr = [r for r in policies if r["label"] == policy]
        aa = [r for r in accounts if r["policy"] == policy]
        pp = [r for r in part if r["policy"] == policy]
        assert len(rr) == len(aa) == len(pp) == 23
        for a in aa:
            equal(
                float(a["level_net"]) - float(a["baseline_entitlement"]),
                float(a["signed_B_net"]),
            )
            equal(
                float(a["signed_B_gross"]) - float(a["signed_B_receipts"]),
                float(a["signed_B_net"]),
            )
            equal(
                float(a["baseline_entitlement"]),
                float(a["beta_B"]) * float(bmap[a["country"]]["B0"]),
            )
            assert float(a["fixed_allocation_choice_gap"]) == 0
            p = pm[(policy, a["country"])]
            cash_country.append(
                dict(
                    country=a["country"],
                    policy=label,
                    level_transfer=float(a["level_net"]),
                    baseline_entitlement=float(a["baseline_entitlement"]),
                    recentered_gross=float(a["signed_B_gross"]),
                    recentered_receipts=float(a["signed_B_receipts"]),
                    recentered_net=float(a["signed_B_net"]),
                    net_invoice_gross=float(a["net_invoice_payments"]),
                    net_invoice_receipts=float(a["net_invoice_receipts"]),
                    units_below_outside_option=int(
                        p["units_below_modeled_outside_option"]
                    ),
                    outside_option_shortfall=float(p["sum_payoff_shortfall"]),
                )
            )
            checks.append("fixed_allocation_accounts:" + policy + ":" + a["country"])
        cash.append(
            dict(
                policy=label,
                countries=len(rr),
                B_complements=sum(float(r["beta_B"]) > 0 for r in rr),
                E=total(rr, "E"),
                C20=total(rr, "C20"),
                resource_opportunity_cost=total(rr, "resource_opportunity_cost"),
                level_gross=total(rr, "gross_component_disbursements"),
                level_receipts=0.0,
                level_net=total(rr, "net_transfer"),
                baseline_entitlement=total(aa, "baseline_entitlement"),
                recentered_gross=total(aa, "signed_B_gross"),
                recentered_receipts=total(aa, "signed_B_receipts"),
                recentered_net=total(aa, "signed_B_net"),
                IR_affected_SDUs=sum(
                    int(r["units_below_modeled_outside_option"]) for r in pp
                ),
                IR_payoff_shortfall=total(pp, "sum_payoff_shortfall"),
            )
        )
    save_csv("financial_accounts.csv", cash_country)
    save_csv("aggregate_accounts.csv", cash)
    groups = read(CERT / "group_accounts.csv")
    for r in groups:
        rr = [
            x
            for x in policies
            if x["label"] == r["policy"]
            and (r["group"] == "EU23" or (r["group"] == "CEE11") == true(x["CEE11"]))
        ]
        for f, field in [
            ("E", "E"),
            ("C20", "C20"),
            ("gross_payments", "gross_component_disbursements"),
            ("net_transfer", "net_transfer"),
        ]:
            equal(float(r[f]), total(rr, field))
    save_csv("subgroup_accounts.csv", groups)
    matched = read(CERT / "verified_matched.csv")
    finalmatched = [
        r for r in matched if r["policy"] in ("national_final", "common_final")
    ]
    assert len(matched) == 19 and len(finalmatched) == 8
    for r in matched:
        lo, hi = float(r["resource_gain_lower"]), float(r["resource_gain_upper"])
        assert lo <= hi
        emax = float(bmap[r["country"]]["Emax"])
        equal(emax - float(r["resource_upper"]) - float(r["reference_E"]), lo)
        equal(emax - float(r["resource_lower"]) - float(r["reference_E"]), hi)
        equal(float(r["planner_E"]) - float(r["reference_E"]), lo)
        for g in ("B", "C20"):
            assert float(r["planner_" + g]) >= float(r["reference_" + g]) - 1e-8 * max(
                1, abs(float(r["reference_" + g]))
            )
        checks.append("matched_resource_ledger:" + r["country"] + ":" + r["policy"])
    save_csv("resource_bounds.csv", matched)
    exclusions = read(SENSITIVITY / "Independent_Sensitivity.csv")
    save_csv("exclusion_sensitivity.csv", exclusions)
    dom = read(CERT / "dominance_summary.csv")
    dd = []
    for r in dom:
        n = int(bmap[r["country"]]["n_sdu"])
        assert int(r["actions_removed"]) + int(r["actions_retained"]) == 13 * n
        dd.append(
            dict(
                country=r["country"],
                climate_base=r["climate"],
                biodiversity_base=r["bio"],
                SDUs=n,
                dominated_actions=int(r["actions_removed"]),
                dominated_baseline_actions=int(r["baseline_actions_removed"]),
                baseline_actions_removed_share=int(r["baseline_actions_removed"]) / n,
            )
        )
    save_csv("dominance_diagnostics.csv", dd)
    save_csv("alternative_policies.csv", read(CERT / "alternative_best_witnesses.csv"))
    save_csv("price_regions.csv", read(CERT / "robust_carbon_outer.csv"))
    save_csv("robustness_history.csv", read(CERT / "old_and_new_robustness.csv"))
    # Identifying an interval endpoint is different from proving an exact feasibility boundary.
    fin = [r for r in exclusions if r["country"] == "Finland"]
    f0 = float(fin[0]["dual_certified_shortfall_lower"])
    f5 = float(fin[-1]["dual_certified_shortfall_lower"])
    cutoff = -f0 / ((f5 - f0) / 0.05)
    facts = dict(
        population=23,
        CEE11=11,
        SDUs=87644,
        actions=13,
        action_rows=87644 * 13,
        source_certificate_files=len(read(CERT / "MANIFEST_SHA256.csv")),
        stored_policy_candidates=382,
        stored_band_checks=1146,
        main_policy_rows=46,
        breadth=breadth,
        stress_logical_bounds=stress_bounds,
        climate_only_exclusions_by_alpha=unimplemented,
        financial_aggregates=cash,
        finland_same_certificate_alpha_threshold=cutoff,
        matched_final=finalmatched,
        verification_layer="Reader tables and numerical accounts; original-coefficient verification is a separate command",
        source_matrix_sha256="f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d",
    )
    (DATA / "facts.json").write_text(json.dumps(facts, indent=2), encoding="utf-8")
    (RESULTS / "consolidation_checks.json").write_text(
        json.dumps(
            dict(
                status="PASS",
                checks=len(checks),
                checks_performed=checks,
                large_optimizations_rerun=0,
                raw_matrix_replay_this_run=False,
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    shutil.copytree(
        ROOT / "data/price_groups", DATA / "price_groups", dirs_exist_ok=True
    )
    shutil.copytree(
        ROOT / "data/fixed_groups", DATA / "fixed_groups", dirs_exist_ok=True
    )
    print(
        json.dumps(
            dict(
                status="PASS",
                checks=len(checks),
                national_excluded=unimplemented,
                accounts=cash,
            ),
            indent=2,
        )
    )
    return facts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs")
    run(parser.parse_args().output)
