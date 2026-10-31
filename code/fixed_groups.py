"""Certify the minimum complement count for two fixed country groups.

CEE-11 countries share one climate coefficient and the other twelve countries
share another. All original national safeguards and payment rules are retained.
The two group coefficients may coincide, so this restriction allows at most two
distinct climate coefficients across all countries.
"""

from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sqlite3
import numpy as np

from price_groups import endpoint, response_checks, write_csv

ROOT = Path(__file__).resolve().parents[1]
SHA = "f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d"


def regional_certificate(country, e, b, c, baseline, cutoff):
    """Exact strict dominance gives an upper bound for every admissible beta."""
    a = np.maximum(c - c[:, [0]], 0.0)
    upper, removed, advantage = F(0), 0, None
    for i in range(len(e)):
        values = [(F(float(e[i, s])), F(float(b[i, s])), F(float(a[i, s])))
                  for s in range(e.shape[1])]
        retained = []
        for s, (es, bs, az) in enumerate(values):
            witness = None
            for t, (et, bt, at) in enumerate(values):
                gap = et - es + cutoff * (at - az)
                if at >= az and bt >= bs and gap > 0:
                    witness = gap
                    break
            if witness is None:
                retained.append(es)
            else:
                removed += 1
                advantage = witness if advantage is None else min(advantage, witness)
        if not retained:
            raise AssertionError("Every unit must retain at least one action")
        upper += max(retained)
    floor = F(float(baseline["E0"]))
    effective = floor - F(1, 100000000) * max(F(1), abs(floor))
    if not upper < effective:
        raise AssertionError("The all-coefficient exclusion is not established")
    return {
        "country": country,
        "climate_cutoff": float(cutoff),
        "climate_cutoff_exact": str(cutoff),
        "removed_actions": removed,
        "retained_actions": int(e.size) - removed,
        "minimum_strict_advantage": float(advantage),
        "land_value_upper": float(upper),
        "land_value_upper_exact": str(upper),
        "original_land_value_floor": float(floor),
        "effective_land_value_floor": float(effective),
        "certified_shortfall": float(effective - upper),
        "certified_shortfall_exact": str(effective - upper),
        "scope": "Every climate coefficient at least cutoff and every nonnegative biodiversity coefficient",
    }


def group_coverage(regions, countries, cutoff):
    """Maximum climate-only coverage strictly below a required group cutoff."""
    spans = {country: [(endpoint(lo), endpoint(hi)) for lo, hi in regions[country]]
             for country in countries}
    breaks = sorted({F(0), cutoff} | {
        v for cells in spans.values() for pair in cells for v in pair
        if v is not None and 0 < v < cutoff
    })
    best, best_countries = -1, []
    for left, right in zip(breaks, breaks[1:]):
        price = (left + right) / 2
        covered = [country for country, cells in spans.items()
                   if any(lo < price and (hi is None or price < hi) for lo, hi in cells)]
        if len(covered) > best:
            best, best_countries = len(covered), sorted(covered)
    return {"countries": sorted(countries), "maximum_climate_only_coverage": best,
            "coverage_witness_countries": best_countries,
            "climate_cutoff_exact": str(cutoff)}


def run(database, intervals, output):
    database, output = Path(database).resolve(), Path(output).resolve()
    if hashlib.sha256(database.read_bytes()).hexdigest() != SHA:
        raise ValueError("The European coefficient database differs")
    output.mkdir(parents=True, exist_ok=True)
    with (ROOT / "certificates/Frozen_Baselines.csv").open(newline="") as stream:
        base = {r["country"]: r for r in csv.DictReader(stream)}
    regions = json.loads(Path(intervals).read_text())["regions_by_ambition"]["0.05"]
    groups = {
        "CEE-11": sorted(c for c, r in base.items() if r["cee11"] == "True"),
        "Other twelve countries": sorted(c for c, r in base.items() if r["cee11"] == "False"),
    }
    if [len(v) for v in groups.values()] != [11, 12]:
        raise AssertionError("The fixed country partition differs")
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as connection:
        arrays = {}
        for country, baseline in base.items():
            data = np.array(connection.execute(
                "SELECT net_econ_value,biodiversity,net_ghg_co2e FROM value "
                "WHERE country=? ORDER BY sdu,scenario_id", (country,)
            ).fetchall()).reshape(-1, 13, 3)
            if len(data) != int(baseline["n_sdu"]):
                raise AssertionError("Source unit count differs")
            arrays[country] = tuple(data[:, :, k] for k in range(3))
    slovenian_cutoff = endpoint(regions["Romania"][0][0])
    cutoffs = {"CEE-11": slovenian_cutoff, "Other twelve countries": F(269, 1000)}
    certificates = []
    coverage = []
    for group, anchor in [("CEE-11", "Slovenia"), ("Other twelve countries", "Belgium")]:
        certificates.append({"group": group, **regional_certificate(
            anchor, *arrays[anchor], base[anchor], cutoffs[group])})
        coverage.append({"group": group, **group_coverage(regions, groups[group], cutoffs[group])})
    lower = len(base) - sum(r["maximum_climate_only_coverage"] for r in coverage)
    with (ROOT / "certificates/price_groups.csv").open(newline="") as stream:
        common = {r["country"]: r for r in csv.DictReader(stream)
                  if float(r["ambition"]) == 0.05 and int(r["maximum_climate_coefficients"]) == 1}
    with (ROOT / "certificates/safeguard_frontier_policies.csv").open(newline="") as stream:
        cee = {r["country"]: r for r in csv.DictReader(stream) if int(r["q"]) == 1}
    policies, checks = [], []
    for country, baseline in base.items():
        group = "CEE-11" if baseline["cee11"] == "True" else "Other twelve countries"
        row = (cee if group == "CEE-11" else common)[country]
        lam, beta = float(row["lambda_climate"]), float(row["beta_B"])
        if not all(np.isfinite(v) and v >= 0 for v in (lam, beta)):
            raise AssertionError("Policy coefficients must be finite and nonnegative")
        policies.append({"country": country, "group": group, "lambda_climate": lam,
                         "beta_B": beta, "biodiversity_complement": beta > 0,
                         "E_floor": float(baseline["E0"]),
                         "B_floor": float(baseline["B0"]) + 0.05 * (float(baseline["Bmax"]) - float(baseline["B0"])),
                         "C20_floor": float(baseline["C0"])})
        for check in response_checks(*arrays[country], baseline, 0.05, lam, beta):
            checks.append({"country": country, "group": group, **check})
    if len(policies) != 23 or len(checks) != 69 or not all(r["compliant"] for r in checks):
        raise AssertionError("A complete robust construction has not been established")
    for group in groups:
        if len({r["lambda_climate"] for r in policies if r["group"] == group}) != 1:
            raise AssertionError("The fixed group uses multiple climate coefficients")
    upper = sum(r["biodiversity_complement"] for r in policies)
    if lower != upper:
        raise AssertionError("The exact fixed-group minimum is unresolved")
    result = {
        "source_sha256": SHA,
        "biodiversity_ambition": 0.05,
        "fixed_groups": groups,
        "regional_certificates": certificates,
        "coverage_bounds": coverage,
        "complement_lower_bound": lower,
        "complement_upper_bound": upper,
        "minimum_identified": lower == upper,
        "distinct_climate_coefficients": len({r["lambda_climate"] for r in policies}),
        "all_original_national_floors_retained": True,
        "country_policies": len(policies),
        "response_checks": len(checks),
        "all_response_checks_pass": all(r["compliant"] for r in checks),
        "minimum_normalized_margin": min(r["minimum_normalized_margin"] for r in checks),
        "source_unchanged": hashlib.sha256(database.read_bytes()).hexdigest() == SHA,
        "interpretation": "A descriptive fixed-partition restriction, not an optimal geographic grouping or a causal regional comparison",
    }
    write_csv(output / "fixed_group_policies.csv", policies)
    write_csv(output / "fixed_group_response_checks.csv", checks)
    (output / "fixed_group_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--intervals", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/fixed_groups")
    args = parser.parse_args()
    run(args.database, args.intervals, args.output)
