"""Verify the tradeoff between shared climate coefficients and biodiversity payments.

The Belgian certificate restricts every feasible policy to at least one climate
coefficient below 0.269. Exact interval coverage then bounds the number of
countries implementable without biodiversity payments. Complete policy
constructions are checked against every original action.
"""

from pathlib import Path
from fractions import Fraction
from itertools import combinations
from collections import defaultdict
import argparse, csv, hashlib, json, math, sqlite3
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DB_SHA = "f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d"
CUTOFF = Fraction(269, 1000)
BANDS = (0.0, 1e-6, 1e-4)


def rational(value):
    return Fraction.from_float(float(value))


def endpoint(value):
    return (
        None
        if value is None
        else Fraction(int(value["numerator"]), int(value["denominator"]))
    )


def write_csv(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def maximal_masks(masks):
    """Remove masks contained in another feasible mask, preserving coverage maxima."""
    return {
        mask: value
        for mask, value in masks.items()
        if not any(mask != other and mask & other == mask for other in masks)
    }


def constrained_coverage(regions, cutoff=CUTOFF):
    """Exact maximum coverage by q points, at least one strictly below cutoff.

    Merged interval interiors may contain isolated noncompliant breakpoints.
    They are safe outer sets. Every finite collection of covered open cells has
    interior witnesses away from all such breakpoints. Endpoint response sets
    cannot improve universal compliance over their neighboring open cells.
    """
    countries = sorted(regions)
    intervals = {
        country: [tuple(endpoint(v) for v in pair) for pair in regions[country]]
        for country in countries
    }
    points = sorted(
        {Fraction(0), cutoff}
        | {
            value
            for cells in intervals.values()
            for cell in cells
            for value in cell
            if value is not None
        }
    )
    masks, low_masks = {}, {}
    for left, right in zip(points, points[1:] + [None]):
        price = (left + right) / 2 if right is not None else left + 1
        mask = sum(
            1 << k
            for k, country in enumerate(countries)
            if any(
                lo < price and (hi is None or price < hi)
                for lo, hi in intervals[country]
            )
        )
        masks.setdefault(mask, price)
        if price < cutoff:
            low_masks.setdefault(mask, price)
    masks, low_masks = maximal_masks(masks), maximal_masks(low_masks)
    individual = sum(bool(cells) for cells in intervals.values())
    result = []
    for q in range(1, len(masks) + 2):
        best, witnesses = -1, None
        for low_mask, low_price in low_masks.items():
            for chosen in combinations(masks, q - 1):
                covered = low_mask
                for mask in chosen:
                    covered |= mask
                count = covered.bit_count()
                if count > best:
                    best, witnesses = count, [low_price] + [
                        masks[mask] for mask in chosen
                    ]
        if witnesses is None:
            break
        result.append(
            {
                "maximum_climate_coefficients": q,
                "maximum_climate_only_coverage_given_anchor": best,
                "complement_lower_bound": len(countries) - best,
                "coverage_witness_prices": [float(price) for price in witnesses],
            }
        )
        if best == individual:
            break
    return result


def anchor_certificate(e, b, c, baseline, sdu_ids, cutoff=CUTOFF):
    """Exact rational regional dominance and retained land-value upper bound.

    A strict shifted-payoff advantage, nonnegative climate-feature difference
    and nonnegative biodiversity difference imply dominance for every
    lambda >= cutoff and every beta >= 0. Retained actions form a superset of
    all actions inducible anywhere in that half-plane.
    """
    a = np.maximum(c - c[:, 0, None], 0)
    witness = np.full(e.shape, -1, dtype=np.int16)
    records, strict_advantages = [], []
    for i in range(len(e)):
        for action in range(e.shape[1]):
            candidates = np.flatnonzero((a[i] >= a[i, action]) & (b[i] >= b[i, action]))
            candidates = sorted(
                candidates,
                key=lambda other: e[i, other]
                - e[i, action]
                + float(cutoff) * (a[i, other] - a[i, action]),
                reverse=True,
            )
            for other in candidates:
                advantage = (
                    rational(e[i, other])
                    - rational(e[i, action])
                    + cutoff * (rational(a[i, other]) - rational(a[i, action]))
                )
                if advantage > 0:
                    witness[i, action] = other
                    strict_advantages.append(advantage)
                    records.append(
                        {
                            "country": "Belgium",
                            "sdu": int(sdu_ids[i]),
                            "action": action,
                            "dominating_action": int(other),
                            "strict_advantage_at_cutoff": float(advantage),
                        }
                    )
                    break
    retained = witness < 0
    if not retained.any(axis=1).all():
        raise AssertionError("A unit has no retained action")
    upper = sum((rational(max(e[i, retained[i]])) for i in range(len(e))), Fraction(0))
    target = rational(baseline["E0"])
    effective = target - Fraction(1, 100000000) * max(Fraction(1), abs(target))
    shortfall = effective - upper
    if shortfall <= 0:
        raise AssertionError("The regional exclusion certificate is not positive")
    result = {
        "country": "Belgium",
        "climate_cutoff": float(cutoff),
        "cutoff_numerator": str(cutoff.numerator),
        "cutoff_denominator": str(cutoff.denominator),
        "excluded_actions": int((witness >= 0).sum()),
        "original_actions": int(e.size),
        "minimum_strict_advantage": float(min(strict_advantages)),
        "retained_land_value_upper": float(upper),
        "land_value_floor": float(target),
        "effective_land_value_floor": float(effective),
        "certified_land_value_shortfall": float(shortfall),
        "upper_numerator": str(upper.numerator),
        "upper_denominator": str(upper.denominator),
        "shortfall_numerator": str(shortfall.numerator),
        "shortfall_denominator": str(shortfall.denominator),
        "scope": "Every finite lambda >= 269/1000 and beta >= 0; even arbitrary mixtures of retained actions",
    }
    return result, records


def response_checks(e, b, c, baseline, alpha, lam, beta):
    a = np.maximum(c - c[:, 0, None], 0)
    carbon_payment, biodiversity_payment = lam * a, beta * b
    payoff = e + carbon_payment + biodiversity_payment
    maximum = payoff.max(axis=1)
    scale = np.maximum(1, np.max(abs(e), axis=1))
    target = np.array(
        [
            float(baseline["E0"]),
            float(baseline["B0"])
            + alpha * (float(baseline["Bmax"]) - float(baseline["B0"])),
            float(baseline["C0"]),
        ]
    )
    effective = target - 1e-8 * np.maximum(1, abs(target))
    # Enlarge every response set conservatively to cover payoff rounding.
    guard = (
        16
        * np.finfo(float).eps
        * np.maximum(
            1, np.max(abs(e) + abs(carbon_payment) + abs(biodiversity_payment), axis=1)
        )
    )
    result = []
    for band in BANDS:
        admitted = (
            maximum[:, None] - payoff <= band * scale[:, None] + 2 * guard[:, None]
        )
        lower = np.array(
            [
                math.fsum(np.min(np.where(admitted, outcome, np.inf), axis=1))
                for outcome in (e, b, c)
            ]
        )
        result.append(
            {
                "response_band": band,
                "compliant": bool(np.all(lower >= effective)),
                "E_min": float(lower[0]),
                "B_min": float(lower[1]),
                "C20_min": float(lower[2]),
                "minimum_normalized_margin": float(
                    np.min((lower - target) / np.maximum(1, abs(target)))
                ),
            }
        )
    return result


def run(database, intervals_path, policies_path, baselines_path, output):
    database, output = Path(database).resolve(), Path(output)
    digest = hashlib.sha256(database.read_bytes()).hexdigest()
    if digest != DB_SHA:
        raise ValueError("The original European coefficient checksum does not match")
    output.mkdir(parents=True, exist_ok=True)
    with Path(baselines_path).open(newline="") as stream:
        baselines = {row["country"]: row for row in csv.DictReader(stream)}
    data = {}
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as connection:
        for country, baseline in baselines.items():
            raw = np.array(
                connection.execute(
                    "SELECT sdu,scenario_id,net_econ_value,biodiversity,net_ghg_co2e FROM value "
                    "WHERE country=? ORDER BY sdu,scenario_id",
                    (country,),
                ).fetchall()
            ).reshape(-1, 13, 5)
            if len(raw) != int(baseline["n_sdu"]) or not np.array_equal(
                raw[:, :, 1], np.tile(np.arange(13), (len(raw), 1))
            ):
                raise AssertionError("Source action or unit order: " + country)
            data[country] = (raw[:, :, 2], raw[:, :, 3], raw[:, :, 4], raw[:, 0, 0])
    e, b, c, ids = data["Belgium"]
    certificate, dominance_records = anchor_certificate(
        e, b, c, baselines["Belgium"], ids
    )
    (output / "climate_ceiling_certificate.json").write_text(
        json.dumps(certificate, indent=2)
    )
    write_csv(output / "climate_ceiling_witnesses.csv", dominance_records)
    interval_data = json.loads(Path(intervals_path).read_text())["regions_by_ambition"]
    lower_bounds = {}
    for alpha, regions in interval_data.items():
        for result in constrained_coverage(regions):
            lower_bounds[(float(alpha), result["maximum_climate_coefficients"])] = (
                result
            )
    with Path(policies_path).open(newline="") as stream:
        policy_rows = list(csv.DictReader(stream))
    checks, constructions = [], defaultdict(list)
    for row in policy_rows:
        alpha, q = float(row["ambition"]), int(row["maximum_climate_coefficients"])
        country, lam, beta = (
            row["country"],
            float(row["lambda_climate"]),
            float(row["beta_B"]),
        )
        if not all(math.isfinite(value) and value >= 0 for value in (lam, beta)):
            raise AssertionError("Coefficients must be finite and nonnegative")
        e, b, c, _ = data[country]
        current = response_checks(e, b, c, baselines[country], alpha, lam, beta)
        for record in current:
            checks.append(
                {
                    "ambition": alpha,
                    "maximum_climate_coefficients": q,
                    "country": country,
                    "lambda_climate": lam,
                    "beta_B": beta,
                    **record,
                }
            )
        constructions[(alpha, q)].append((country, lam, beta, current))
    frontier = []
    for key, bound in sorted(lower_bounds.items()):
        alpha, q = key
        records = constructions.get(key, [])
        if records and (
            len(records) != len(baselines)
            or {record[0] for record in records} != set(baselines)
        ):
            raise AssertionError("Incomplete country construction")
        if records and len({record[1] for record in records}) > q:
            raise AssertionError("Climate coefficient cap exceeded")
        for band_index, band in enumerate(BANDS):
            feasible_counts = []
            # A lower ambition can use any construction at a higher ambition.
            for (other_alpha, other_q), candidate in constructions.items():
                if (
                    other_alpha >= alpha
                    and other_q <= q
                    and all(record[3][band_index]["compliant"] for record in candidate)
                ):
                    feasible_counts.append(sum(record[2] > 0 for record in candidate))
            upper = min(feasible_counts) if feasible_counts else None
            low = bound["complement_lower_bound"]
            if upper is not None and upper < low:
                raise AssertionError(
                    "Constructive upper bound contradicts necessary lower bound"
                )
            frontier.append(
                {
                    "ambition": alpha,
                    "maximum_climate_coefficients": q,
                    "response_band": band,
                    "maximum_climate_only_coverage_given_anchor": bound[
                        "maximum_climate_only_coverage_given_anchor"
                    ],
                    "complement_lower_bound": low,
                    "complement_upper_bound": upper,
                    "minimum_identified": upper == low,
                }
            )
    primary = [
        row
        for row in frontier
        if row["ambition"] == 0.05 and row["response_band"] == 1e-4
    ]
    if [
        (
            row["maximum_climate_coefficients"],
            row["complement_lower_bound"],
            row["complement_upper_bound"],
        )
        for row in primary
    ] != [(q, k, k) for q, k in enumerate([20, 14, 11, 9, 7, 5, 4], 1)]:
        raise AssertionError("Primary policy frontier changed")
    write_csv(output / "price_group_frontier.csv", frontier)
    write_csv(output / "price_group_policy_checks.csv", checks)
    write_csv(output / "price_group_policies.csv", policy_rows)
    summary = {
        "source_sha256": digest,
        "countries": len(data),
        "units": sum(len(v[0]) for v in data.values()),
        "actions_per_unit": 13,
        "complete_constructions": len(constructions),
        "country_policy_rows": len(policy_rows),
        "response_checks": len(checks),
        "primary_frontier": primary,
        "minimum_primary_distinct_climate_coefficients_for_four_complements": min(
            row["maximum_climate_coefficients"]
            for row in primary
            if row["minimum_identified"]
            and row["complement_lower_bound"]
            == min(item["complement_lower_bound"] for item in primary)
        ),
        "minimum_primary_coefficient_count": min(
            row["maximum_climate_coefficients"] + row["complement_lower_bound"]
            for row in primary
        ),
        "primary_common_climate_coefficient": next(
            iter({row[1] for row in constructions[(0.05, 1)]})
        ),
        "coefficient_count_definition": "Distinct climate coefficient values plus positive national biodiversity coefficients; group assignment and administrative costs excluded",
        "primary_response_checks_pass": all(
            row["compliant"] for row in checks if row["ambition"] == 0.05
        ),
        "arithmetic": "Exact rational regional dominance and land-value separation; exact rational interval endpoints; conservatively enlarged float64 policy response sets",
        "source_unchanged": hashlib.sha256(database.read_bytes()).hexdigest() == digest,
    }
    (output / "price_group_verification.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--intervals", type=Path, required=True)
    parser.add_argument(
        "--policies", type=Path, default=ROOT / "certificates/price_groups.csv"
    )
    parser.add_argument(
        "--baselines", type=Path, default=ROOT / "certificates/Frozen_Baselines.csv"
    )
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/price_groups")
    args = parser.parse_args()
    run(args.database, args.intervals, args.policies, args.baselines, args.output)
