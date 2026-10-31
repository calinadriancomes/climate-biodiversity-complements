"""Verify target-sensitive exclusions, participation and coefficient neighborhoods."""

from pathlib import Path
from fractions import Fraction
import argparse, csv, hashlib, json, math, sqlite3
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DB_SHA = "f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d"


def total(a):
    return math.fsum(np.asarray(a).ravel().tolist())


def records(p):
    with p.open(newline="") as f:
        return list(csv.DictReader(f))


def close(a, b):
    assert abs(a - b) <= max(0.02, 5e-11 * max(1, abs(a), abs(b))), (a, b)


def run(database, output):
    database = Path(database).resolve()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    assert hashlib.sha256(database.read_bytes()).hexdigest() == DB_SHA
    db = sqlite3.connect(database.as_uri() + "?mode=ro", uri=True)
    base = {
        r["country"]: r for r in records(ROOT / "certificates/Frozen_Baselines.csv")
    }
    policies = records(ROOT / "certificates/final_EU23.csv")
    reported = {
        (r["policy"], r["country"]): r
        for r in records(ROOT / "certificates/final_recentered_participation.csv")
    }
    certs = json.loads((ROOT / "certificates/screen_bounds.json").read_text())
    dom = np.load(ROOT / "certificates/dominance_witnesses.npz", allow_pickle=False)
    sensitivity = []
    screen100 = []
    participation = []
    radii = []
    measurement = []
    supplied_sensitivity = {
        (r["country"], float(r["biodiversity_ambition_fraction"])): r
        for r in records(ROOT / "certificates/sensitivity/Independent_Sensitivity.csv")
    }
    for country, row in base.items():
        x = np.array(
            db.execute(
                "SELECT net_econ_value,biodiversity,net_ghg_co2e,cropland_value,grazing_value,forestry_value,transition_cost,carbon,grazing_methane FROM value WHERE country=? ORDER BY sdu,scenario_id",
                (country,),
            ).fetchall()
        ).reshape(-1, 13, 9)
        e, b, c = x[:, :, 0], x[:, :, 1], x[:, :, 2]
        a = np.maximum(c - c[:, 0, None], 0)
        h = np.maximum(b - b[:, 0, None], 0)
        g = np.stack([e, b, c])
        q = np.maximum(1, np.max(abs(e), axis=1))
        recomputed_e = x[:, :, 3] + x[:, :, 4] + x[:, :, 5] - x[:, :, 6]
        recomputed_c = (44 / 12) * x[:, :, 7] - 0.02 * x[:, :, 8]
        error_e = float(np.max(abs(e - recomputed_e)))
        error_c = float(np.max(abs(c - recomputed_c)))
        assert np.allclose(e, recomputed_e, rtol=1e-12, atol=1e-7)
        assert np.allclose(c, recomputed_c, rtol=1e-12, atol=1e-7)
        configuration = db.execute(
            "SELECT contents FROM config WHERE country=?", (country,)
        ).fetchone()[0]
        menu = configuration.split("optimization_scenarios:", 1)[1]
        source_names = [
            line.strip().removeprefix("-").strip().strip('"')
            for line in menu.splitlines()
            if line.strip().startswith("-")
        ]
        assert set(source_names) == {
            r[0] for r in db.execute("SELECT name FROM scenario")
        }
        measurement.append(
            {
                "country": country,
                "unit_action_records": len(e) * 13,
                "maximum_E_reconstruction_error": error_e,
                "maximum_C20_reconstruction_error": error_c,
                "original_configuration_actions": len(source_names),
            }
        )
        assert b.min() >= 0
        for p in [r for r in policies if r["country"] == country]:
            lam, beta = float(p["lambda_climate"]), float(p["beta_B"])
            payment = lam * a + beta * b
            assert payment.min() >= 0
            umax = np.max(e + payment, axis=1)
            outside = np.max(e, axis=1)
            assert np.all(umax >= outside)
            recentered = np.max(e + lam * a + beta * (b - b[:, 0, None]), axis=1)
            shortfall = np.maximum(outside - recentered, 0)
            rr = reported[(p["label"], country)]
            count = int(np.sum(shortfall > 0))
            claimed = int(rr["units_below_modeled_outside_option"])
            assert count == claimed, (country, p["label"], count, claimed)
            close(total(shortfall), float(rr["sum_payoff_shortfall"]))
            close(float(np.max(shortfall)), float(rr["maximum_unit_shortfall"]))
            indifferent = umax == outside
            outside_actions = e == outside[:, None]
            assert np.all(
                (~indifferent[:, None])
                | (~outside_actions)
                | (e + payment == umax[:, None])
            )
            participation.append(
                {
                    "country": country,
                    "policy": p["label"],
                    "minimum_program_surplus_over_outside": float(
                        np.min(umax - outside)
                    ),
                    "indifferent_units": int(indifferent.sum()),
                    "recentered_units_below_outside_option": count,
                    "recentered_payoff_shortfall": total(shortfall),
                }
            )
        diameter = np.ptp(a, axis=1) + np.ptp(b, axis=1)
        radius = np.divide(
            1e-4 * q, diameter, out=np.full(len(q), np.inf), where=diameter > 0
        )
        k = int(np.argmin(radius))
        radii.append(
            {
                "country": country,
                "radius": float(radius[k]),
                "limiting_unit_row": k,
                "q": float(q[k]),
                "paid_feature_L1_diameter_bound": float(diameter[k]),
            }
        )
        if country in ["Finland", "Sweden"]:
            cert = next(
                r
                for r in certs
                if r["country"] == country
                and r["climate"] == "positive"
                and r["bio"] == "gain"
            )
            mu = np.array(cert["mu"])
            sc = np.array(cert["scales"])
            keep = dom[country + "__positive__gain"] < 0
            weighted = np.einsum("k,kij->ij", mu / sc, g)
            bound = total(np.max(weighted, where=keep, initial=-np.inf, axis=1))
            exact_weights = [
                Fraction.from_float(float(m)) / Fraction.from_float(float(d))
                for m, d in zip(mu, sc)
            ]
            active_weights = [k for k, w in enumerate(exact_weights) if w]
            exact_bound = Fraction(0)
            for i in range(len(e)):
                exact_bound += max(
                    sum(
                        (
                            exact_weights[k] * Fraction.from_float(float(g[k, i, s]))
                            for k in active_weights
                        ),
                        Fraction(0),
                    )
                    for s in np.flatnonzero(keep[i])
                )
            for alpha in [0.0, 0.025, 0.05]:
                target = np.array(
                    [
                        float(row["E0"]),
                        float(row["B0"])
                        + alpha * (float(row["Bmax"]) - float(row["B0"])),
                        float(row["C0"]),
                    ]
                )
                eff = target - 1e-8 * np.maximum(1, abs(target))
                margin = float(np.dot(mu, eff / sc) - bound) - 2e-9
                assert (
                    abs(
                        margin
                        - float(
                            supplied_sensitivity[(country, alpha)][
                                "dual_certified_shortfall_lower"
                            ]
                        )
                    )
                    < 2e-12
                )
                exact_margin = (
                    sum(
                        (
                            w * Fraction.from_float(float(t))
                            for w, t in zip(exact_weights, eff)
                        ),
                        Fraction(0),
                    )
                    - exact_bound
                    - Fraction.from_float(2e-9)
                )
                assert (exact_margin > 0) == (margin > 0) and abs(
                    float(exact_margin) - margin
                ) < 2e-12
                sensitivity.append(
                    {
                        "country": country,
                        "biodiversity_ambition_fraction": alpha,
                        "dual_certified_shortfall_lower": margin,
                        "interpretation": (
                            "ALL_FINITE_NONNEGATIVE_PRICES_EXCLUDED"
                            if margin > 0
                            else "INCONCLUSIVE"
                        ),
                        "exact_rational_margin": float(exact_margin),
                        "exact_rational_numerator": str(exact_margin.numerator),
                        "exact_rational_denominator": str(exact_margin.denominator),
                    }
                )
            witness = dom[country + "__positive__gain"]
            advantage = np.take_along_axis(e, np.maximum(witness, 0), axis=1) - e
            # Relax the delivered sufficient screen: return any removed action
            # whose recorded witness improves land value by at most USD100.
            # A different witness might still dominate it, so this is generous.
            keep100 = (witness < 0) | (advantage <= 100)
            target = np.array(
                [
                    float(row["E0"]),
                    float(row["B0"]) + 0.05 * (float(row["Bmax"]) - float(row["B0"])),
                    float(row["C0"]),
                ]
            )
            eff = target - 1e-8 * np.maximum(1, abs(target))
            margin = (
                float(
                    np.dot(mu, eff / sc)
                    - total(np.max(weighted, where=keep100, initial=-np.inf, axis=1))
                )
                - 2e-9
            )
            assert margin > 0
            screen100.append(
                {
                    "country": country,
                    "minimum_land_value_advantage_USD2015": 100,
                    "retained_actions": int(keep100.sum()),
                    "margin": margin,
                }
            )
    db.close()
    fin = [r for r in sensitivity if r["country"] == "Finland"]
    a0 = fin[0]["dual_certified_shortfall_lower"]
    a5 = fin[-1]["dual_certified_shortfall_lower"]
    crossing = -a0 / ((a5 - a0) / 0.05)
    with (output / "verified_exclusion_sensitivity.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(sensitivity[0]))
        w.writeheader()
        w.writerows(sensitivity)
    assert 4.51e-10 <= min(r["radius"] for r in radii)
    summary = {
        "status": "PASS",
        "source_sha256": DB_SHA,
        "source_unchanged": hashlib.sha256(database.read_bytes()).hexdigest() == DB_SHA,
        "exclusion_target_tolerance": 1e-8,
        "exclusion_normalized_arithmetic_allowance": 2e-9,
        "exclusion_sensitivity": sensitivity,
        "finland_same_certificate_crossing": crossing,
        "land_value_advantage_screen_100": screen100,
        "participation": participation,
        "source_measurement_checks": measurement,
        "coefficient_neighborhood": {
            "epsilon": 1e-4,
            "coordinate_convention": "Separate absolute changes in climate and national biodiversity coefficients, each bounded by r in their original numeric units, restricted to nonnegative coefficients; zero biodiversity coefficients stay zero when preserving instrument counts",
            "sufficient_condition": "r * (range_s a_is + range_s B_is) <= epsilon * q_i for every unit",
            "global_radius": min(r["radius"] for r in radii),
            "reported_conservative_radius": 4.51e-10,
            "country_bounds": radii,
        },
    }
    (output / "additional_verification.json").write_text(json.dumps(summary, indent=2))
    print(
        json.dumps(
            {
                k: v
                for k, v in summary.items()
                if k not in ["participation", "coefficient_neighborhood"]
            },
            indent=2,
        )
    )
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--database", type=Path, required=True)
    ap.add_argument("--output", type=Path, default=ROOT / "outputs/matrix_replay")
    args = ap.parse_args()
    run(args.database, args.output)
