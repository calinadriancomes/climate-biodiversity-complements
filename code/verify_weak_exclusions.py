"""Check primary climate-only exclusions with optimistic outcomes at every tie.

Separate coordinate maxima are necessary, not sufficient, for weak compliance.
Their failure at all prices therefore excludes even favorable action selection.
All finite breakpoints, isolated ties, zero and the unbounded tail are examined.
"""

from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse, csv, hashlib, json, sqlite3, time
import numpy as np
from verify_breadth import upper_envelope, rational

ROOT = Path(__file__).resolve().parents[1]
DB_SHA = "f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d"


def run(database, output):
    t0 = time.time()
    database = Path(database).resolve()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    assert hashlib.sha256(database.read_bytes()).hexdigest() == DB_SHA
    with (ROOT / "certificates/Frozen_Baselines.csv").open(newline="") as f:
        base = {r["country"]: r for r in csv.DictReader(f)}
    db = sqlite3.connect(database.as_uri() + "?mode=ro", uri=True)
    reports = []
    for country in ["Belgium", "Finland", "Germany", "Sweden"]:
        row = base[country]
        x = np.array(
            db.execute(
                "SELECT net_econ_value,biodiversity,net_ghg_co2e FROM value WHERE country=? ORDER BY sdu,scenario_id",
                (country,),
            ).fetchall()
        ).reshape(-1, 13, 3)
        e, b, c = x[:, :, 0], x[:, :, 1], x[:, :, 2]
        a = np.maximum(c - c[:, 0, None], 0)
        ef = [[rational(v) for v in rr] for rr in e]
        af = [[rational(v) for v in rr] for rr in a]
        bf = [[rational(v) for v in rr] for rr in b]
        cf = [[rational(v) for v in rr] for rr in c]
        g = [ef, bf, cf]
        target = np.array(
            [
                float(row["E0"]),
                float(row["B0"]) + 0.05 * (float(row["Bmax"]) - float(row["B0"])),
                float(row["C0"]),
            ]
        )
        effective = [rational(v) for v in target - 1e-8 * np.maximum(1, abs(target))]
        scales = [rational(max(1, abs(v))) for v in target]
        events = defaultdict(list)
        state = []
        totals = [F(0), F(0), F(0)]
        at_zero = [F(0), F(0), F(0)]

        def maxima(i, actions):
            return tuple(max(g[k][i][s] for s in actions) for k in range(3))

        for i in range(len(e)):
            hull = upper_envelope(e[i], a[i])
            v = maxima(i, hull[0][1])
            state.append(v)
            for k in range(3):
                totals[k] += v[k]
            q0 = [s for s in range(13) if ef[i][s] == max(ef[i])]
            v0 = maxima(i, q0)
            for k in range(3):
                at_zero[k] += v0[k]
            for price, actions in hull[1:]:
                events[price].append((i, actions))
        possible = 0
        tested = 0
        best_margin = -float("inf")
        best_location = None

        def examine(values, location):
            nonlocal possible, tested, best_margin, best_location
            tested += 1
            valid = all(v >= t for v, t in zip(values, effective))
            possible += int(valid)
            margin = min(
                float((v - t) / s) for v, t, s in zip(values, effective, scales)
            )
            if margin > best_margin:
                best_margin = margin
                best_location = location

        examine(at_zero, "zero")
        for price, changes in sorted(events.items()):
            examine(totals, "open cell before " + str(float(price)))
            point = totals.copy()
            for i, new_actions in changes:
                u = [ev + price * av for ev, av in zip(ef[i], af[i])]
                top = max(u)
                all_ties = [s for s, v in enumerate(u) if v == top]
                vp = maxima(i, all_ties)
                for k in range(3):
                    point[k] += vp[k] - state[i][k]
            examine(point, "breakpoint " + str(float(price)))
            for i, actions in changes:
                vn = maxima(i, actions)
                for k in range(3):
                    totals[k] += vn[k] - state[i][k]
                state[i] = vn
        examine(totals, "unbounded tail")
        assert possible == 0, (country, possible, best_margin, best_location)
        reports.append(
            {
                "country": country,
                "alpha": 0.05,
                "optimistic_compliant_price_states": possible,
                "price_states_checked": tested,
                "best_optimistic_minimum_normalized_margin": best_margin,
                "best_location": best_location,
                "weak_climate_only_implementation_excluded": True,
            }
        )
        print("Verified optimistic climate-only exclusion:", country, flush=True)
    db.close()
    report = {
        "status": "PASS",
        "source_sha256": DB_SHA,
        "source_unchanged": hashlib.sha256(database.read_bytes()).hexdigest() == DB_SHA,
        "arithmetic": "Exact rational full affine payoff sets and optimistic coordinate maxima; effective target allowance 1e-8",
        "includes_isolated_breakpoint_ties": True,
        "optimizer_calls": 0,
        "seconds": time.time() - t0,
        "countries": reports,
    }
    (output / "weak_climate_only_exclusions.json").write_text(
        json.dumps(report, indent=2)
    )
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--database", type=Path, required=True)
    ap.add_argument("--output", type=Path, default=ROOT / "outputs/matrix_replay")
    args = ap.parse_args()
    run(args.database, args.output)
