"""Exact fixed-climate-coefficient economic ceiling for Belgium.

Sweep every open biodiversity-price cell under universal best response, retaining
all persistent affine ties. Maximize guaranteed original land value subject to
the original biodiversity and climate safeguards. Every breakpoint has a weakly
larger response set, so cannot improve this objective while meeting safeguards.
"""

from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import sqlite3, csv, json, math, hashlib, argparse
import numpy as np

F = lambda x: Fraction.from_float(float(x))


def envelope(intercepts, slopes):
    lines = {}
    for s, (e, b) in enumerate(zip(intercepts, slopes)):
        if b not in lines or e > lines[b][0]:
            lines[b] = (e, [s])
        elif e == lines[b][0]:
            lines[b][1].append(s)
    stack = []
    starts = []
    for b, (e, ss) in sorted(lines.items()):
        while stack:
            cross = (stack[-1][1] - e) / (b - stack[-1][0])
            if len(stack) == 1 or cross > starts[-1]:
                break
            stack.pop()
            starts.pop()
        starts.append(None if not stack else cross)
        stack.append((b, e, ss))
    first = 0
    for k in range(1, len(stack)):
        if starts[k] <= 0:
            first = k
        else:
            break
    return [(Fraction(0), stack[first][2])] + [
        (starts[k], stack[k][2]) for k in range(first + 1, len(stack))
    ]


def run(database, baselines, climate_witnesses, output):
    database = Path(database).resolve()
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    assert (
        hashlib.sha256(database.read_bytes()).hexdigest()
        == "f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d"
    )
    r = next(
        row
        for row in csv.DictReader(Path(baselines).open())
        if row["country"] == "Belgium"
    )
    lam_string = next(
        row["lambda"]
        for row in csv.DictReader(Path(climate_witnesses).open())
        if float(row["alpha"]) == 0.05 and row["test_kind"] == "common_witness"
    )
    lam = F(lam_string)
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as db:
        x = np.array(
            db.execute(
                'select net_econ_value,biodiversity,net_ghg_co2e from value where country="Belgium" order by sdu,scenario_id'
            ).fetchall()
        ).reshape(-1, 13, 3)
    e, b, c = x[:, :, 0], x[:, :, 1], x[:, :, 2]
    a = np.maximum(c - c[:, 0, None], 0)
    targets = [
        F(r["E0"]),
        F(float(r["B0"]) + 0.05 * (float(r["Bmax"]) - float(r["B0"]))),
        F(r["C0"]),
    ]
    effective = [t - Fraction(1, 100000000) * max(Fraction(1), abs(t)) for t in targets]
    events = defaultdict(lambda: [Fraction(0)] * 3)
    total = [Fraction(0)] * 3
    segments = 0
    for i in range(len(e)):
        intercepts = [F(e[i, s]) + lam * F(a[i, s]) for s in range(13)]
        slopes = [F(b[i, s]) for s in range(13)]
        cells = envelope(intercepts, slopes)
        segments += len(cells)
        prev = [min(F(g[i, s]) for s in cells[0][1]) for g in (e, b, c)]
        total = [u + v for u, v in zip(total, prev)]
        for beta, ss in cells[1:]:
            values = [min(F(g[i, s]) for s in ss) for g in (e, b, c)]
            events[beta] = [u + v - w for u, v, w in zip(events[beta], values, prev)]
            prev = values
    best = None
    feasible = 0
    left = Fraction(0)
    for right in sorted(events) + [None]:
        if total[1] >= effective[1] and total[2] >= effective[2]:
            feasible += 1
            if best is None or total[0] > best["values"][0]:
                best = {"left": left, "right": right, "values": total.copy()}
        if right is not None:
            total = [u + v for u, v in zip(total, events[right])]
            left = right
    assert best is not None
    beta = (
        (best["left"] + best["right"]) / 2
        if best["right"] is not None
        else best["left"] + 1
    )
    relax = effective[0] - best["values"][0]
    beta_f = float(beta)
    u = e + float(lam) * a + beta_f * b
    mx = u.max(axis=1)
    q = np.maximum(1, abs(e).max(axis=1))
    check = []
    for eps in [0.0, 1e-6, 1e-4]:
        mask = (u == mx[:, None]) if eps == 0 else mx[:, None] - u <= eps * q[:, None]
        vals = [math.fsum(np.min(np.where(mask, g, np.inf), axis=1)) for g in (e, b, c)]
        check.append(
            {
                "band": eps,
                "minimum_E": vals[0],
                "minimum_B": vals[1],
                "minimum_C20": vals[2],
                "B_and_C_safeguards": vals[1] >= float(effective[1])
                and vals[2] >= float(effective[2]),
            }
        )
    result = {
        "country": "Belgium",
        "climate_coefficient": float(lam),
        "biodiversity_ambition": 0.05,
        "original_E_floor": float(targets[0]),
        "effective_E_floor": float(effective[0]),
        "maximum_guaranteed_E_subject_to_B_C": float(best["values"][0]),
        "minimum_required_E_floor_relaxation": float(relax),
        "fraction_of_original_E_floor": float(relax / targets[0]),
        "biodiversity_coefficient_witness": beta_f,
        "beta_cell_left": float(best["left"]),
        "beta_cell_right": None if best["right"] is None else float(best["right"]),
        "attained_B_guarantee": float(best["values"][1]),
        "attained_C20_guarantee": float(best["values"][2]),
        "B_target_effective": float(effective[1]),
        "C20_target_effective": float(effective[2]),
        "exact_E_numerator": str(best["values"][0].numerator),
        "exact_E_denominator": str(best["values"][0].denominator),
        "exact_relaxation_numerator": str(relax.numerator),
        "exact_relaxation_denominator": str(relax.denominator),
        "unit_biodiversity_price_segments": segments,
        "national_breakpoints": len(events),
        "B_C_compliant_open_cells": feasible,
        "float64_witness_checks": check,
        "interpretation": "Global maximum of guaranteed original E at this fixed climate coefficient under exact best response, subject to B and C safeguards; not a welfare optimum or a comparison of matched complete EU policies.",
    }
    (out / "belgium_economic_ceiling.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--database", type=Path, required=True)
    ap.add_argument("--baselines", type=Path, required=True)
    ap.add_argument("--climate-witnesses", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    run(args.database, args.baselines, args.climate_witnesses, args.output)
