"""Exact finite examples for the model's interpretive propositions.
These are mathematical unit tests, not additional European observations.
"""

from fractions import Fraction as F
from itertools import product
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def choices(pay):
    best = max(pay)
    return tuple(i for i, x in enumerate(pay) if x == best)


def main(output=None):
    checks = []

    def check(name, cond):
        if not cond:
            raise AssertionError(name)
        checks.append(name)

    # A single unit: maintaining the reference has lower return than degrading it.
    E = [F(4), F(6), F(3)]
    B = [F(5), F(2), F(7)]
    C = [F(8), F(4), F(9)]
    a = [max(c - C[0], F(0)) for c in C]
    h = [max(b - B[0], F(0)) for b in B]
    check(
        "preservation_and_degradation_have_same_zero_credits",
        a[0] == a[1] == h[0] == h[1] == 0,
    )
    for l, b in product(map(F, [0, 1, 3, 10]), repeat=2):
        pg = [E[s] + l * a[s] + b * h[s] for s in range(3)]
        check("preservation_strictly_dominated_for_gains_" + str((l, b)), pg[1] > pg[0])
        level = [E[s] + l * a[s] + b * B[s] for s in range(3)]
        signed = [E[s] + l * a[s] + b * (B[s] - B[0]) for s in range(3)]
        check(
            "recentring_preserves_entire_argmax_" + str((l, b)),
            choices(level) == choices(signed),
        )
        check(
            "recentring_constant_entitlement_" + str((l, b)),
            all(level[s] - signed[s] == b * B[0] for s in range(3)),
        )
        check(
            "nonnegative_menu_gives_weak_participation_" + str((l, b)),
            max(level) >= max(E),
        )
        k = F(17, 3)
        scaled = [E[s] + l * a[s] + (b / k) * (k * B[s]) for s in range(3)]
        check("national_B_scale_invariance_" + str((l, b)), scaled == level)
    # Two actions with preservation optimum under level payments, below outside option after recentering.
    e = [F(4), F(6)]
    b = [F(5), F(2)]
    beta = F(1)
    level = [e[s] + beta * b[s] for s in range(2)]
    signed = [e[s] + beta * (b[s] - b[0]) for s in range(2)]
    check(
        "same_choice_not_same_participation",
        choices(level) == choices(signed) == (0,)
        and max(signed) < max(e) <= max(level),
    )
    # A first-argmax witness can satisfy a goal even though another tied response fails.
    check(
        "weak_is_not_universal_implementation",
        choices([F(5), F(5)]) == (0, 1) and min([F(2), F(0)]) < 1 <= max([F(2), F(0)]),
    )
    # Independent additive outcomes: each coordinate minimum equals the true exhaustive minimum.
    gs = [[(F(4), F(0)), (F(0), F(4))], [(F(1), F(2)), (F(3), F(0))]]
    alloc = list(product(range(2), repeat=2))
    for k in (0, 1):
        exhaustive = min(sum(gs[i][x[i]][k] for i in range(2)) for x in alloc)
        check(
            "universal_additive_coordinate_minimum_" + str(k),
            exhaustive == sum(min(z[k] for z in g) for g in gs),
        )
    # A clipped credit is not an action-independent recentering.
    check(
        "clipping_changes_rankings",
        choices([e[s] + beta * max(b[s] - b[0], 0) for s in range(2)])
        != choices(level),
    )
    # Numerical count bounds are inequalities, not point identification.
    for a, national, common in [(0, 2, 12), (F(1, 40), 3, 14), (F(1, 20), 4, 14)]:
        check("exact_choice_count_gap_" + str(a), common - national >= 10)
        check("nested_band_bound_" + str(a), common - 4 in (8, 10))
    report = {
        "status": "PASS",
        "synthetic_checks": len(checks),
        "scope": "Exact illustrative propositions; not observations or new EU matrix simulations",
        "checks": checks,
    }
    if output is not None:
        Path(output).mkdir(parents=True, exist_ok=True)
        (Path(output) / "theory_checks.json").write_text(json.dumps(report, indent=2))
    print(
        json.dumps(
            {k: report[k] for k in ["status", "synthetic_checks", "scope"]}, indent=2
        )
    )

    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    main(args.output)
