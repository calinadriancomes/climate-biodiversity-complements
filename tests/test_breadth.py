"""Independent exhaustive-cell checks for the exact all-price coverage verifier."""

from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse, json, random, sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from verify_breadth import country_intervals, maximum_overlap, rational


def brute_cells(e, b, c, targets):
    a = np.maximum(c - c[:, 0, None], 0)
    breaks = {F(0)}
    for er, ar in zip(e, a):
        for s, t in combinations(range(len(er)), 2):
            if ar[s] != ar[t]:
                x = (rational(er[t]) - rational(er[s])) / (
                    rational(ar[s]) - rational(ar[t])
                )
                if x >= 0:
                    breaks.add(x)
    points = sorted(breaks)
    samples = []
    for i, left in enumerate(points):
        right = points[i + 1] if i + 1 < len(points) else None
        samples.append(
            (left, (left + right) / 2 if right is not None else left + 1, right)
        )

    def compliance(price, target):
        totals = [F(0), F(0), F(0)]
        for i in range(len(e)):
            u = [rational(ev) + price * rational(av) for ev, av in zip(e[i], a[i])]
            mx = max(u)
            ss = [j for j, v in enumerate(u) if v == mx]
            for k, g in enumerate([e, b, c]):
                totals[k] += min(rational(g[i, s]) for s in ss)
        return all(v >= rational(t) for v, t in zip(totals, target))

    return samples, compliance


def main(output=None):
    rng = random.Random(601)
    checks = 0
    # Repeated slopes, identical affine utilities with different ecological values,
    # negative intercepts and multiway endpoint ties occur in this finite ensemble.
    for trial in range(60):
        n = rng.randint(1, 4)
        m = rng.randint(2, 6)
        e = np.array([[rng.randint(-5, 8) / 2 for _ in range(m)] for _ in range(n)])
        b = np.array([[rng.randint(0, 10) / 2 for _ in range(m)] for _ in range(n)])
        c = np.array([[rng.randint(-4, 12) / 2 for _ in range(m)] for _ in range(n)])
        targets = [
            np.array(
                [
                    rng.randint(-3, 6) * n / 2,
                    rng.randint(0, 8) * n / 2,
                    rng.randint(-2, 8) * n / 2,
                ]
            )
            for _ in range(3)
        ]
        intervals, _ = country_intervals(e, b, c, targets)
        samples, compliance = brute_cells(e, b, c, targets)
        for k, t in enumerate(targets):
            for i, (left, price, right) in enumerate(samples):
                expected = compliance(price, t)
                got = any(
                    lo < price and (hi is None or price < hi) for lo, hi in intervals[k]
                )
                assert got == expected, (trial, k, price, intervals[k])
                checks += 1
                if compliance(left, t):
                    assert expected
                    if i:
                        assert compliance(samples[i - 1][1], t)
                    checks += 1
    # Closing incompatible endpoints would falsely increase common coverage.
    ii = {"one": [(F(0), F(1))], "two": [(F(1), F(2))]}
    assert maximum_overlap(ii)[0] == 1
    checks += 1
    ii["three"] = [(F(1, 2), F(3, 2))]
    assert maximum_overlap(ii)[0] == 2
    checks += 1
    # A third action can tie only at an internal breakpoint and violate a
    # safeguard although both neighboring cells comply. The interval export
    # merges cells for coverage counting; it does not certify such breakpoints.
    e = np.array([[0.0, -1.0, -0.5]])
    b = np.array([[1.0, 1.0, 0.0]])
    c = np.array([[0.0, 1.0, 0.5]])
    targets = [np.array([-1.0, 0.5, 0.0])]
    intervals, _ = country_intervals(e, b, c, targets)
    _, compliance = brute_cells(e, b, c, targets)
    assert compliance(F(1, 2), targets[0])
    checks += 1
    assert compliance(F(3, 2), targets[0])
    checks += 1
    assert not compliance(F(1), targets[0])
    checks += 1
    assert maximum_overlap({"one": intervals[0]})[0] == 1
    checks += 1
    # Nonnegative transfers preserve outside-option choices whenever enrollment
    # is indifferent. Check all cases in a small complete integer grid.
    for e0, e1, t0, t1 in product(range(-2, 3), range(-2, 3), range(3), range(3)):
        outside = max(e0, e1)
        u = [e0 + t0, e1 + t1]
        top = max(u)
        assert top >= outside
        if top == outside:
            assert all(u[s] == top for s, ev in enumerate([e0, e1]) if ev == outside)
        checks += 1
    report = {
        "status": "PASS",
        "checks": checks,
        "scope": "Independent exhaustive comparisons on synthetic finite technologies and exact costless opt-out property; no additional empirical observations",
    }
    if output is not None:
        Path(output).mkdir(parents=True, exist_ok=True)
        (Path(output) / "exact_breadth_tests.json").write_text(
            json.dumps(report, indent=2)
        )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    main(args.output)
