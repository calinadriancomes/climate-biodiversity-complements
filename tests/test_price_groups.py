"""Independent finite examples for shared-price coverage and the regional certificate."""

from pathlib import Path
from fractions import Fraction as F
from itertools import product
import random
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from price_groups import constrained_coverage, anchor_certificate, response_checks
from belgium_economic_diagnostic import envelope


def encode(value):
    if value is None:
        return None
    value = F(value)
    return {
        "decimal": float(value),
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def encoded(regions):
    return {
        country: [[encode(a), encode(b)] for a, b in cells]
        for country, cells in regions.items()
    }


def exhaustive_coverage(regions, cutoff, q):
    points = sorted(
        {F(0), cutoff}
        | {
            v
            for cells in regions.values()
            for pair in cells
            for v in pair
            if v is not None
        }
    )
    samples = [(a + b) / 2 for a, b in zip(points, points[1:])] + [points[-1] + 1]
    best = 0
    for prices in product(samples, repeat=q):
        if not any(price < cutoff for price in prices):
            continue
        count = sum(
            any(
                a < price and (b is None or price < b)
                for price in prices
                for a, b in cells
            )
            for cells in regions.values()
        )
        best = max(best, count)
    return best


class SharedPriceTests(unittest.TestCase):
    def test_anchor_changes_complete_policy_bound(self):
        regions = {"a": [(F(0), F(1))], "b": [(F(3), F(4))], "c": [(F(3), F(5))]}
        got = constrained_coverage(encoded(regions), F(2))
        self.assertEqual(
            [row["maximum_climate_only_coverage_given_anchor"] for row in got], [1, 3]
        )

    def test_independent_exhaustive_price_tuples(self):
        rng = random.Random(2726)
        for _ in range(25):
            regions = {}
            for country in ["a", "b", "c", "d"]:
                cells = []
                for _ in range(rng.randint(0, 2)):
                    left = F(rng.randint(0, 4))
                    right = left + rng.randint(1, 3)
                    cells.append((left, None if rng.randrange(7) == 0 else right))
                regions[country] = cells
            got = constrained_coverage(encoded(regions), F(5, 2))
            for row in got[:3]:
                q = row["maximum_climate_coefficients"]
                self.assertEqual(
                    row["maximum_climate_only_coverage_given_anchor"],
                    exhaustive_coverage(regions, F(5, 2), q),
                )

    def test_rational_regional_dominance(self):
        e, b, c = (
            np.array([[10.0, 9.0]]),
            np.array([[1.0, 2.0]]),
            np.array([[0.0, 10.0]]),
        )
        certificate, witnesses = anchor_certificate(
            e, b, c, {"E0": 10}, np.array([1]), F(1, 5)
        )
        self.assertEqual(certificate["excluded_actions"], 1)
        self.assertEqual(certificate["retained_land_value_upper"], 9)
        self.assertGreater(certificate["certified_land_value_shortfall"], 0.999)
        self.assertEqual(witnesses[0]["action"], 0)
        self.assertEqual(witnesses[0]["dominating_action"], 1)
        for lam, beta in product([F(1, 5), F(1), F(100)], [F(0), F(1), F(100)]):
            self.assertGreater(F(9) + lam * 10 + beta * 2, F(10) + beta)

    def test_biodiversity_envelope_against_all_actions(self):
        rng = random.Random(1126)
        for trial in range(40):
            e = [F(rng.randint(-8, 8), 2) for _ in range(6)]
            b = [F(rng.randint(0, 8), 2) for _ in range(6)]
            cells = envelope(e, b)
            for index, (left, selected) in enumerate(cells):
                right = cells[index + 1][0] if index + 1 < len(cells) else None
                price = (left + right) / 2 if right is not None else left + 1
                utility = [ee + price * bb for ee, bb in zip(e, b)]
                expected = [
                    i for i, value in enumerate(utility) if value == max(utility)
                ]
                self.assertEqual(sorted(selected), expected, (trial, price))

    def test_wider_choice_band_can_fail(self):
        e, b, c = (
            np.array([[10.0, 10.0]]),
            np.array([[1.0, 2.0]]),
            np.array([[0.0, 0.0]]),
        )
        result = response_checks(
            e, b, c, {"E0": 10, "B0": 1, "Bmax": 2, "C0": 0}, 0.5, 0.0, 0.0005
        )
        self.assertEqual([r["compliant"] for r in result], [True, True, False])


if __name__ == "__main__":
    unittest.main(verbosity=2)
