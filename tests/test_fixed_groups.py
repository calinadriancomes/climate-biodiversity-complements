"""Independent finite examples and released fixed-group record checks."""
from pathlib import Path
from fractions import Fraction as F
import csv
import json
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from fixed_groups import group_coverage, regional_certificate


def encode(value):
    if value is None:
        return None
    value = F(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


class FixedGroupTests(unittest.TestCase):
    def test_required_low_rate_limits_its_entire_group(self):
        regions = {"anchor": [], "low": [[encode(0), encode(1)]],
                   "high1": [[encode(3), encode(5)]],
                   "high2": [[encode(3), None]]}
        result = group_coverage(regions, list(regions), F(2))
        self.assertEqual(result["maximum_climate_only_coverage"], 1)
        self.assertEqual(result["coverage_witness_countries"], ["low"])

    def test_strict_dominance_covers_zero_biodiversity_and_unbounded_rates(self):
        e, b, c = (np.array([[10., 9.]]), np.array([[1., 2.]]), np.array([[0., 10.]]))
        certificate = regional_certificate("test", e, b, c, {"E0": 10}, F(1, 5))
        self.assertEqual(certificate["land_value_upper_exact"], "9")
        self.assertEqual(certificate["removed_actions"], 1)
        self.assertGreater(certificate["certified_shortfall"], .999)
        for climate in [F(1, 5), F(1), F(1000000)]:
            for biodiversity in [F(0), F(1), F(1000000)]:
                self.assertGreater(F(9) + climate * 10 + biodiversity * 2,
                                   F(10) + biodiversity)

    def test_complete_fixed_group_record(self):
        folder = ROOT / "data/fixed_groups"
        report = json.loads((folder / "fixed_group_verification.json").read_text())
        with (folder / "fixed_group_policies.csv").open(newline="") as stream:
            policies = list(csv.DictReader(stream))
        with (folder / "fixed_group_response_checks.csv").open(newline="") as stream:
            checks = list(csv.DictReader(stream))
        self.assertEqual(len(policies), 23)
        self.assertEqual(len({r["country"] for r in policies}), 23)
        self.assertEqual(len(checks), 69)
        self.assertTrue(all(r["compliant"] == "True" for r in checks))
        self.assertEqual(sum(float(r["beta_B"]) > 0 for r in policies), 19)
        self.assertEqual(sum(r["maximum_climate_only_coverage"] for r in report["coverage_bounds"]), 4)
        self.assertEqual(report["complement_lower_bound"], report["complement_upper_bound"])
        for group, countries in report["fixed_groups"].items():
            block = [r for r in policies if r["group"] == group]
            self.assertEqual({r["country"] for r in block}, set(countries))
            self.assertEqual(len({r["lambda_climate"] for r in block}), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
