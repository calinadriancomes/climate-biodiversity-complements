"""Compare regenerated scientific records with the independent numerical reference."""

from pathlib import Path
import argparse
import csv
import json
import math
from PIL import Image
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VARIABLE_FIELDS = {"seconds", "command", "source_certificate_files"}


def equivalent(expected, actual, where):
    if isinstance(expected, dict):
        for key, value in expected.items():
            if key not in VARIABLE_FIELDS:
                if key not in actual:
                    raise AssertionError(f"Missing field: {where}.{key}")
                equivalent(value, actual[key], where + "." + key)
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise AssertionError(f"Record count differs: {where}")
        for index, (left, right) in enumerate(zip(expected, actual)):
            equivalent(left, right, f"{where}[{index}]")
    elif isinstance(expected, bool) or expected is None:
        if expected is not actual:
            raise AssertionError(f"Logical value differs: {where}")
    elif isinstance(expected, int):
        if expected != actual:
            raise AssertionError(
                f"Integer value differs: {where}: {expected} != {actual}"
            )
    elif isinstance(expected, float):
        if not math.isclose(expected, actual, rel_tol=5e-12, abs_tol=2e-12):
            raise AssertionError(
                f"Numerical value differs: {where}: {expected} != {actual}"
            )
    else:
        if expected != actual:
            # CSV fields retain exact integers and compare decimal values using
            # a numerical tolerance substantially smaller than displayed rounding.
            try:
                a, b = float(expected), float(actual)
            except (ValueError, TypeError):
                raise AssertionError(
                    f"Text differs: {where}: {expected!r} != {actual!r}"
                )
            if not math.isclose(a, b, rel_tol=5e-12, abs_tol=2e-12):
                raise AssertionError(f"Numerical value differs: {where}: {a} != {b}")


def load(path):
    if path.suffix == ".json":
        return json.loads(path.read_text())
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main(output, full=False):
    output = Path(output)
    compared = []
    for expected in sorted((ROOT / "data").glob("*")):
        if expected.suffix not in (".csv", ".json"):
            continue
        if expected.name.startswith("safeguard_"):
            # Source-derived sensitivity is calculated separately in full mode.
            continue
        actual = output / "data" / expected.name
        if not actual.is_file():
            raise AssertionError("Missing regenerated table: " + expected.name)
        equivalent(load(expected), load(actual), "data/" + expected.name)
        compared.append("data/" + expected.name)
    if full:
        reference = json.loads((ROOT / "results/reference.json").read_text())[
            "matrix_replay"
        ]
        for name, expected in reference.items():
            actual = output / "matrix_replay" / name
            if not actual.is_file():
                raise AssertionError("Missing coefficient verification: " + name)
            equivalent(expected, load(actual), "matrix_replay/" + name)
            compared.append("matrix_replay/" + name)
    if full:
        sensitivity = output / "safeguard/safeguard_frontier.csv"
        assert sensitivity.exists(), "Missing original-coefficient safeguard frontier"
        equivalent(load(ROOT / "data/safeguard_frontier.csv"), load(sensitivity),
                   "safeguard/frontier")
        independent = json.loads((output / "safeguard/independent_verification.json").read_text())
        assert independent["status"] == "verified"
        assert independent["long_double_response_checks"] == 483
        assert independent["source_unchanged"]
        compared.append("safeguard/frontier")
    for directory in ["price_groups", "fixed_groups"]:
        for expected in sorted((ROOT / "data" / directory).glob("*")):
            actual = output / "data" / directory / expected.name
            if not actual.is_file():
                raise AssertionError("Missing policy record: " + expected.name)
            equivalent(load(expected), load(actual), "data/" + directory + "/" + expected.name)
            compared.append("data/" + directory + "/" + expected.name)
    parity = []
    for png in sorted((output / "figures").glob("*.png")):
        tiff = png.with_suffix(".tiff")
        with Image.open(png) as left, Image.open(tiff) as right:
            if not np.array_equal(
                np.asarray(left.convert("RGB")), np.asarray(right.convert("RGB"))
            ):
                raise AssertionError("PNG and TIFF pixels differ: " + png.stem)
        parity.append(png.stem)
    report = {
        "reference": "Numerical records from the preceding complete coefficient verification",
        "coefficient_replay_compared": full,
        "records_compared": compared,
        "relative_numeric_tolerance": 5e-12,
        "absolute_numeric_tolerance": 2e-12,
        "PNG_TIFF_pixel_parity": parity,
    }
    (output / "checks").mkdir(parents=True, exist_ok=True)
    (output / "checks/regression.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"Regression: {len(compared)} numerical record files agree; "
        f"{len(parity)} figure pairs have identical pixels."
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs")
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()
    main(args.output, args.full)
