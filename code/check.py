"""Check scientific input identities and finite mathematical examples."""

from pathlib import Path
import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def check_inputs():
    entries = (ROOT / "inputs/checksums.sha256").read_text().splitlines()
    for entry in entries:
        digest, name = entry.split("  ", 1)
        path = ROOT / name
        with path.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        if actual != digest:
            raise ValueError(f"Scientific input checksum differs: {name}")
    with (ROOT / "certificates/Frozen_Baselines.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 23 or sum(int(r["n_sdu"]) for r in rows) != 87644:
        raise ValueError("The declared country or spatial-unit coverage differs.")
    return len(entries)


def main(output):
    started = time.perf_counter()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    count = check_inputs()
    for script in ["test_theory.py", "test_breadth.py"]:
        subprocess.run(
            [
                sys.executable,
                "-B",
                str(ROOT / "tests" / script),
                "--output",
                str(output),
            ],
            check=True,
        )
    subprocess.run(
        [sys.executable, "-B", str(ROOT / "tests/test_price_groups.py")], check=True
    )
    subprocess.run(
        [sys.executable, "-B", str(ROOT / "tests/test_fixed_groups.py")], check=True
    )
    report = {
        "scientific_input_checksums": count,
        "countries": 23,
        "spatial_units": 87644,
        "scope": "Input identities and synthetic mathematical examples",
        "coefficient_replay_performed": False,
        "seconds": time.perf_counter() - started,
    }
    (output / "quick_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Checked {count} scientific inputs and the finite mathematical examples.")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/checks")
    main(parser.parse_args().output)
