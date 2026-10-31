"""Recalculate tables and figures, with optional full source-coefficient replay."""

from pathlib import Path
import argparse
import json
import subprocess
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main(output, full=False, database=None, figures=True, quick=True):
    output = Path(output).resolve()
    if output in [ROOT, ROOT / "data", ROOT / "certificates", ROOT / "inputs"]:
        raise ValueError(
            "Select a separate output directory to preserve the scientific inputs."
        )
    output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()

    def run(script, *args):
        subprocess.run(
            [sys.executable, "-B", str(ROOT / script), *map(str, args)],
            check=True,
            cwd=ROOT,
        )

    if quick:
        run("code/check.py", "--output", output / "checks")
    run("code/consolidate.py", "--output", output)
    if full:
        args = ["--output", output / "matrix_replay"]
        if database is not None:
            args += ["--database", database]
        run("code/verify_matrices.py", *args)
        extra = ["--repository", ROOT, "--output", output / "safeguard",
                 "--intervals", output / "matrix_replay/exact_climate_only_intervals.json"]
        if database is not None:
            extra += ["--database", database]
        run("code/safeguard_sensitivity.py", *extra)
        check = ["--repository", ROOT, "--results", output / "safeguard",
                 "--output", output / "safeguard/independent_verification.json"]
        if database is not None:
            check += ["--database", database]
        run("code/verify_safeguard.py", *check)
    else:
        run("tests/test_safeguard.py", "--output", output / "safeguard_quick.json")
    run("code/tables.py", "--data", output / "data", "--output", output / "tables")
    if figures:
        run(
            "code/figures.py", "--data", output / "data", "--output", output / "figures"
        )
    run("tests/test_regression.py", "--output", output, *(["--full"] if full else []))
    report = {
        "source_coefficient_replay": full,
        "figure_generation": figures,
        "seconds": time.perf_counter() - start,
        "scope": (
            "Country calculations and certificate replay from original coefficients"
            if full
            else "Reader calculations from released numerical certificates"
        ),
    }
    (output / "reproduction.json").write_text(json.dumps(report, indent=2) + "\n")
    archive = output.parent / (output.name + ".zip")
    with zipfile.ZipFile(
        archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6
    ) as target:
        for path in sorted(output.rglob("*")):
            if path.is_file():
                target.write(path, path.relative_to(output))
    print(f'Reproduction completed in {report["seconds"]:.1f} seconds.')
    print(f"Output archive: {archive.name}")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs")
    parser.add_argument(
        "--full",
        action="store_true",
        help="Replay all source coefficients and certificates",
    )
    parser.add_argument(
        "--database",
        type=Path,
        help="Optional exact uncompressed European coefficient database",
    )
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument(
        "--skip-quick",
        action="store_true",
        help="Use after code/check.py has completed",
    )
    args = parser.parse_args()
    main(
        args.output, args.full, args.database, not args.no_figures, not args.skip_quick
    )
