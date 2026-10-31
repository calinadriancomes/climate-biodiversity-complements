"""Verify policies and certificates against the exact European coefficients."""

from pathlib import Path
import argparse, gzip, hashlib, os, shutil, subprocess, sys, tempfile, json, time

ROOT = Path(__file__).resolve().parents[1]
DB_SHA = "f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d"


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--database", type=Path, help="Optional uncompressed exact European SQLite"
    )
    ap.add_argument("--output", type=Path, default=ROOT / "outputs/matrix_replay")
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    start = time.time()
    with tempfile.TemporaryDirectory(prefix="eu23_coefficients_") as temp:
        db = a.database
        if db is None:
            db = Path(temp) / "eu23.sqlite"
            compressed = ROOT / "inputs/eu23_coefficients.sqlite.gz"
            if not compressed.is_file():
                chunks = [
                    ROOT / "inputs" / f"eu23_coefficients.sqlite.gz.part{i:02d}"
                    for i in (1, 2)
                ]
                if not all(x.is_file() for x in chunks):
                    raise FileNotFoundError(
                        "The two European coefficient data parts are required."
                    )
                compressed = Path(temp) / "eu23_coefficients.sqlite.gz"
                with compressed.open("wb") as output:
                    for part in chunks:
                        with part.open("rb") as source:
                            shutil.copyfileobj(source, output)
            with gzip.open(compressed, "rb") as inp, db.open("wb") as out:
                shutil.copyfileobj(inp, out)
        db = db.resolve()
        if sha(db) != DB_SHA:
            raise ValueError("Exact European coefficient checksum does not match")
        env = dict(
            os.environ,
            CBC_MATRIX_DATABASE=str(db),
            CBC_VERIFY_OUTPUT=str(a.output.resolve()),
            OMP_NUM_THREADS="2",
            OPENBLAS_NUM_THREADS="2",
        )
        subprocess.run(
            [sys.executable, "-B", str(ROOT / "code/verify_coefficients.py")],
            env=env,
            check=True,
        )
        for script in [
            "verify_breadth.py",
            "verify_additional.py",
            "verify_weak_exclusions.py",
        ]:
            subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(ROOT / "code" / script),
                    "--database",
                    str(db),
                    "--output",
                    str(a.output.resolve()),
                ],
                env=env,
                check=True,
            )
        subprocess.run(
            [
                sys.executable,
                "-B",
                str(ROOT / "code/price_groups.py"),
                "--database",
                str(db),
                "--intervals",
                str(a.output / "exact_climate_only_intervals.json"),
                "--output",
                str(a.output.parent / "price_groups"),
            ],
            env=env,
            check=True,
        )
        subprocess.run(
            [
                sys.executable,
                "-B",
                str(ROOT / "code/belgium_economic_diagnostic.py"),
                "--database",
                str(db),
                "--baselines",
                str(ROOT / "certificates/Frozen_Baselines.csv"),
                "--climate-witnesses",
                str(ROOT / "certificates/sensitivity/ClimateOnly_Alpha_Witnesses.csv"),
                "--output",
                str(a.output.parent / "price_groups"),
            ],
            env=env,
            check=True,
        )
        subprocess.run(
            [sys.executable, "-B", str(ROOT / "code/fixed_groups.py"),
             "--database", str(db),
             "--intervals", str(a.output / "exact_climate_only_intervals.json"),
             "--output", str(a.output.parent / "fixed_groups")],
            env=env, check=True,
        )
        fixed_extension = a.output.parent / "data/fixed_groups"
        fixed_extension.mkdir(parents=True, exist_ok=True)
        for path in sorted((a.output.parent / "fixed_groups").iterdir()):
            if path.is_file():
                shutil.copyfile(path, fixed_extension / path.name)
        reader_extension = a.output.parent / "data/price_groups"
        reader_extension.mkdir(parents=True, exist_ok=True)
        for name in [
            "price_group_frontier.csv",
            "price_group_policy_checks.csv",
            "price_group_verification.json",
            "climate_ceiling_certificate.json",
            "belgium_economic_ceiling.json",
        ]:
            shutil.copyfile(
                a.output.parent / "price_groups" / name, reader_extension / name
            )
        summary = json.loads((a.output / "verification.json").read_text())
        breadth = json.loads((a.output / "exact_breadth_verification.json").read_text())
        additional = json.loads((a.output / "additional_verification.json").read_text())
        weak = json.loads((a.output / "weak_climate_only_exclusions.json").read_text())
        primary = next(r for r in breadth["breadth"] if r["alpha"] == 0.05)
        summary["common_complement_lower_from_support_bounds"] = primary[
            "common_complement_lower_bound"
        ]
        summary["global_common_complement_upper_constructed"] = summary[
            "robust_constructions"
        ]["common_final"]["positive_B"]
        summary["national_minimum_complements_certified"] = primary[
            "national_complement_lower_bound"
        ]
        assert (
            summary["national_minimum_complements_certified"]
            == summary["robust_constructions"]["national_final"]["positive_B"]
        )
        summary["all_ambition_breadth_recomputed_from_source"] = breadth["breadth"]
        summary["coverage_bounds"] = (
            "Exact rational reconstruction from all original actions and cross-country price-cell aggregation"
        )
        summary["exclusion_sensitivity_recomputed"] = additional[
            "exclusion_sensitivity"
        ]
        summary["weak_climate_only_exclusions_recomputed"] = weak["countries"]
        (a.output / "verification.json").write_text(json.dumps(summary, indent=2))
        if sha(db) != DB_SHA:
            raise AssertionError("Input database changed")
    (a.output / "execution.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "coefficient_sha256": DB_SHA,
                "source_unchanged": True,
                "optimizers_called": 0,
                "seconds": time.time() - start,
                "command": "python code/verify_matrices.py",
            },
            indent=2,
        )
    )
    print(
        "Source verification and exact shared-price analysis completed; source unchanged."
    )


if __name__ == "__main__":
    main()
