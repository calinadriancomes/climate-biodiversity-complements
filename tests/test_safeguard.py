"""Check the released Belgian safeguard frontier and complete policy records."""
from pathlib import Path
import argparse, csv, json, hashlib

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = [20, 14, 11, 9, 7, 5, 4]
RELAXED = [19, 13, 10, 8, 6, 5, 4]
SOURCE_SHA = "f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d"

def main(output=None):
    rows = list(csv.DictReader((ROOT / "data/safeguard_frontier.csv").open(newline="")))
    policies = list(csv.DictReader((ROOT / "certificates/safeguard_frontier_policies.csv").open(newline="")))
    checks = list(csv.DictReader((ROOT / "certificates/safeguard_frontier_checks.csv").open(newline="")))
    summary = json.loads((ROOT / "data/safeguard_analysis.json").read_text())
    assert summary["source_sha256"] == SOURCE_SHA
    assert (len(rows), len(policies), len(checks)) == (35, 161, 483)
    for i, delta in enumerate([0.0, 0.005, 0.01, 0.03, 0.05]):
        block = rows[i*7:(i+1)*7]
        assert all(float(r["Belgian_floor_relaxation"]) == delta for r in block)
        assert [int(r["q"]) for r in block] == list(range(1,8))
        lo, hi = [int(r["lower_bound"]) for r in block], [int(r["upper_bound"]) for r in block]
        assert lo == (PRIMARY if i <= 1 else RELAXED)
        assert hi == (RELAXED if i >= 3 else PRIMARY)
        assert all((r["exact"] == "True") == (a==b) for r,a,b in zip(block,lo,hi))
        assert all(r["all_three_response_bands"] == "True" for r in block)
    assert all(r["compliant"] == "True" for r in checks)
    for q in range(1, 8):
        block = [p for p in policies if int(p["q"]) == q]
        assert len(block) == 23 and len({p["country"] for p in block}) == 23
        assert len({p["lambda_climate"] for p in block}) <= q
        assert sum(float(p["beta_B"]) > 0 for p in block) == RELAXED[q-1]
        assert len([r for r in checks if int(r["q"]) == q]) == 69
    assert min(range(1,8),key=lambda q:q+RELAXED[q-1]) == 5
    assert [q for q in range(1,8) if q+RELAXED[q-1] == 11] == [5,6,7]
    result={"release_consistency":True,"scenarios":5,"country_policy_records":161,"response_checks":483,
            "source_sha256":SOURCE_SHA,"scope":"Published-record checks; full coefficient verification is separate"}
    if output:
        path=Path(output);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    main(a.output)
