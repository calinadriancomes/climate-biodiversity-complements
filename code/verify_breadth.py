"""Recompute exact climate-only coverage from all original actions without a solver.

Affine payoff envelopes use exact rational representations of the stored float64
coefficients. Coordinate sums and target comparisons also use exact arithmetic.
Only open price cells are needed to maximize robust coverage: at a breakpoint
the full best-response set contains the neighboring cell responses, so universal
compliance cannot improve there. Persistent affine ties are retained in every cell.
"""

from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import argparse, csv, hashlib, json, sqlite3, time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DB_SHA = "f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d"


def rational(x):
    return Fraction.from_float(float(x))


def upper_envelope(e, a):
    """Return (left breakpoint, persistent maximizers) for cells right of zero."""
    by_slope = {}
    for s, (ev, av) in enumerate(zip(e, a)):
        ev = float(ev)
        av = float(av)
        if av not in by_slope or ev > by_slope[av][0]:
            by_slope[av] = (ev, [s])
        elif ev == by_slope[av][0]:
            by_slope[av][1].append(s)
    hull = []
    starts = []
    for av, (ev, actions) in sorted(by_slope.items()):
        line = (rational(av), rational(ev), tuple(actions))
        while hull:
            cross = (hull[-1][1] - line[1]) / (line[0] - hull[-1][0])
            if len(hull) == 1 or cross > starts[-1]:
                break
            hull.pop()
            starts.pop()
        starts.append(None if not hull else cross)
        hull.append(line)
    first = 0
    for k in range(1, len(hull)):
        if starts[k] <= 0:
            first = k
        else:
            break
    out = [(Fraction(0), hull[first][2])]
    out.extend((starts[k], hull[k][2]) for k in range(first + 1, len(hull)))
    return out


def outcome_minima(g, actions):
    return tuple(rational(min(float(row[s]) for s in actions)) for row in g)


def append_interval(intervals, left, right):
    """Merge compliant open cells; interior breakpoints need not be compliant."""
    if intervals and intervals[-1][1] == left:
        intervals[-1] = (intervals[-1][0], right)
    else:
        intervals.append((left, right))


def country_intervals(e, b, c, targets):
    """Exact sums of outcome minima across independent full response sets."""
    a = np.maximum(c - c[:, 0, None], 0)
    events = defaultdict(lambda: [Fraction(0), Fraction(0), Fraction(0)])
    total = [Fraction(0), Fraction(0), Fraction(0)]
    segment_count = 0
    for i in range(len(e)):
        hull = upper_envelope(e[i], a[i])
        segment_count += len(hull)
        g = (e[i], b[i], c[i])
        prev = outcome_minima(g, hull[0][1])
        for k in range(3):
            total[k] += prev[k]
        for price, actions in hull[1:]:
            new = outcome_minima(g, actions)
            delta = events[price]
            for k in range(3):
                delta[k] += new[k] - prev[k]
            prev = new
    thresholds = [[rational(v) for v in t] for t in targets]
    intervals = [[] for _ in targets]
    left = Fraction(0)
    for right, delta in sorted(events.items()):
        for k, t in enumerate(thresholds):
            if all(v >= bound for v, bound in zip(total, t)):
                append_interval(intervals[k], left, right)
        for k in range(3):
            total[k] += delta[k]
        left = right
    for k, t in enumerate(thresholds):
        if all(v >= bound for v, bound in zip(total, t)):
            append_interval(intervals[k], left, None)
    return intervals, segment_count


def maximum_overlap(intervals_by_country):
    """Maximum coverage of positive-length price cells, excluding endpoints."""
    events = defaultdict(lambda: [set(), set()])
    for country, intervals in intervals_by_country.items():
        for left, right in intervals:
            events[left][0].add(country)
            if right is not None:
                events[right][1].add(country)
    active = set()
    best = set()
    best_interval = None
    points = sorted(events)
    for i, left in enumerate(points):
        starts, ends = events[left]
        active.difference_update(ends)
        active.update(starts)
        right = points[i + 1] if i + 1 < len(points) else None
        if len(active) > len(best):
            best = active.copy()
            best_interval = (left, right)
    return len(best), sorted(best), best_interval


def encode_endpoint(x):
    if x is None:
        return None
    return {
        "decimal": float(x),
        "numerator": str(x.numerator),
        "denominator": str(x.denominator),
    }


def exact_source_total(values):
    return float(sum((rational(v) for v in np.asarray(values).ravel()), Fraction(0)))


def check_close(actual, expected, label):
    if abs(actual - expected) > max(0.02, 5e-11 * max(1, abs(actual), abs(expected))):
        raise AssertionError((label, actual, expected))


def run(database, output):
    started = time.time()
    database = Path(database).resolve()
    output = Path(output)
    digest = hashlib.sha256(database.read_bytes()).hexdigest()
    if digest != DB_SHA:
        raise ValueError("Exact European coefficient checksum does not match")
    output.mkdir(parents=True, exist_ok=True)
    with (ROOT / "certificates/Frozen_Baselines.csv").open(newline="") as f:
        baselines = list(csv.DictReader(f))
    all_intervals = {0.0: {}, 0.025: {}, 0.05: {}}
    cee = set()
    source_checks = []
    segments = 0
    db = sqlite3.connect(database.as_uri() + "?mode=ro", uri=True)
    scenarios = db.execute("select * from scenario order by scenario_id").fetchall()
    if len(scenarios) != 13 or scenarios[0] != (0, "sustainable_current"):
        raise AssertionError("Source scenario order")
    for row in baselines:
        country = row["country"]
        x = np.array(
            db.execute(
                "SELECT sdu,scenario_id,net_econ_value,biodiversity,net_ghg_co2e FROM value WHERE country=? ORDER BY sdu,scenario_id",
                (country,),
            ).fetchall()
        )
        if x.shape != (int(row["n_sdu"]) * 13, 5):
            raise AssertionError(country)
        x = x.reshape(-1, 13, 5)
        if not np.array_equal(x[:, :, 1], np.tile(np.arange(13), (len(x), 1))):
            raise AssertionError(country)
        e, b, c = x[:, :, 2], x[:, :, 3], x[:, :, 4]
        if not (np.isfinite(e).all() and np.isfinite(b).all() and np.isfinite(c).all()):
            raise AssertionError("Nonfinite source value: " + country)
        if float(np.min(b)) < 0:
            raise AssertionError(
                "Negative biodiversity level invalidates nonnegative-transfer premise"
            )
        observed = {
            "E0": exact_source_total(e[:, 0]),
            "B0": exact_source_total(b[:, 0]),
            "C0": exact_source_total(c[:, 0]),
            "Emax": exact_source_total(np.max(e, axis=1)),
            "Bmax": exact_source_total(np.max(b, axis=1)),
        }
        for key, value in observed.items():
            check_close(value, float(row[key]), country + ":" + key)
        targets = []
        for alpha in all_intervals:
            target = np.array(
                [
                    float(row["E0"]),
                    float(row["B0"]) + alpha * (float(row["Bmax"]) - float(row["B0"])),
                    float(row["C0"]),
                ]
            )
            targets.append(target - 1e-8 * np.maximum(1, abs(target)))
        intervals, count = country_intervals(e, b, c, targets)
        segments += count
        for alpha, ii in zip(all_intervals, intervals):
            all_intervals[alpha][country] = ii
        if row["cee11"].lower() in ("true", "1"):
            cee.add(country)
        source_checks.append(
            {
                "country": country,
                "units": len(e),
                "minimum_B": float(np.min(b)),
                **observed,
                "exact_scalar_segments": count,
            }
        )
        print("Verified climate-only coverage:", country, flush=True)
    db.close()
    if (
        len(source_checks) != 23
        or sum(r["units"] for r in source_checks) != 87644
        or len(cee) != 11
    ):
        raise AssertionError("Source coverage")
    rows = []
    detailed = {}
    expected = {0.0: (21, 11, 6), 0.025: (20, 9, 5), 0.05: (19, 9, 5)}
    for alpha, intervals in all_intervals.items():
        individual = sum(bool(ii) for ii in intervals.values())
        maximum, countries, witness = maximum_overlap(intervals)
        cee_max, cee_countries, cee_witness = maximum_overlap(
            {c: v for c, v in intervals.items() if c in cee}
        )
        if (individual, maximum, cee_max) != expected[alpha]:
            raise AssertionError(
                (
                    "Breadth disagreement",
                    alpha,
                    (individual, maximum, cee_max),
                    expected[alpha],
                )
            )
        rows.append(
            {
                "alpha": alpha,
                "individual_climate_only_countries": individual,
                "maximum_common_climate_only_countries": maximum,
                "national_complement_lower_bound": 23 - individual,
                "common_complement_lower_bound": 23 - maximum,
                "CEE11_common_climate_only_maximum": cee_max,
                "national_complement_minimum_with_existing_remedies": 23 - individual,
                "minimum_additional_jurisdictions": individual - maximum,
                "exact_cell_witness_countries": countries,
                "exact_cell_left": encode_endpoint(witness[0]),
                "exact_cell_right": encode_endpoint(witness[1]),
            }
        )
        detailed[str(alpha)] = {
            country: [[encode_endpoint(a), encode_endpoint(b)] for a, b in ii]
            for country, ii in intervals.items()
        }
    summary = {
        "status": "PASS",
        "arithmetic": "Exact rational scalar envelopes and outcome sums from stored float64 technology; original 1e-8 effective-target allowance",
        "source_sha256": digest,
        "source_unchanged": hashlib.sha256(database.read_bytes()).hexdigest() == digest,
        "countries": 23,
        "units": 87644,
        "actions_per_unit": 13,
        "exact_scalar_segments": segments,
        "optimizer_calls": 0,
        "seconds": time.time() - started,
        "baseline_and_maximum_definition_checks": source_checks,
        "breadth": rows,
        "endpoint_logic": "Every breakpoint response set contains neighboring open-cell responses. Universal compliance cannot increase at a breakpoint, including zero. Maximum robust coverage is therefore attained on an open cell.",
    }
    (output / "exact_breadth_verification.json").write_text(
        json.dumps(summary, indent=2)
    )
    interval_export = {
        "scope": "Merged regions of compliant open price cells",
        "breakpoint_convention": "Internal breakpoints are not certified by this export. Additional tied actions can invalidate compliance at a breakpoint even when both adjacent open cells comply. Maximum common coverage over open cells remains exact.",
        "regions_by_ambition": detailed,
    }
    (output / "exact_climate_only_intervals.json").write_text(
        json.dumps(interval_export, indent=2)
    )
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--database", type=Path, required=True)
    ap.add_argument("--output", type=Path, default=ROOT / "outputs/matrix_replay")
    args = ap.parse_args()
    run(args.database, args.output)
