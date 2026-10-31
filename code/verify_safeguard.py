"""Independent checks of the Belgian safeguard sensitivity.

Uses pairwise support inequalities instead of an upper-envelope constructor,
exact rational dual evaluation without an optimizer and long-double responses.
"""

from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse
import csv
import hashlib
import json
import sqlite3
import gzip
import tempfile
import numpy as np


def main(args):
    report = json.loads((args.results / 'safeguard_analysis.json').read_text())
    assert hashlib.sha256(args.database.read_bytes()).hexdigest() == report['source_sha256']
    base = {r['country']: r for r in csv.DictReader((args.repository / 'certificates/Frozen_Baselines.csv').open())}
    with sqlite3.connect(args.database.resolve().as_uri() + '?mode=ro', uri=True) as db:
        source = db.execute('select sdu,scenario_id,net_econ_value,biodiversity,net_ghg_co2e from value where country="Belgium" order by sdu,scenario_id').fetchall()
        arrays = {}
        for country in base:
            values = np.array(db.execute('select net_econ_value,biodiversity,net_ghg_co2e from value where country=? order by sdu,scenario_id', (country,)).fetchall()).reshape(-1, 13, 3)
            a = np.maximum(values[:, :, 2] - values[:, [0], 2], 0.0)
            arrays[country] = [values[:, :, k].astype(np.longdouble) for k in range(3)] + [a.astype(np.longdouble)]
    units = defaultdict(list)
    for unit, action, e, b, c in source:
        units[unit].append((action, e, b, c))
    values, paid, starts, ends = {}, {}, defaultdict(list), defaultdict(list)
    for unit, rows in units.items():
        assert [r[0] for r in rows] == list(range(13))
        values[unit] = [(F(e), F(b), F(c)) for s, e, b, c in rows]
        paid[unit] = [F(max(c - rows[0][3], 0.0)) for s, e, b, c in rows]
        for s in range(13):
            lo, hi = F(0), None
            for t in range(13):
                de = values[unit][s][0] - values[unit][t][0]
                da = paid[unit][s] - paid[unit][t]
                if da == 0 and de < 0:
                    hi = F(-1)
                    break
                if da > 0:
                    lo = max(lo, -de / da)
                elif da < 0:
                    hi = -de / da if hi is None else min(hi, -de / da)
            if hi is None or lo <= hi:
                starts[lo].append((unit, s))
                if hi is not None:
                    ends[hi].append((unit, s))
    b0 = base['Belgium']
    raw = [F(float(b0['E0'])), F(float(b0['B0']) + .05 * (float(b0['Bmax']) - float(b0['B0']))), F(float(b0['C0']))]
    effective = [v - F(1, 100000000) * max(F(1), abs(v)) for v in raw]
    active = {unit: set() for unit in units}
    old = {unit: [F(0)] * 3 for unit in units}
    total = [F(0)] * 3

    def update(unit):
        new = [min(values[unit][s][k] for s in active[unit]) for k in range(3)]
        for k in range(3):
            total[k] += new[k] - old[unit][k]
        old[unit] = new

    best = None
    point_count, cell_count = 0, 0
    for point in sorted(set(starts) | set(ends)):
        changed = set()
        for unit, s in starts[point]:
            active[unit].add(s)
            changed.add(unit)
        for unit in changed:
            update(unit)
        assert all(active.values())
        point_count += 1
        if total[1] >= effective[1] and total[2] >= effective[2]:
            best = total[0] if best is None else max(best, total[0])
        changed = set()
        for unit, s in ends[point]:
            active[unit].remove(s)
            changed.add(unit)
        for unit in changed:
            update(unit)
        assert all(active.values())
        cell_count += 1
        if total[1] >= effective[1] and total[2] >= effective[2]:
            best = total[0] if best is None else max(best, total[0])
    assert best == F(report['Belgian_climate_only_ceiling']['exact_maximum_E'])
    dual_checks = []
    for entry in report['regional_bounds']:
        cut, yb, yc = F(entry['cutoff_exact']), F(entry['dual_B_exact']), F(entry['dual_C20_exact'])
        assert min(yb, yc) >= 0
        upper = F(0)
        count = 0
        for unit in units:
            retained = []
            for s, (es, bs, cs) in enumerate(values[unit]):
                dominated = any(at >= paid[unit][s] and bt >= bs and et - es + cut * (at - paid[unit][s]) > 0
                                for (et, bt, ct), at in zip(values[unit], paid[unit]))
                if not dominated:
                    retained.append(es + yb * bs + yc * cs)
            assert retained
            count += len(retained)
            upper += max(retained)
        upper -= yb * effective[1] + yc * effective[2]
        assert upper == F(entry['upper_E_exact']) and count == entry['retained_actions']
        dual_checks.append({'threshold_country': entry['threshold_country'], 'exact_upper_verified': True, 'upper_E': float(upper)})
    rows = list(csv.DictReader((args.results / 'safeguard_frontier_policies.csv').open()))
    passed, min_margin = 0, float('inf')
    ld = np.longdouble
    for row in rows:
        country = row['country']
        e, b, c, a = arrays[country]
        lam, beta = ld(float(row['lambda_climate'])), ld(float(row['beta_B']))
        lc, bb = lam * a, beta * b
        u = e + lc + bb
        guard = 64 * np.finfo(ld).eps * np.maximum(ld(1), (abs(e) + abs(lc) + abs(bb)).max(1))
        scale = np.maximum(ld(1), abs(e).max(1))
        ref = base[country]
        target = np.array([float(ref['E0']) * float(row['country_floor_ratio']),
                           float(ref['B0']) + .05 * (float(ref['Bmax']) - float(ref['B0'])), float(ref['C0'])], dtype=ld)
        eff = target - ld(1e-8) * np.maximum(ld(1), abs(target))
        for eps in (0, 1e-6, 1e-4):
            mask = u.max(1)[:, None] - u <= ld(eps) * scale[:, None] + guard[:, None]
            low = np.array([np.where(mask, g, np.inf).min(1).sum(dtype=ld) for g in (e, b, c)])
            assert (low >= eff).all(), (country, row['q'], eps)
            passed += 1
            min_margin = min(min_margin, float(((low - target) / np.maximum(ld(1), abs(target))).min()))
    for q, entry in enumerate(report['relaxed_floor_exact_frontier'], 1):
        block = [row for row in rows if int(row['q']) == q]
        assert len(block) == 23
        assert len({row['lambda_climate'] for row in block}) <= q
        assert sum(float(row['beta_B']) > 0 for row in block) == entry['lower_bound']
    result = {'status': 'verified', 'Belgian_climate_only_method': 'All pairwise closed support intervals; every boundary and every open cell including zero and the unbounded tail',
              'climate_only_boundaries': point_count, 'climate_only_open_cells': cell_count,
              'climate_only_exact_maximum_E': str(best), 'dual_bound_checks': dual_checks,
              'complete_policies': 7, 'long_double_response_checks': passed,
              'minimum_normalized_margin': min_margin, 'source_unchanged': hashlib.sha256(args.database.read_bytes()).hexdigest() == report['source_sha256']}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--database', type=Path, help='Optional exact uncompressed source database')
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='safeguard_independent_') as temporary:
        if args.database is None:
            args.database = Path(temporary) / 'eu23.sqlite'
            parts = [args.repository / 'inputs' / f'eu23_coefficients.sqlite.gz.part{i:02d}' for i in (1, 2)]
            args.database.write_bytes(gzip.decompress(b''.join(p.read_bytes() for p in parts)))
        main(args)
