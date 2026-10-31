"""Exact Belgian safeguard bounds and a complete EU23 policy sensitivity.

The original technology and ecological targets stay fixed. Only Belgium's
aggregate land-value floor is relaxed in the explicitly labeled scenarios.
This script does not edit the repository or its publication data.
"""

from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import gzip
import subprocess
import sys
import sqlite3
import tempfile
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

F = lambda x: Fraction.from_float(float(x))
SHA = 'f82ee0744c2eb30c6b72f5fbc993a42645ebb09264917f9b940198a02c36497d'


def write_csv(path, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def targets(base, alpha=0.05, floor_ratio=1.0):
    raw = [F(float(base['E0']) * floor_ratio),
           F(float(base['B0']) + alpha * (float(base['Bmax']) - float(base['B0']))),
           F(base['C0'])]
    return raw, [v - Fraction(1, 100000000) * max(Fraction(1), abs(v)) for v in raw]


def exact_cells(e, b, c, envelope, fixed_climate=None):
    """All open cells of the remaining coefficient, with persistent ties."""
    a = np.maximum(c - c[:, [0]], 0.0)
    events = defaultdict(lambda: [Fraction(0)] * 3)
    total = [Fraction(0)] * 3
    for i in range(len(e)):
        if fixed_climate is None:
            intercepts = [F(v) for v in e[i]]
            slopes = [F(v) for v in a[i]]
        else:
            intercepts = [F(e[i, s]) + F(fixed_climate) * F(a[i, s]) for s in range(13)]
            slopes = [F(v) for v in b[i]]
        pieces = envelope(intercepts, slopes)
        previous = [min(F(g[i, s]) for s in pieces[0][1]) for g in (e, b, c)]
        total = [x + y for x, y in zip(total, previous)]
        for coefficient, actions in pieces[1:]:
            current = [min(F(g[i, s]) for s in actions) for g in (e, b, c)]
            events[coefficient] = [x + y - z for x, y, z in zip(events[coefficient], current, previous)]
            previous = current
    left = Fraction(0)
    cells = []
    for right in sorted(events) + [None]:
        cells.append((left, right, total.copy()))
        if right is not None:
            total = [x + y for x, y in zip(total, events[right])]
            left = right
    return cells


def band_checks(e, b, c, lam, beta, raw, effective):
    a = np.maximum(c - c[:, [0]], 0.0)
    lc, bb = lam * a, beta * b
    payoff = e + lc + bb
    gap = payoff.max(1)[:, None] - payoff
    scale = np.maximum(1, abs(e).max(1))
    guard = 32 * np.finfo(float).eps * np.maximum(1, (abs(e) + abs(lc) + abs(bb)).max(1))
    result = []
    for eps in (0.0, 1e-6, 1e-4):
        admitted = gap <= eps * scale[:, None] + guard[:, None]
        low = [math.fsum(np.where(admitted, g, np.inf).min(1)) for g in (e, b, c)]
        result.append({'response_band': eps, 'compliant': all(v >= float(t) for v, t in zip(low, effective)),
                       'E_min': low[0], 'B_min': low[1], 'C20_min': low[2],
                       'minimum_normalized_margin': min((v - float(t)) / max(1, abs(float(t))) for v, t in zip(low, raw))})
    return result


def regional_upper(e, b, c, cutoff, effective):
    """Valid rational land-value upper bound throughout lambda >= cutoff.

    The optimizer proposes nonnegative dual weights. The published upper bound
    is then recomputed directly with exact rational arithmetic and is valid
    independently of solver success, its primal values or its tolerances.
    """
    a = np.maximum(c - c[:, [0]], 0.0)
    retained = np.ones(e.shape, dtype=bool)
    er = [[F(v) for v in row] for row in e]
    br = [[F(v) for v in row] for row in b]
    cr = [[F(v) for v in row] for row in c]
    ar = [[F(v) for v in row] for row in a]
    for i in range(len(e)):
        for s in range(13):
            for t in range(13):
                if ar[i][t] >= ar[i][s] and br[i][t] >= br[i][s]:
                    if er[i][t] - er[i][s] + cutoff * (ar[i][t] - ar[i][s]) > 0:
                        retained[i, s] = False
                        break
    assert retained.any(1).all()
    ii, ss = np.where(retained)
    eq = csr_matrix((np.ones(len(ii)), (ii, np.arange(len(ii)))), shape=(len(e), len(ii)))
    escale = max(1.0, abs(float(effective[0])))
    scales = np.maximum(1.0, np.abs([float(effective[1]), float(effective[2])]))
    sol = linprog(-e[ii, ss] / escale,
                  A_ub=-np.vstack([b[ii, ss], c[ii, ss]]) / scales[:, None],
                  b_ub=-np.array([float(effective[1]), float(effective[2])]) / scales,
                  A_eq=eq, b_eq=np.ones(len(e)), bounds=(0, None), method='highs')
    assert sol.success, sol.message
    weights = [F(max(0.0, -float(v) * escale / scale)) for v, scale in zip(sol.ineqlin.marginals, scales)]
    upper = sum((max(er[i][s] + weights[0] * br[i][s] + weights[1] * cr[i][s]
                     for s in np.flatnonzero(retained[i])) for i in range(len(e))), Fraction(0))
    upper -= weights[0] * effective[1] + weights[1] * effective[2]
    return {'climate_cutoff': float(cutoff), 'cutoff_exact': str(cutoff),
            'retained_actions': int(retained.sum()), 'upper_E': float(upper),
            'upper_E_exact': str(upper), 'dual_B_exact': str(weights[0]),
            'dual_C20_exact': str(weights[1]),
            'minimum_effective_floor_relaxation': float(effective[0] - upper),
            'scope': 'All finite climate coefficients at least cutoff and all nonnegative biodiversity coefficients; fractional retained-action relaxation'}, upper


def main(args):
    repo, db, out = args.repository.resolve(), args.database.resolve(), args.output.resolve()
    if out == repo or any(out == repo / name or (repo / name) in out.parents for name in ('inputs', 'data', 'certificates')):
        raise ValueError('Use a separate output directory.')
    out.mkdir(parents=True, exist_ok=True)
    assert hashlib.sha256(db.read_bytes()).hexdigest() == SHA
    spec = importlib.util.spec_from_file_location('fixed_coefficient', repo / 'code/belgium_economic_diagnostic.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    bases = {r['country']: r for r in csv.DictReader((repo / 'certificates/Frozen_Baselines.csv').open())}
    mats = {}
    with sqlite3.connect(db.as_uri() + '?mode=ro', uri=True) as conn:
        for country in bases:
            v = np.array(conn.execute('select net_econ_value,biodiversity,net_ghg_co2e from value where country=? order by sdu,scenario_id', (country,)).fetchall()).reshape(-1, 13, 3)
            mats[country] = tuple(v[:, :, k] for k in range(3))
    e, b, c = mats['Belgium']
    raw, eff = targets(bases['Belgium'])
    cells = exact_cells(e, b, c, module.envelope)
    eligible = [cell for cell in cells if cell[2][1] >= eff[1] and cell[2][2] >= eff[2]]
    assert eligible
    best = max(eligible, key=lambda cell: cell[2][0])
    climate_only = {'maximum_E_subject_to_B_C': float(best[2][0]), 'exact_maximum_E': str(best[2][0]),
                    'coefficient_cell': [str(best[0]), str(best[1])],
                    'minimum_effective_floor_relaxation': float(eff[0] - best[2][0]),
                    'fraction_of_original_floor': float((eff[0] - best[2][0]) / raw[0])}
    print('Belgian climate-only ceiling', climate_only, flush=True)
    regions = json.loads(args.intervals.read_text())['regions_by_ambition']['0.05']
    bounds = []
    for country in ['Slovenia', 'Austria', 'Bulgaria', 'Portugal', 'Lithuania', 'Czechia']:
        endpoint = regions[country][0][0]
        cutoff = Fraction(int(endpoint['numerator']), int(endpoint['denominator']))
        record, upper = regional_upper(e, b, c, cutoff, eff)
        record['threshold_country'] = country
        record['fraction_of_original_floor'] = float((eff[0] - upper) / raw[0])
        bounds.append(record)
        print('Regional bound', country, record['upper_E'], record['fraction_of_original_floor'], flush=True)
    lam = 0.27
    policies, checks, unresolved = [], [], []
    for country, (e, b, c) in mats.items():
        ratio = 0.97 if country == 'Belgium' else 1.0
        raw, eff = targets(bases[country], floor_ratio=ratio)
        tests = band_checks(e, b, c, lam, 0.0, raw, eff)
        beta = 0.0
        if not all(x['compliant'] for x in tests):
            cells = exact_cells(e, b, c, module.envelope, fixed_climate=lam)
            candidates = [cell for cell in cells if all(v >= t for v, t in zip(cell[2], eff))]
            candidates.sort(key=lambda cell: cell[2][0], reverse=True)
            found = False
            for left, right, guaranteed in candidates:
                candidate = float((left + right) / 2 if right is not None else left + 1)
                tests = band_checks(e, b, c, lam, candidate, raw, eff)
                if all(x['compliant'] for x in tests):
                    beta, found = candidate, True
                    break
            if not found:
                unresolved.append({'country': country, 'exact_compliant_cells': len(candidates)})
                print('No band-compliant candidate', country, len(candidates), flush=True)
                continue
        policies.append({'country': country, 'lambda_climate': lam, 'beta_B': beta, 'Belgian_floor_relaxation': .03, 'country_floor_ratio': ratio})
        checks.extend({'country': country, **row} for row in tests)
        print('Policy', country, beta, min(row['minimum_normalized_margin'] for row in tests), flush=True)
    assert len(policies) == 23 and not unresolved
    # Coverage bounds use the original regions only after independently
    # excluding every Belgian climate-only policy at the relaxed floors.
    assert targets(bases['Belgium'], floor_ratio=.95)[1][0] > best[2][0]
    by_country = {row['country']: row for row in policies}
    gspec = importlib.util.spec_from_file_location('price_groups', repo / 'code/price_groups.py')
    groups = importlib.util.module_from_spec(gspec)
    gspec.loader.exec_module(groups)
    austria = regions['Austria'][0][0]
    cutoff = Fraction(int(austria['numerator']), int(austria['denominator']))
    new_frontier = groups.constrained_coverage(regions, cutoff=cutoff)
    exact_regions = {country: [(groups.endpoint(lo), groups.endpoint(hi)) for lo, hi in spans]
                     for country, spans in regions.items()}
    endpoints = sorted({Fraction(0), cutoff} | {v for spans in exact_regions.values() for pair in spans for v in pair if v is not None})
    names = sorted(exact_regions)
    masks, low_masks = set(), set()
    for left, right in zip(endpoints, endpoints[1:] + [None]):
        point = (left + right) / 2 if right is not None else left + 1
        mask = sum(1 << k for k, country in enumerate(names)
                   if any(lo < point and (hi is None or point < hi) for lo, hi in exact_regions[country]))
        masks.add(mask)
        if point < cutoff:
            low_masks.add(mask)
    states = set(low_masks)
    full_policies, full_checks, frontier_results = [], [], []
    for row in new_frontier:
        q = row['maximum_climate_coefficients']
        if q > 1:
            states = {previous | new for previous in states for new in masks}
        independent_bound = 23 - max(mask.bit_count() for mask in states)
        assert independent_bound == row['complement_lower_bound']
        prices = [lam] + row['coverage_witness_prices'][1:]
        block = []
        for country, (e, b, c) in mats.items():
            raw_c, eff_c = targets(bases[country], floor_ratio=.97 if country == 'Belgium' else 1.0)
            selected = None
            for price in prices:
                tests = band_checks(e, b, c, price, 0.0, raw_c, eff_c)
                if all(r['compliant'] for r in tests):
                    selected = {'q': q, **by_country[country], 'lambda_climate': price, 'beta_B': 0.0}
                    break
            if selected is None:
                selected = {'q': q, **by_country[country]}
                tests = band_checks(e, b, c, selected['lambda_climate'], selected['beta_B'], raw_c, eff_c)
            assert all(r['compliant'] for r in tests)
            block.append(selected)
            full_checks.extend({'q': q, 'country': country, **r} for r in tests)
        count = sum(r['beta_B'] > 0 for r in block)
        assert count == independent_bound
        assert len({r['lambda_climate'] for r in block}) <= q
        full_policies.extend(block)
        frontier_results.append({'q': q, 'lower_bound': independent_bound, 'attained_upper_bound': count,
                                 'exact': True, 'all_three_response_bands': True})
    # A relaxed floor leaves the original feasible constructions feasible.
    # The low-relaxation certificate and all-price Belgian climate-only ceiling
    # retain their original lower bounds at 0.5 percent. At 3 and 5 percent,
    # the Austrian cutoff and Belgian climate-only ceiling give the new bounds.
    sensitivity = []
    for relaxation in (0.0, .005, .01, .03, .05):
        eff_E = targets(bases['Belgium'], floor_ratio=1-relaxation)[1][0]
        assert eff_E > best[2][0]
        bound_row = next(r for r in bounds if r['threshold_country'] == ('Slovenia' if relaxation <= .005 else 'Austria'))
        assert eff_E > Fraction(bound_row['upper_E_exact'])
        for q, original_k, revised in zip(range(1, 8), [20, 14, 11, 9, 7, 5, 4], frontier_results):
            lower = original_k if relaxation <= .005 else revised['lower_bound']
            upper = original_k if relaxation < .03 else revised['attained_upper_bound']
            sensitivity.append({'Belgian_floor_relaxation': relaxation, 'q': q, 'lower_bound': lower,
                                'upper_bound': upper, 'exact': lower == upper, 'all_three_response_bands': True})
    report = {'source_sha256': SHA, 'biodiversity_ambition': .05,
              'Belgian_climate_only_ceiling': climate_only,
              'regional_bounds': bounds, 'fixed_common_rate': lam,
              'complete_country_count': len(policies),
              'complement_count': sum(row['beta_B'] > 0 for row in policies),
              'all_bands_compliant': len(policies) == 23 and all(row['compliant'] for row in checks),
              'unresolved_countries': unresolved,
              'relaxed_floor_exact_frontier': frontier_results,
              'relaxed_frontier_policy_rows': len(full_policies),
              'relaxed_frontier_response_checks': len(full_checks),
              'scope_of_exact_relaxed_frontier': 'Belgian floor relaxation 3 percent or 5 percent, all other national safeguards unchanged',
              'source_unchanged': hashlib.sha256(db.read_bytes()).hexdigest() == SHA}
    (out / 'safeguard_analysis.json').write_text(json.dumps(report, indent=2) + '\n')
    write_csv(out / 'safeguard_policies.csv', policies)
    write_csv(out / 'safeguard_policy_checks.csv', checks)
    write_csv(out / 'safeguard_bounds.csv', bounds)
    write_csv(out / 'safeguard_frontier.csv', sensitivity)
    write_csv(out / 'safeguard_frontier_policies.csv', full_policies)
    write_csv(out / 'safeguard_frontier_checks.csv', full_checks)
    print(json.dumps({key: value for key, value in report.items() if key not in ['regional_bounds']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--database', type=Path, help='Optional exact uncompressed source database')
    parser.add_argument('--intervals', type=Path, help='Optional freshly reproduced exact interval JSON')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='safeguard_source_') as temporary:
        temp = Path(temporary)
        if args.database is None:
            args.database = temp / 'eu23.sqlite'
            parts = [args.repository / 'inputs' / f'eu23_coefficients.sqlite.gz.part{i:02d}' for i in (1, 2)]
            args.database.write_bytes(gzip.decompress(b''.join(p.read_bytes() for p in parts)))
        if args.intervals is None:
            source_output = temp / 'intervals'
            with (temp / 'interval_replay.log').open('w') as log:
                subprocess.run([sys.executable, str(args.repository / 'code/verify_breadth.py'),
                                '--database', str(args.database), '--output', str(source_output)],
                               check=True, stdout=log, stderr=subprocess.STDOUT)
            args.intervals = source_output / 'exact_climate_only_intervals.json'
        main(args)
