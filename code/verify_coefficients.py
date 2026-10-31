"""Verify fixed-coefficient policies, dominance, support and resource certificates.
No optimization routines are imported or executed.
"""

from pathlib import Path
import sqlite3, zipfile, io, json, math, time, hashlib, os
import numpy as np, pandas as pd

P = Path(__file__).resolve().parents[1]
R = P / "certificates"
O = Path(os.environ.get("CBC_VERIFY_OUTPUT", str(P / "outputs/matrix_replay")))


def S(x):
    return math.fsum(np.asarray(x).ravel().tolist())


def read(n):
    return pd.read_csv(R / n, float_precision="round_trip")


def close(a, b):
    assert abs(a - b) <= max(0.02, 5e-11 * max(1, abs(a), abs(b))), (a, b)


def dump(n, v):
    (O / n).write_text(json.dumps(v, indent=2, ensure_ascii=False), encoding="utf-8")


def main():
    start = time.time()
    base = pd.read_csv(
        R / "Frozen_Baselines.csv", float_precision="round_trip"
    ).set_index("country")
    assert len(base) == 23 and base.cee11.sum() == 11
    db = sqlite3.connect(
        Path(os.environ["CBC_MATRIX_DATABASE"]).resolve().as_uri() + "?mode=ro",
        uri=True,
    )
    order = db.execute("select * from scenario order by scenario_id").fetchall()
    assert order[0] == (0, "sustainable_current")
    dump("scenario_order.json", order)
    data = {}
    for c, b in base.iterrows():
        x = np.array(
            db.execute(
                "SELECT sdu,scenario_id,net_econ_value,biodiversity,net_ghg_co2e FROM value WHERE country=? ORDER BY sdu,scenario_id",
                (c,),
            ).fetchall()
        ).reshape(-1, 13, 5)
        assert len(x) == b.n_sdu and np.array_equal(
            x[:, :, 1], np.tile(np.arange(13), (len(x), 1))
        )
        e, bb, cc = x[:, :, 2], x[:, :, 3], x[:, :, 4]
        d = dict(
            E=e,
            B=bb,
            C=cc,
            G=np.stack([e, bb, cc]),
            A=np.maximum(cc - cc[:, 0, None], 0),
            D=cc - cc[:, 0, None],
            H=np.maximum(bb - bb[:, 0, None], 0),
            scale=np.maximum(1, np.max(abs(e), axis=1)),
            target=np.array([b.E0, b.B0 + 0.05 * (b.Bmax - b.B0), b.C0]),
        )
        data[c] = d
    db.close()
    arrays = dict(np.load(R / "assignments.npz"))
    rows = read("policies.csv")
    max_cash_error = 0.0
    max_gap = 0.0
    checks = []
    for rr in rows.itertuples():
        d = data[rr.country]
        s = arrays[rr.assignment].astype(int)
        ii = np.arange(len(s))
        a = d["A"] if rr.climate == "positive" else d["D"]
        h = (
            d["B"]
            if rr.bio == "level"
            else d["H"] if rr.bio == "gain" else d["B"] - d["B"][:, 0, None]
        )
        pc = rr.lambda_climate * a
        pb = rr.beta_B * h
        u = d["E"] + pc + pb
        ma = np.max(u, axis=1)
        gap = float(np.max(ma - u[ii, s]))
        assert gap <= 1e-8 * max(1, float(np.max(abs(u))))
        max_gap = max(max_gap, gap)
        v = np.array([S(g[ii, s]) for g in d["G"]])
        for val, stated in zip(v, [rr.E, rr.B, rr.C20]):
            close(val, stated)
        gross = S(np.maximum(pc[ii, s], 0)) + S(np.maximum(pb[ii, s], 0))
        receipt = S(np.maximum(-pc[ii, s], 0)) + S(np.maximum(-pb[ii, s], 0))
        net = S(pc[ii, s] + pb[ii, s])
        close(gross, rr.gross_component_disbursements)
        close(receipt, rr.gross_component_receipts)
        close(net, rr.net_transfer)
        max_cash_error = max(max_cash_error, abs(gross - receipt - net))
        close(S(np.max(d["E"], axis=1)) - v[0], rr.resource_opportunity_cost)
        for tag, eps in [("exact", 0.0), ("stress6", 1e-6), ("stress4", 1e-4)]:
            mask = (
                (u == ma[:, None])
                if eps == 0
                else ma[:, None] - u <= eps * d["scale"][:, None]
            )
            low = np.array(
                [S(np.min(g, where=mask, initial=np.inf, axis=1)) for g in d["G"]]
            )
            high = np.array(
                [S(np.max(g, where=mask, initial=-np.inf, axis=1)) for g in d["G"]]
            )
            t = d["target"]
            ok = bool(np.all(low >= t - 1e-8 * np.maximum(1, abs(t))))
            assert ok == getattr(rr, tag + "_robust"), (
                rr.country,
                rr.label,
                tag,
                low - t,
            )
            for k, g in enumerate(["E", "B", "C20"]):
                close(low[k], getattr(rr, tag + "_min_" + g))
                close(high[k], getattr(rr, tag + "_max_" + g))
            checks.append(
                dict(
                    policy=rr.label,
                    country=rr.country,
                    assignment=rr.assignment,
                    band=eps,
                    robust=ok,
                    minimum_normalized_margin=float(
                        np.min((low - t) / np.maximum(1, abs(t)))
                    ),
                    net_min=S(np.min(pc + pb, where=mask, initial=np.inf, axis=1)),
                    net_max=S(np.max(pc + pb, where=mask, initial=-np.inf, axis=1)),
                )
            )
    pd.DataFrame(checks).to_csv(O / "verified_policy_bands.csv", index=False)
    # Verify every dominance witness, using the original three features, not the
    # optimization helper's mask implementation.
    dom = dict(np.load(R / "dominance_witnesses.npz"))
    domcount = 0
    for key, m in dom.items():
        c, cl, bi = key.split("__")
        d = data[c]
        a = d["A"] if cl == "positive" else d["D"]
        h = (
            d["B"]
            if bi == "level"
            else d["H"] if bi == "gain" else d["B"] - d["B"][:, 0, None]
        )
        ii, ss = np.where(m >= 0)
        tt = m[ii, ss]
        assert np.all(d["E"][ii, tt] > d["E"][ii, ss])
        assert np.all(a[ii, tt] >= a[ii, ss]) and np.all(h[ii, tt] >= h[ii, ss])
        domcount += len(ss)
    certs = []
    for rr in json.loads((R / "screen_bounds.json").read_text()):
        if "mu" not in rr:
            continue
        c = rr["country"]
        d = data[c]
        keep = dom[c + "__" + rr["climate"] + "__" + rr["bio"]] < 0
        mu = np.array(rr["mu"])
        assert np.min(mu) >= 0 and mu.sum() <= 1 + 1e-12
        sc = np.array(rr["scales"])
        eff = d["target"] - 1e-8 * np.maximum(1, abs(d["target"]))
        weighted = np.einsum("k,kij->ij", mu / sc, d["G"])
        maxima = np.max(weighted, where=keep, initial=-np.inf, axis=1)
        lower = float(np.dot(mu, eff / sc) - S(maxima)) - 2e-9
        certs.append(
            dict(
                country=c,
                climate=rr["climate"],
                bio=rr["bio"],
                verified_shortfall_lower=lower,
                all_finite_prices_excluded=bool(lower > 1e-8),
            )
        )
    pd.DataFrame(certs).to_csv(O / "verified_exclusions.csv", index=False)
    # Validate the support of each delivered interior region in homogeneous
    # coordinates, including unchosen original alternatives.
    regionchecked = 0
    for rr in json.loads((R / "support_regions.json").read_text()):
        if not rr.get("region_assignment"):
            continue
        c = rr["country"]
        d = data[c]
        s = arrays[rr["region_assignment"]].astype(int)
        ii = np.arange(len(s))
        a = d["A"] if rr["climate"] == "positive" else d["D"]
        h = (
            d["B"]
            if rr["bio"] == "level"
            else d["H"] if rr["bio"] == "gain" else d["B"] - d["B"][:, 0, None]
        )
        q = np.array([rr["r"], rr["u"], rr["v"]])
        assert min(q) >= -1e-8 and abs(q.sum() - 1) < 1e-7
        D = np.stack(
            [
                d["E"][ii, s, None] - d["E"],
                rr["L"] * (a[ii, s, None] - a),
                rr["H"] * (h[ii, s, None] - h),
            ],
            axis=2,
        ).reshape(-1, 3)
        nz = np.max(abs(D), axis=1) > 0
        D = D[nz]
        D /= np.max(abs(D), axis=1)[:, None]
        viol = np.max(-D @ q + rr["margin"] * np.sum(abs(D), axis=1))
        assert viol < 5e-8, (c, rr["label"], viol)
        regionchecked += 1
    # Certify supplied scalar hulls against ALL 13 alternatives. At either endpoint
    # an affine support inequality implies it throughout each open interval; at the
    # unbounded tail, check slopes explicitly. Persistent ecological ties are complete.
    hulls = dict(np.load(R / "hull_certificates.npz"))
    hrows = 0
    worst_norm = 0.0
    assert set(hulls) == {
        c + "__" + kind for c in data for kind in ["carbon", "biodiversity"]
    }
    for key, hh in hulls.items():
        c, kind = key.split("__")
        d = data[c]
        e = d["E"] if kind == "carbon" else d["E"] + 0.1 * d["A"]
        z = d["A"] if kind == "carbon" else d["B"]
        starts = np.r_[np.flatnonzero(np.r_[True, hh[1:, 0] != hh[:-1, 0]]), len(hh)]
        assert hh.ndim == 2 and hh.shape[1] == 5 and np.isfinite(hh).all()
        assert np.array_equal(np.unique(hh[:, 0]), np.arange(len(e)))
        assert np.all(np.diff(hh[:, 0]) >= 0)
        for lo, hi in zip(starts[:-1], starts[1:]):
            vv = hh[lo:hi]
            i = int(vv[0, 0])
            assert np.array_equal(vv[:, 1], np.arange(len(vv)))
            assert vv[0, 2] == 0 and np.all(np.diff(vv[:, 2]) > 0)
            assert np.all(vv[:, 3] == vv[:, 3].astype(int)) and np.all(
                (vv[:, 3] >= 0) & (vv[:, 3] < 13)
            )
            for k, (un, ind, q, s, mask) in enumerate(vv):
                s = int(s)
                opts = np.flatnonzero((e[i] == e[i, s]) & (z[i] == z[i, s]))
                assert int(mask) == sum(1 << int(j) for j in opts)
                points = [float(q)]
                if k + 1 < len(vv):
                    points.append(float(vv[k + 1, 2]))
                else:
                    assert np.all(z[i, s] >= z[i])
                for price in points:
                    diff = (e[i, s] - e[i]).astype(np.longdouble) + np.longdouble(
                        price
                    ) * (z[i, s] - z[i]).astype(np.longdouble)
                    guard = 1e-9 * max(1, float(np.max(abs(e[i] + price * z[i]))))
                    assert np.min(diff) >= -guard, (
                        c,
                        kind,
                        i,
                        price,
                        float(np.min(diff)),
                    )
                    worst_norm = min(worst_norm, float(np.min(diff)) / guard)
                hrows += 1
    # Matched resource certificates: source x is a verified implemented integer
    # vector. Upper/lower resource certificates are recomputed without the solver.
    matched = []
    for name in ["matched_alternative.csv", "matched_final_alternative.csv"]:
        for rr in read(name).itertuples():
            d = data[rr.country]
            s = arrays[rr.assignment].astype(int)
            ref = arrays[rr.reference_assignment].astype(int)
            ii = np.arange(len(s))
            rv = np.array([S(g[ii, ref]) for g in d["G"]])
            target = np.maximum(d["target"], [d["target"][0], rv[1], rv[2]])
            v = np.array([S(g[ii, s]) for g in d["G"]])
            assert np.min((v - target) / np.maximum(1, abs(target))) >= -1.0001e-8
            mu = np.array([rr.mu_E, rr.mu_B, rr.mu_C])
            assert np.min(mu) >= 0
            cost = np.max(d["E"], axis=1)[:, None] - d["E"]
            eff = target - 1e-8 * np.maximum(1, abs(target))
            w = np.einsum("k,kij->ij", mu, d["G"])
            guard = max(
                30.0,
                256
                * np.finfo(float).eps
                * (
                    S(
                        np.max(
                            abs(cost) + np.einsum("k,kij->ij", abs(mu), abs(d["G"])),
                            axis=1,
                        )
                    )
                    + S(abs(mu * eff))
                ),
            )
            lower = max(0, S(mu * eff) + S(np.min(cost - w, axis=1)) - guard)
            upper = S(cost[ii, s])
            assert lower <= upper + 1
            matched.append(
                dict(
                    country=rr.country,
                    policy=rr.label,
                    resource_lower=lower,
                    resource_upper=upper,
                    resource_gain_lower=v[0] - rv[0],
                    resource_gain_upper=S(np.max(d["E"], axis=1)) - lower - rv[0],
                    reference_E=rv[0],
                    planner_E=v[0],
                    reference_B=rv[1],
                    planner_B=v[1],
                    reference_C20=rv[2],
                    planner_C20=v[2],
                )
            )
    pd.DataFrame(matched).to_csv(O / "verified_matched.csv", index=False)
    # Fixed-allocation recentering checks component receipts BEFORE netting.
    for rr in read("fixed_allocation_accounts.csv").itertuples():
        d = data[rr.country]
        x = arrays[rr.assignment].astype(int)
        ii = np.arange(len(x))
        pc = rr.lambda_climate * d["A"][ii, x]
        pb = rr.beta_B * (d["B"][ii, x] - d["B"][:, 0])
        gross = S(np.maximum(pc, 0)) + S(np.maximum(pb, 0))
        receipt = S(np.maximum(-pc, 0)) + S(np.maximum(-pb, 0))
        close(gross, rr.signed_B_gross)
        close(receipt, rr.signed_B_receipts)
        close(gross - receipt, rr.signed_B_net)
        close(rr.level_net - rr.signed_B_net, rr.beta_B * S(d["B"][:, 0]))
    # Final robust constructions: count instruments from coefficients, not labels.
    final = read("final_EU23.csv")
    counts = {}
    for label, g in final.groupby("label"):
        assert (
            len(g) == 23
            and g.exact_robust.all()
            and g.stress6_robust.all()
            and g.stress4_robust.all()
        )
        counts[label] = dict(
            countries=23,
            positive_B=int(np.sum(g.beta_B > 0)),
            CEE_positive_B=int(np.sum((g.beta_B > 0) & g.CEE11)),
            guaranteed_band=1e-4,
        )
    assert (
        counts["national_final"]["positive_B"] == 4
        and counts["common_final"]["positive_B"] == 21
    )
    summary = dict(
        status="PASS",
        seconds=time.time() - start,
        optimizer_calls=0,
        policy_rows_checked=len(rows),
        policy_band_checks=len(checks),
        dominance_witness_pairs_verified=domcount,
        all_price_exclusions=[r for r in certs if r["all_finite_prices_excluded"]],
        region_certificates_checked=regionchecked,
        scalar_support_segments_checked=hrows,
        support_domain_completeness_checked=True,
        largest_relative_endpoint_guard_use=worst_norm,
        resource_certificates_checked=len(matched),
        maximum_cash_identity_error=max_cash_error,
        maximum_canonical_choice_gap=max_gap,
        robust_constructions=counts,
        common_complements_constructed=counts["common_final"]["positive_B"],
        national_complements_constructed=counts["national_final"]["positive_B"],
        coverage_bounds="Computed separately from original coefficients by verify_breadth.py",
    )
    dump("verification.json", summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
