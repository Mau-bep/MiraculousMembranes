#!/usr/bin/env python
"""
extract_lines.py -- phase lines of the planar bead simulations in Deserno's variables, versus the theory.

Reads the PostProcessing output (Coverage_data*.txt) of the planar-excess campaigns (flat start and full start),
converts to Deserno's variables and extracts, for every KA (sigma~ = 2 KA a^2/KB), as functions of
w~ = 4 KI a^2/KB

    S1_sim : first-order binding jump of the FLAT start (lowest w~ at which the flat start has jumped to the
             enveloped state, z = 2*CoverageUnion above a threshold for all larger KI);  ~ spinodal S1
    S2_sim : lowest w~ at which the FULL start stays enveloped (it has unwrapped below);   ~ spinodal S2
    E_sim  : w~ at which Total_E(flat start) = Total_E(full start) (both states distinct), linear (and quadratic)
             interpolation in KI

and compares them with the theory lines of Scripts/deserno_theory.py (spinodal_S1, spinodal_S2, w_E).

Importable functions (all pure numpy):
    load_coverage(path)                  -> dict of columns + bookkeeping (duplicates, NaNs, constants)
    build_grid(flat, full)               -> {KAkey: series dict on the KI grid}
    lowest_persistent(KI, z, thr)        -> bracket of the lowest KI from which z >= thr for all larger KI
    largest_increment(KI, z)             -> bracket of the largest increment of z between neighbouring KI
    energy_crossing(KI, Ef, Eu, zf, zu)  -> E-line estimate with bracket and interpolation error
    extract_all(flat, full)              -> list of row dicts (one per KA) = baseline_lines.csv
    hellmann_feynman(series)             -> consistency of Total_E(KI) with Bead/KI (envelope theorem)
    enveloped_identity(flat, full)       -> energy of the enveloped state against -2(w~-4)+4 sigma~
    partial_branch(...)                  -> energy and z of the partially wrapped state against the theory
Command line:   python extract_lines.py [--out-dir DIR]     (needs ../data/theory_cache.npz for the energy part,
                                                              built by theory_cache_build.py)
Conventions: kappa = KB/2, sigma = KA, w = KI, E~ = E/(pi kappa) = 2 E/(pi KB), z = 1 - cos(alpha) = 2 CoverageUnion.
"""
import os
import sys
import csv
import argparse
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
P = STUDY.parents[1]                                   # projects/geometric-flow
OLD = Path("/home/mrojasve/Documents/DDG/MiraculousMembranes/projects/geometric-flow/Results/Wrapping_planar_excess")
FLAT_FILE = OLD / "Coverage_data.txt"
FULL_FILE = OLD / "Coverage_data_full.txt"
sys.path.insert(0, str(P / "Scripts"))

COLUMNS = ["DIR", "KA", "KB", "KI", "BeadRadius", "rc", "CoveredArea", "CoverageUnion", "MultilayerFrac", "BeadX",
           "Area", "Excess_tension", "Bending_tan", "Bead", "Total_E"]

# ---- analysis parameters -------------------------------------------------------------------------------------------
Z_THR = 1.75        # z above which a state counts as "enveloped" (z = 2*CoverageUnion; enveloped plateau is 1.94-1.99)
DZ_SHARP = 0.30     # a jump counts as sharp (first-order-like) if the largest increment of z between neighbouring KI >= this
DZ_DISTINCT = 0.20  # flat/full states count as different states if |z_full - z_flat| > this (same-branch scatter is below)
A0_FLAT = 314.138   # reference area in both campaigns (checked in check_reference_area)


# =====================================================================================================================
# 1. Loading
# =====================================================================================================================
def load_coverage(path):
    """Read a Coverage_data file written by PostProcessing.  Returns dict with arrays (one entry per kept row) and
    a 'info' dict: n_rows, n_duplicates (later rows replace earlier ones, as PostProcessing appends), n_nonfinite
    (dropped), constants (set of KB, BeadRadius, rc)."""
    rows, nonfinite = [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            if len(t) != len(COLUMNS):
                nonfinite += 1
                continue
            v = [float(x) for x in t[1:]]
            if not np.all(np.isfinite(v)):
                nonfinite += 1
                continue
            rows.append((t[0], v))
    seen, ndup = {}, 0
    for i, (d, v) in enumerate(rows):
        key = (round(v[0], 3), round(v[2], 2))
        if key in seen:
            ndup += 1
        seen[key] = i                                      # keep the last occurrence
    keep = sorted(seen.values())
    out = {"DIR": np.array([rows[i][0] for i in keep])}
    for j, name in enumerate(COLUMNS[1:]):
        out[name] = np.array([rows[i][1][j] for i in keep])
    out["key"] = np.round(out["KA"], 3)
    out["KIr"] = np.round(out["KI"], 2)
    out["info"] = dict(n_rows=len(rows), n_kept=len(keep), n_duplicates=ndup, n_nonfinite=nonfinite,
                       constants=dict(KB=sorted(set(out["KB"])), a=sorted(set(out["BeadRadius"])),
                                      rc=sorted(set(out["rc"]))))
    return out


def deserno_vars(ds):
    """w~, sigma~, E~ (units pi*kappa), z (cov), z_ad (adhesion-weighted covered solid angle / 2 pi)."""
    a, KB = ds["BeadRadius"], ds["KB"]
    w = 4.0 * ds["KI"] * a ** 2 / KB
    s = 2.0 * ds["KA"] * a ** 2 / KB
    E = 2.0 * ds["Total_E"] / (np.pi * KB)
    z = 2.0 * ds["CoverageUnion"]
    zad = -ds["Bead"] / (2.0 * np.pi * ds["KI"] * a ** 2)
    return w, s, E, z, zad


def build_grid(flat, full):
    """Series per KA on the common KI grid; missing (KA, KI) combinations (unfinished runs) are NaN."""
    KIs = np.array(sorted(set(flat["KIr"]) | set(full["KIr"])))
    keys = sorted(set(flat["key"]) | set(full["key"]))
    names = ["CoverageUnion", "Total_E", "Bead", "Bending_tan", "Excess_tension", "Area", "BeadX", "KA"]
    grid = {}
    for k in keys:
        g = {"key": k, "KI": KIs}
        for tag, ds in (("f", flat), ("u", full)):
            m = ds["key"] == k
            idx = np.searchsorted(KIs, ds["KIr"][m])
            for n in names:
                arr = np.full(len(KIs), np.nan)
                arr[idx] = ds[n][m]
                g[n + "_" + tag] = arr
            g["a"] = ds["BeadRadius"][m][0]
            g["KB"] = ds["KB"][m][0]
        grid[k] = g
    return grid


# =====================================================================================================================
# 2. Line extraction
# =====================================================================================================================
def lowest_persistent(KI, z, thr=Z_THR):
    """Lowest KI from which z >= thr for ALL larger KI (non-finite entries skipped).  Returns dict(status, lo, hi,
    mid): status 'ok' (bracket (lo, hi) = last KI below thr, first KI of the persistent enveloped run), 'low' (already
    above thr at the smallest KI), 'high' (never persistently above thr up to the largest KI; lo = largest KI)."""
    f = np.isfinite(z)
    K, Z = KI[f], z[f]
    n = len(Z)
    i = n
    while i > 0 and Z[i - 1] >= thr:
        i -= 1
    if i == n:
        return dict(status="high", lo=K[-1], hi=np.nan, mid=np.nan)
    if i == 0:
        return dict(status="low", lo=np.nan, hi=K[0], mid=np.nan)
    return dict(status="ok", lo=K[i - 1], hi=K[i], mid=0.5 * (K[i - 1] + K[i]))


def largest_increment(KI, z):
    """Largest increment of z between neighbouring (finite) KI: returns dict(lo, hi, mid, inc)."""
    f = np.isfinite(z)
    K, Z = KI[f], z[f]
    d = np.diff(Z)
    i = int(np.argmax(d))
    return dict(lo=K[i], hi=K[i + 1], mid=0.5 * (K[i] + K[i + 1]), inc=d[i], i=i)


def energy_crossing(KI, Ef, Eu, zf, zu, s_lo, s_hi, sigE=0.05, dz_distinct=DZ_DISTINCT):
    """Line E: KI at which Total_E(flat) = Total_E(full).  Only points where the two runs are in different states
    (z_full - z_flat > dz_distinct) carry information.  dE = Ef - Eu: negative where the partially wrapped (flat
    start) state is lower, positive where the enveloped one is lower.  Returns dict with
      method 'interp' (sign change between two adjacent distinct points: linear interpolation, quadratic from up to
              three consecutive distinct points for the curvature error), or 'bracket' (no sign change seen: the
              line is only bracketed by [lo, hi] from the last negative / first positive distinct point, else by the
              hysteresis window s_lo (lowest KI where full start is enveloped, bracket lower edge) .. s_hi),
      est, lo, hi (bracket in KI), err (1 sigma-like uncertainty in KI),  n_distinct."""
    ok = np.isfinite(Ef) & np.isfinite(Eu) & np.isfinite(zf) & np.isfinite(zu)
    dist = ok & ((zu - zf) > dz_distinct)
    dE = np.where(ok, Ef - Eu, np.nan)
    idx = np.where(dist)[0]
    res = dict(n_distinct=int(len(idx)), method="none", est=np.nan, lo=np.nan, hi=np.nan, err=np.nan, quad=np.nan)
    if len(idx) == 0:
        return res
    # adjacent distinct pairs with a sign change negative -> positive
    for a, b in zip(idx[:-1], idx[1:]):
        if b == a + 1 and dE[a] < 0.0 <= dE[b]:
            lin = KI[a] + (KI[b] - KI[a]) * (-dE[a]) / (dE[b] - dE[a])
            # quadratic through up to 3 consecutive distinct points containing (a, b)
            run = [j for j in (a - 1, a, b, b + 1) if j in set(idx)]
            quad = np.nan
            if len(run) >= 3:
                c = np.polyfit(KI[run], dE[run], 2)
                r = np.roots(c)
                r = r[np.isreal(r)].real
                r = r[(r >= KI[a] - 1e-9) & (r <= KI[b] + 1e-9)]
                if len(r):
                    quad = float(r[0])
            slope = (dE[b] - dE[a]) / (KI[b] - KI[a])
            curv = abs(quad - lin) if np.isfinite(quad) else 0.0
            noise = sigE / max(abs(slope), 1e-12)
            res.update(method="interp", est=float(lin), lo=float(KI[a]), hi=float(KI[b]), quad=quad,
                       err=float(np.hypot(curv, noise)))
            return res
    # no sign change: bracket
    neg = [j for j in idx if dE[j] < 0]
    pos = [j for j in idx if dE[j] >= 0]
    lo = KI[max(neg)] if neg else s_lo
    hi = KI[min(pos)] if pos else s_hi
    if np.isfinite(lo) and np.isfinite(hi) and hi >= lo:
        res.update(method="bracket", est=0.5 * (lo + hi), lo=float(lo), hi=float(hi), err=0.5 * (hi - lo))
    return res


def extract_all(flat, full, z_thr=Z_THR, dz_sharp=DZ_SHARP, sigE=0.05):
    """One row per KA with the extracted simulated lines (KI and w~ units) and the theory lines (deserno_theory)."""
    import deserno_theory as dt
    grid = build_grid(flat, full)
    rows = []
    for k, g in grid.items():
        KI, a, KB = g["KI"], g["a"], g["KB"]
        KA = float(np.nanmean(np.r_[g["KA_f"], g["KA_u"]]))
        sig = 2.0 * KA * a ** 2 / KB
        w2 = lambda x: 4.0 * x * a ** 2 / KB                        # KI -> w~
        dw_res = w2(np.median(np.diff(KI)))                          # bracket width in w~ (0.4)
        zf, zu = 2.0 * g["CoverageUnion_f"], 2.0 * g["CoverageUnion_u"]
        zadf = -g["Bead_f"] / (2 * np.pi * KI * a ** 2)
        zadu = -g["Bead_u"] / (2 * np.pi * KI * a ** 2)
        Ef, Eu = g["Total_E_f"], g["Total_E_u"]
        S1 = lowest_persistent(KI, zf, z_thr)
        S2 = lowest_persistent(KI, zu, z_thr)
        S1b = largest_increment(KI, zf)
        S2b = largest_increment(KI, zu)
        S1c = lowest_persistent(KI, zadf, z_thr)
        S2c = lowest_persistent(KI, zadu, z_thr)
        # hysteresis present in the simulation?  flat and full start differ somewhere
        both = np.isfinite(zf) & np.isfinite(zu)
        dist = both & ((zu - zf) > DZ_DISTINCT)
        # window edges in KI used to bracket E when no sign change is seen (lower edge of S2 bracket, upper of S1)
        s_lo = S2["lo"] if S2["status"] == "ok" else np.nan
        s_hi = S1["hi"] if S1["status"] == "ok" else (KI[-1] if S1["status"] == "high" else np.nan)
        Ec = energy_crossing(KI, Ef, Eu, zf, zu, s_lo, s_hi, sigE=sigE)
        if sig < 1e-4:                  # sigma~ -> 0: all three lines emanate from the triple point w~ = 4 (outside the table)
            th = dict(S1=4.0, S2=4.0, E=4.0)
        else:
            th = dict(S1=float(dt.spinodal_S1(sig)), S2=float(dt.spinodal_S2(sig)), E=float(dt.w_E(sig)))
        r = dict(KA=KA, sigma_t=sig, n_distinct=int(dist.sum()),
                 theory_E=th["E"], theory_S1=th["S1"], theory_S2=th["S2"])
        # S1 (flat start)
        r.update(S1_status=S1["status"], S1_lo=w2(S1["lo"]), S1_hi=w2(S1["hi"]), S1_sim=w2(S1["mid"]),
                 S1_inc_lo=w2(S1b["lo"]), S1_inc_hi=w2(S1b["hi"]), S1_inc=S1b["inc"],
                 S1_sharp=bool(S1["status"] == "ok" and S1b["inc"] >= dz_sharp and abs(S1b["mid"] - S1["mid"]) <= 0.1 + 1e-9),
                 S1_sim_zad=w2(S1c["mid"]), S1_err=dw_res / 2)
        r.update(S2_status=S2["status"], S2_lo=w2(S2["lo"]), S2_hi=w2(S2["hi"]), S2_sim=w2(S2["mid"]),
                 S2_inc_lo=w2(S2b["lo"]), S2_inc_hi=w2(S2b["hi"]), S2_inc=S2b["inc"],
                 S2_sharp=bool(S2["status"] == "ok" and S2b["inc"] >= dz_sharp and abs(S2b["mid"] - S2["mid"]) <= 0.1 + 1e-9),
                 S2_sim_zad=w2(S2c["mid"]), S2_err=dw_res / 2)
        r.update(E_method=Ec["method"], E_sim=w2(Ec["est"]), E_lo=w2(Ec["lo"]), E_hi=w2(Ec["hi"]),
                 E_quad=w2(Ec["quad"]), E_err=w2(Ec["err"]) if np.isfinite(Ec["err"]) else np.nan,
                 E_half_bracket=(w2(Ec["hi"]) - w2(Ec["lo"])) / 2 if np.isfinite(Ec["hi"]) else np.nan)
        # hysteresis window in the simulation / theory and E inside?
        r["win_sim_lo"], r["win_sim_hi"] = r["S2_sim"], r["S1_sim"]
        r["win_sim_width"] = r["S1_sim"] - r["S2_sim"] if (S1["status"] == "ok" and S2["status"] == "ok") else np.nan
        r["win_th_width"] = th["S1"] - th["S2"]
        r["E_in_window_sim"] = bool(np.isfinite(r["E_sim"]) and np.isfinite(r["S1_sim"]) and np.isfinite(r["S2_sim"])
                                    and r["S2_sim"] <= r["E_sim"] <= r["S1_sim"])
        # offsets
        for nme in ("E", "S1", "S2"):
            r["d" + nme] = r[nme + "_sim"] - th[nme]
        rows.append(r)
    return rows


def write_csv(rows, path):
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(keys)
        for r in rows:
            wr.writerow([("%.6g" % r[k]) if isinstance(r[k], (float, np.floating)) else r[k] for k in keys])


# =====================================================================================================================
# 3. Data-quality diagnostics
# =====================================================================================================================
def check_reference_area(flat, full):
    """Implied A0 = Area - Excess/KA for rows with Excess > 0 (KA rounded to 4 digits in the files gives the scatter)."""
    out = {}
    for nm, ds in (("flat", flat), ("full", full)):
        m = (ds["Excess_tension"] > 0) & (ds["KA"] > 0)
        A0 = ds["Area"][m] - ds["Excess_tension"][m] / ds["KA"][m]
        out[nm] = (float(A0.min()), float(np.median(A0)), float(A0.max()))
        out[nm + "_sumcheck"] = float(np.abs(ds["Total_E"] - ds["Excess_tension"] - ds["Bending_tan"] - ds["Bead"]).max())
    return out


def same_branch_scatter(grid):
    """Scatter between flat and full start where both are in the same state (|dz| small): robust sigma of
    dz = z_flat - z_full and dE = Total_E_flat - Total_E_full (code units and E~ units)."""
    dz, dE, dEt, tags = [], [], [], []
    for k, g in grid.items():
        zf, zu = 2 * g["CoverageUnion_f"], 2 * g["CoverageUnion_u"]
        ok = np.isfinite(zf) & np.isfinite(zu) & (np.abs(zf - zu) < DZ_DISTINCT)
        dz += list(zf[ok] - zu[ok])
        d = (g["Total_E_f"] - g["Total_E_u"])[ok]
        dE += list(d)
        dEt += list(2 * d / (np.pi * g["KB"]))
        tags += [(k, KI) for KI in g["KI"][ok]]
    dz, dE, dEt = map(np.array, (dz, dE, dEt))
    mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
    return dict(n=len(dz), dz_mad=mad(dz), dz_p95=np.percentile(np.abs(dz), 95), dE_mad=mad(dE), dE_p95=np.percentile(np.abs(dE), 95),
                dEt_mad=mad(dEt), dE_max=np.abs(dE).max(), dEt_p95=np.percentile(np.abs(dEt), 95))


def hellmann_feynman(grid, which="f"):
    """At a stationary state dE/dKI = Bead/KI (KI enters only through the adhesion prefactor).  Compare
    Delta Total_E with the trapezoid of Bead/KI between neighbouring KI on the same branch (both z < thr or both >=
    thr and |dz| < DZ_SHARP).  Returns residual array (code units) and the ratio to |Delta Total_E|."""
    res, rel, info = [], [], []
    for k, g in grid.items():
        E, B, KI = g["Total_E_" + which], g["Bead_" + which], g["KI"]
        z = 2 * g["CoverageUnion_" + which]
        for i in range(len(KI) - 1):
            if not np.all(np.isfinite([E[i], E[i + 1], B[i], B[i + 1], z[i], z[i + 1]])):
                continue
            same = (z[i] < Z_THR) == (z[i + 1] < Z_THR) and abs(z[i + 1] - z[i]) < DZ_SHARP
            if not same:
                continue
            dEobs = E[i + 1] - E[i]
            dEhf = 0.5 * (B[i] / KI[i] + B[i + 1] / KI[i + 1]) * (KI[i + 1] - KI[i])
            res.append(dEobs - dEhf)
            rel.append((dEobs - dEhf) / abs(dEobs))
            info.append((g["key"], KI[i], z[i]))
    return np.array(res), np.array(rel), info


# =====================================================================================================================
# 4. Energy of the enveloped state
# =====================================================================================================================
def enveloped_rows(flat, full, zmin=1.9):
    """All rows of both files in the enveloped state (z_cov >= zmin) with the decomposition in E~ units."""
    out = []
    for nm, ds in (("flat", flat), ("full", full)):
        w, s, E, z, zad = deserno_vars(ds)
        a, KB = ds["BeadRadius"], ds["KB"]
        m = z >= zmin
        for i in np.where(m)[0]:
            Bend_t = 2 * ds["Bending_tan"][i] / (np.pi * KB[i])
            Bead_t = 2 * ds["Bead"][i] / (np.pi * KB[i])
            Exc_t = 2 * ds["Excess_tension"][i] / (np.pi * KB[i])
            dA = ds["Area"][i] - A0_FLAT
            th = -2 * (w[i] - 4) + 4 * s[i]
            out.append(dict(file=nm, KA=ds["KA"][i], KI=ds["KI"][i], w=w[i], sigma=s[i], z=z[i], zad=zad[i], E=E[i], Eth=th,
                            resid=E[i] - th, bend=Bend_t, adh=Bead_t, tens=Exc_t,
                            wavg=-ds["Bead"][i] / (4 * np.pi * ds["KI"][i] * a[i] ** 2), dA_over_pi=dA / np.pi,
                            resid_bend=Bend_t - 8.0, resid_adh=Bead_t + 2 * w[i],
                            resid_tens=Exc_t - 4 * s[i]))
    return out


def shell_offset_from_weight(wavg, s_frac=0.25, a=1.0):
    """Equivalent uniform radial offset delta of a surface sitting at r = a + delta in the shell
    w(r) = cos^2(pi delta / (2 s a)):  delta = (2 s a / pi) arccos(sqrt(w))."""
    return (2 * s_frac * a / np.pi) * np.arccos(np.sqrt(np.clip(wavg, 0, 1)))


# =====================================================================================================================
# 5. Partially wrapped branch against the theory (needs data/theory_cache.npz)
# =====================================================================================================================
def load_theory_cache(path=None):
    path = path or (STUDY / "data" / "theory_cache.npz")
    if not Path(path).exists():
        return None
    c = np.load(path)
    return {k: c[k] for k in c.files}


def theory_minima(cache, k_idx, w):
    """Local minima of E~(z) = -(w-4) z + sigma z^2 + F(z; sigma) over z in [0, 2] (cache at sigma = 2 k/15).  Returns
    list of (z, E) sorted by z; endpoints z = 0 and z = 2 count when they are minima.  The lowest-energy branch
    that free_energy() returns beyond the energy crossing is branch 3; the partially wrapped (branch 1) minimum is
    the first one."""
    z = cache["z"]
    key = "s%02d" % k_idx
    s = float(cache[key + "_sigma"])
    E = -(w - 4.0) * z + s * z ** 2 + cache[key + "_F"]
    mins = []
    n = len(z)
    for i in range(n):
        l = E[i - 1] if i > 0 else np.inf
        r = E[i + 1] if i < n - 1 else np.inf
        if E[i] <= l and E[i] <= r:
            # refine parabolically
            if 0 < i < n - 1:
                den = E[i - 1] - 2 * E[i] + E[i + 1]
                dz = 0.5 * (E[i - 1] - E[i + 1]) / den * (z[1] - z[0]) if den > 0 else 0.0
                mins.append((z[i] + dz, E[i]))
            else:
                mins.append((z[i], E[i]))
    return mins, E


def partial_branch(flat, full, cache, zmax_partial=1.85):
    """For all rows in neither the enveloped state (z_cov < zmax_partial): E~_sim, z_sim against the theoretical
    partially wrapped minimum (first local minimum of E~(z), z_th_p, E_th_p) at the same (w~, sigma~), and against
    E~_th(z_sim).  The free state (z = 0, E~ = 0) is returned as the first minimum for w~ < 4."""
    if cache is None:
        return []
    out = []
    for nm, ds in (("flat", flat), ("full", full)):
        w, s, E, z, zad = deserno_vars(ds)
        for i in range(len(w)):
            k = int(round(ds["KA"][i] * 15))
            if k < 1 or z[i] >= zmax_partial:
                continue
            mins, Eth = theory_minima(cache, k, w[i])
            zz = cache["z"]
            zp, Ep = mins[0]
            Ez = np.interp(z[i], zz, Eth)
            out.append(dict(file=nm, KA=ds["KA"][i], KI=ds["KI"][i], w=w[i], sigma=s[i], z=z[i], zad=zad[i], E=E[i],
                            z_th=zp, E_th=Ep, E_th_at_zsim=Ez, E_th_env=-2 * (w[i] - 4) + 4 * s[i]))
    return out



def offset_models(rows, key, sigma_min=0.5):
    """Compare simple descriptions of the offset of one line (key 'E', 'S1' or 'S2') for sigma~ >= sigma_min, using
    only well defined points (E: interpolated; S1/S2: sharp jump).  Models for w_sim given w_th and sigma~:
      const : w_sim = w_th + c                  prop : w_sim = c w_th                  prop4 : w_sim = 4 + c (w_th - 4)
      lin_s : w_sim = w_th + c0 + c1 sigma~     flat : w_sim = c (independent of sigma~)
    returns {model: (parameters, rms residual)}; the rms of a uniform +-0.2 resolution error is 0.115."""
    sel = []
    for r in rows:
        if r["sigma_t"] < sigma_min or not np.isfinite(r[key + "_sim"]):
            continue
        if key == "E" and r["E_method"] != "interp":
            continue
        if key != "E" and not r[key + "_sharp"]:
            continue
        sel.append(r)
    ws = np.array([r[key + "_sim"] for r in sel]); wt = np.array([r["theory_" + key] for r in sel])
    sg = np.array([r["sigma_t"] for r in sel])
    out = {"n": len(sel)}
    c = np.mean(ws - wt); out["const"] = ((c,), np.sqrt(np.mean((ws - wt - c) ** 2)))
    c = np.sum(ws * wt) / np.sum(wt * wt); out["prop"] = ((c,), np.sqrt(np.mean((ws - c * wt) ** 2)))
    c = np.sum((ws - 4) * (wt - 4)) / np.sum((wt - 4) ** 2); out["prop4"] = ((c,), np.sqrt(np.mean((ws - 4 - c * (wt - 4)) ** 2)))
    A = np.vstack([np.ones_like(sg), sg]).T
    cc, *_ = np.linalg.lstsq(A, ws - wt, rcond=None); out["lin_s"] = (tuple(cc), np.sqrt(np.mean((ws - wt - A @ cc) ** 2)))
    c = np.mean(ws); out["flat"] = ((c,), np.sqrt(np.mean((ws - c) ** 2)))
    cc, *_ = np.linalg.lstsq(A, ws, rcond=None); out["abs_lin_s"] = (tuple(cc), np.sqrt(np.mean((ws - A @ cc) ** 2)))
    return out


def partial_offset_table(part, rows, edges=(2, 3, 4, 5, 6, 7, 8, 9, 10)):
    """Mean and spread of E~_sim - E~_th,partial (theory partially wrapped minimum, or the free state for w~ < 4) per w~ bin,
    for flat-start rows with sigma~ >= 0.5 whose theoretical partial minimum exists (w~ < S1_th)."""
    s1 = {round(r["KA"] * 15): r["theory_S1"] for r in rows}
    sel = [p for p in part if p["file"] == "flat" and p["sigma"] >= 0.5 and p["z_th"] < 1.9
           and p["w"] < s1[int(round(p["KA"] * 15))]]
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        d = np.array([p["E"] - p["E_th"] for p in sel if lo <= p["w"] < hi])
        zs = np.array([p["z"] for p in sel if lo <= p["w"] < hi])
        zt = np.array([p["z_th"] for p in sel if lo <= p["w"] < hi])
        zadd = np.array([p["zad"] for p in sel if lo <= p["w"] < hi])
        if len(d):
            out.append((lo, hi, len(d), d.mean(), d.std(), zs.mean(), zadd.mean(), zt.mean()))
    return out


def threshold_sensitivity(flat, full, thrs=(1.5, 1.75, 1.9)):
    """S1_sim / S2_sim (bracket mid, w~) for several z thresholds: {thr: (sigma, S1, S2)}."""
    grid = build_grid(flat, full)
    out = {}
    for t in thrs:
        L = []
        for k, g in grid.items():
            a, KB = g["a"], g["KB"]
            s1 = lowest_persistent(g["KI"], 2 * g["CoverageUnion_f"], t)
            s2 = lowest_persistent(g["KI"], 2 * g["CoverageUnion_u"], t)
            L.append((2 * float(np.nanmean(g["KA_f"])) * a ** 2 / KB, 4 * s1["mid"] * a ** 2 / KB, 4 * s2["mid"] * a ** 2 / KB))
        out[t] = np.array(L)
    return out


# =====================================================================================================================
# 6. Command line: tables and figures
# =====================================================================================================================
def _fit(x, y, yerr=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return None
    A = np.vstack([np.ones_like(x), x]).T
    c, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c
    s2 = (r @ r) / max(len(x) - 2, 1)
    cov = s2 * np.linalg.inv(A.T @ A)
    return dict(c0=c[0], c1=c[1], e0=np.sqrt(cov[0, 0]), e1=np.sqrt(cov[1, 1]), rms=np.sqrt(np.mean(r ** 2)), n=len(x))


def make_figures(rows, fig_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import deserno_theory as dt
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 120, "savefig.dpi": 200})
    cE, c1, c2 = "#222222", "#d95f02", "#1b7f8c"
    sig = np.array([r["sigma_t"] for r in rows])
    d = dt.fig2_line_data((0.0, 2.1), n=300)

    # ---------------- Figure 1: phase diagram in Deserno's axes
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.axvline(4, color="0.8", lw=1, ls=":", zorder=0)
    ax.plot(d["E"][0], d["E"][1], color=cE, lw=2.2, label=r"theory E", zorder=2)
    ax.plot(d["S1"][0], d["S1"][1], color=c1, lw=1.5, ls=(0, (4, 2)), label=r"theory S$_1$", zorder=2)
    ax.plot(d["S2"][0], d["S2"][1], color=c2, lw=1.5, ls=(0, (4, 2)), label=r"theory S$_2$", zorder=2)

    def pts(key, color, marker, label, errkey):
        x = np.array([r[key + "_sim"] for r in rows], float)
        sharp = np.array([r.get(key + "_sharp", r["E_method"] == "interp") for r in rows], bool)
        xe = np.array([r[errkey] for r in rows], float)
        stat = np.array([r.get(key + "_status", "ok") for r in rows])
        okm = np.isfinite(x)
        ax.errorbar(x[okm & sharp], sig[okm & sharp], xerr=xe[okm & sharp], fmt=marker, color=color, mfc=color, ms=6.5,
                    elinewidth=1.0, capsize=2, label=label, zorder=3)
        ax.errorbar(x[okm & ~sharp], sig[okm & ~sharp], xerr=xe[okm & ~sharp], fmt=marker, color=color, mfc="white", ms=6.5,
                    elinewidth=1.0, capsize=2, zorder=3)
        if key == "S1":                                   # flat start never jumps within the KI range: lower bound
            for r in rows:
                if r["S1_status"] == "high":
                    ax.annotate("", xy=(r["S1_lo"] + 2.2, r["sigma_t"]), xytext=(r["S1_lo"], r["sigma_t"]),
                                arrowprops=dict(arrowstyle="->", color=color, lw=1.2))
                    ax.plot([r["S1_lo"]], [r["sigma_t"]], marker=marker, color=color, mfc="white", ms=6.5, zorder=3)

    pts("S1", c1, "^", r"sim S$_1$ (flat start jumps)", "S1_err")
    pts("S2", c2, "v", r"sim S$_2$ (full start unwraps)", "S2_err")
    pts("E", cE, "s", r"sim E ($E_\mathrm{flat}=E_\mathrm{full}$)", "E_err")
    ax.set_xlim(2.5, 14.2)
    ax.set_ylim(0, 2.08)
    ax.set_xlabel(r"$\tilde w = 4K_I a^2/K_B$")
    ax.set_ylabel(r"$\tilde\sigma = 2K_A a^2/K_B$")
    ax.set_title("Planar bead: simulation (old data, shell 0.25) vs. Deserno theory", fontsize=11)
    ax.legend(loc="lower right", fontsize=9, frameon=False)
    fig.text(0.01, 0.005, r"horizontal bars: S$_{1,2}$ $\pm$ half a KI step (0.2 in $\tilde w$); E: interpolation error; open "
             r"markers: soft crossover, no sharp jump; arrows: no jump up to $\tilde w$=12", fontsize=7, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(fig_dir / "baseline_vs_deserno.png")
    plt.close(fig)

    # ---------------- Figure 2: offsets
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, xk in zip(axs, ("sigma_t", "theory")):
        for key, color, marker, lab in (("E", cE, "s", "E"), ("S1", c1, "^", r"S$_1$"), ("S2", c2, "v", r"S$_2$")):
            x = np.array([r["sigma_t"] if xk == "sigma_t" else r["theory_" + key] for r in rows], float)
            y = np.array([r["d" + key] for r in rows], float)
            e = np.array([r[key + "_err"] for r in rows], float)
            if key in ("S1", "S2"):
                sh = np.array([r[key + "_sharp"] for r in rows], bool)
            else:
                sh = np.array([r["E_method"] == "interp" for r in rows], bool)
            m = np.isfinite(y)
            ax.errorbar(x[m & sh], y[m & sh], yerr=e[m & sh], fmt=marker + "-", color=color, ms=6, lw=0.8, capsize=2, label=lab)
            ax.errorbar(x[m & ~sh], y[m & ~sh], yerr=e[m & ~sh], fmt=marker, color=color, mfc="white", ms=6, capsize=2)
            if key == "S1":
                for r in rows:
                    if r["S1_status"] == "high":
                        xx = r["sigma_t"] if xk == "sigma_t" else r["theory_S1"]
                        yy = r["S1_lo"] - r["theory_S1"]
                        ax.annotate("", xy=(xx, yy + 1.8), xytext=(xx, yy), arrowprops=dict(arrowstyle="->", color=color, lw=1.2))
                        ax.plot([xx], [yy], marker="^", color=color, mfc="white", ms=6)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_ylabel(r"$\Delta\tilde w=\tilde w_\mathrm{sim}-\tilde w_\mathrm{theory}$")
        ax.set_xlabel(r"$\tilde\sigma$" if xk == "sigma_t" else r"$\tilde w_\mathrm{theory}$ of the same line")
    axs[0].legend(frameon=False)
    axs[0].set_title("offset vs surface tension", fontsize=11)
    axs[1].set_title("offset vs theoretical position of the line", fontsize=11)
    fig.tight_layout()
    fig.savefig(fig_dir / "baseline_offsets.png")
    plt.close(fig)


def make_energy_figure(env, part, fig_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 120, "savefig.dpi": 200})
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    w = np.array([e["w"] for e in env]); s = np.array([e["sigma"] for e in env])
    sc = axs[0].scatter(w, [e["resid"] for e in env], c=s, cmap="viridis", s=14)
    axs[0].set_xlabel(r"$\tilde w$"); axs[0].set_ylabel(r"$\tilde E_\mathrm{sim}-[-2(\tilde w-4)+4\tilde\sigma]$")
    axs[0].set_title("enveloped state: energy residual", fontsize=11)
    plt.colorbar(sc, ax=axs[0], label=r"$\tilde\sigma$")
    for nm, key, c in (("adhesion", "resid_adh", "#d95f02"), ("bending", "resid_bend", "#1b7f8c"), ("tension", "resid_tens", "#7570b3")):
        axs[1].scatter(w, [e[key] for e in env], s=10, color=c, label=nm)
    axs[1].axhline(0, color="0.5", lw=0.8)
    axs[1].set_xlabel(r"$\tilde w$"); axs[1].set_ylabel(r"term minus ideal (units $\pi\kappa$)")
    axs[1].set_title("which term carries the residual", fontsize=11); axs[1].legend(frameon=False)
    part = [p for p in part if p["z_th"] < 1.9 and p["file"] == "flat"]   # theory partial minimum exists
    if part:
        w2 = np.array([p["w"] for p in part]); s2 = np.array([p["sigma"] for p in part])
        sc2 = axs[2].scatter(w2, [p["E"] - p["E_th"] for p in part], c=s2, cmap="viridis", s=12)
        axs[2].axhline(0, color="0.5", lw=0.8)
        axs[2].set_xlabel(r"$\tilde w$"); axs[2].set_ylabel(r"$\tilde E_\mathrm{sim}-\tilde E_\mathrm{th}$ (partial minimum)")
        axs[2].set_title("flat start vs. theory partial minimum", fontsize=11)
        plt.colorbar(sc2, ax=axs[2], label=r"$\tilde\sigma$")
    fig.tight_layout()
    fig.savefig(fig_dir / "baseline_energy.png")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--flat", default=str(FLAT_FILE))
    ap.add_argument("--full", default=str(FULL_FILE))
    ap.add_argument("--data-dir", default=str(STUDY / "data"))
    ap.add_argument("--fig-dir", default=str(STUDY / "figures"))
    a = ap.parse_args(argv)
    flat, full = load_coverage(a.flat), load_coverage(a.full)
    for nm, d in (("flat", flat), ("full", full)):
        print(nm, d["info"])
    print("reference area", check_reference_area(flat, full))
    rows = extract_all(flat, full)
    data_dir, fig_dir = Path(a.data_dir), Path(a.fig_dir)
    write_csv(rows, data_dir / "baseline_lines.csv")
    grid = build_grid(flat, full)
    print("same-branch scatter", same_branch_scatter(grid))
    for which in ("f", "u"):
        r, rel, info = hellmann_feynman(grid, which)
        print("Hellmann-Feynman", which, "n=%d" % len(r), "rms resid %.4f  median |rel| %.4f  p95 |rel| %.4f" %
              (np.sqrt(np.mean(r ** 2)), np.median(np.abs(rel)), np.percentile(np.abs(rel), 95)))
    make_figures(rows, fig_dir)
    for key in ("E", "S1", "S2"):
        om = offset_models(rows, key)
        print("offset models", key, {k: (tuple(round(float(x), 3) for x in v[0]), round(float(v[1]), 3)) if k != "n" else v for k, v in om.items()})
    env = enveloped_rows(flat, full)
    cache = load_theory_cache()
    part = partial_branch(flat, full, cache)
    make_energy_figure(env, part, fig_dir)
    if part:
        with open(data_dir / "baseline_partial_branch.csv", "w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(list(part[0].keys()))
            for p in part:
                wr.writerow([("%.6g" % v) if isinstance(v, float) else v for v in p.values()])
    return rows, env, part, grid


if __name__ == "__main__":
    main()
