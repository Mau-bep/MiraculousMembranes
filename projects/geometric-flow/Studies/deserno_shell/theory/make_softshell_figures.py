#!/usr/bin/env python
"""Figures of the continuum soft-shell study (workstream T1).  Reads ../data/softshell_lines.csv, softshell_vs_sim_*.csv,
softshell_branches.npz and writes ../figures/softshell_*.png

    micromamba run -n mir_mem python make_softshell_figures.py
"""
import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FIG = os.path.join(HERE, "..", "figures")
COL = dict(S1="tab:red", E="k", S2="tab:blue")
MARK = dict(S1="o", E="s", S2="^")
SIG = [0.13, 0.27, 0.53, 1.0]


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def rows_main():
    rs = read(os.path.join(DATA, "softshell_lines.csv"))
    for r in rs:
        for k in ("sigma_t", "s", "p", "R", "S1", "E", "S2", "G_S1", "G_S2", "G_E_partial", "G_E_env",
                  "deserno_S1", "deserno_E", "deserno_S2"):
            r[k] = fl(r.get(k))
    return rs


def fig_vs_s(rs):
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.8), sharey=False)
    for ax, sg in zip(axs, SIG):
        sel = [r for r in rs if abs(r["sigma_t"] - sg) < 1e-6 and r["p"] == 1 and r["R"] == 10.0 and r["set"] in ("main", "s0")]
        d = {}
        for r in sel:
            d[r["s"]] = r
        ss = sorted(d)
        for key in ("S1", "E", "S2"):
            y = np.array([d[s][key] for s in ss])
            ax.plot(ss, y, MARK[key] + "-", color=COL[key], label=key + " soft shell", ms=4)
            ax.axhline(d[ss[0]]["deserno_" + key], color=COL[key], ls="--", lw=1)
        nb = [s for s in ss if not np.isfinite(d[s]["S1"])]
        for s in nb:
            ax.axvline(s, color="0.8", lw=6, zorder=0)
        ax.set_xscale("log")
        ax.set_xlabel("shell half width s")
        ax.set_title(r"$\tilde\sigma$ = %.2f" % sg)
        if ax is axs[0]:
            ax.set_ylabel(r"$\tilde w$ of the line")
            ax.legend(fontsize=7)
        ax.text(0.02, 0.02, "grey band: no barrier\n(smooth crossover)" if nb else "", transform=ax.transAxes, fontsize=6, color="0.4")
    fig.suptitle("continuum soft shell (p = 1, R = 10) vs shell width; dashed = Deserno (zero range, R = inf)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "softshell_lines_vs_s.png"), dpi=130)
    plt.close(fig)


def fig_vs_p(rs):
    ps = [1, 2, 4, 8]
    svals = sorted({r["s"] for r in rs if r["set"] == "main"}, reverse=True)
    cmap = plt.get_cmap("viridis")
    fig, axs = plt.subplots(3, 4, figsize=(15, 9), sharex=True)
    for j, sg in enumerate(SIG):
        for i, key in enumerate(("S1", "E", "S2")):
            ax = axs[i, j]
            for k, s in enumerate(svals):
                y = []
                for p in ps:
                    m = [r for r in rs if r["set"] == "main" and abs(r["sigma_t"] - sg) < 1e-6 and r["s"] == s and r["p"] == p]
                    y.append(m[0][key] if m else np.nan)
                ax.plot(ps, y, "o-", color=cmap(k / max(1, len(svals) - 1)), ms=3.5, label="s = %g" % s)
            dd = [r for r in rs if abs(r["sigma_t"] - sg) < 1e-6]
            if dd:
                ax.axhline(dd[0]["deserno_" + key], color="k", ls="--", lw=1.2, label="Deserno")
            ax.set_xscale("log", base=2)
            ax.set_xticks(ps)
            ax.set_xticklabels([str(p) for p in ps])
            if i == 0:
                ax.set_title(r"$\tilde\sigma$ = %.2f" % sg)
            if j == 0:
                ax.set_ylabel(r"$\tilde w_{%s}$" % key)
            if i == 2:
                ax.set_xlabel("shell power p  (w = cos$^{2p}$($\\pi x/2$))")
    axs[0, 0].legend(fontsize=6, ncol=2)
    fig.suptitle("lines of the continuum soft shell vs shell steepness p, colour = shell half width s")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "softshell_lines_vs_p.png"), dpi=130)
    plt.close(fig)


def fig_s0_scaling(rs):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    for ax, key in zip(axs, ("S1", "E", "S2")):
        for sg, c in zip(SIG, ("tab:green", "tab:orange", "tab:purple", "tab:brown")):
            sel = {r["s"]: r for r in rs if abs(r["sigma_t"] - sg) < 1e-6 and r["p"] == 1 and r["R"] == 10.0 and r["set"] in ("main", "s0")}
            ss = np.array(sorted(sel))
            sh = np.array([sel[s][key] - sel[s]["deserno_" + key] for s in ss])
            ok = np.isfinite(sh) & (sh > 0)
            ax.plot(ss[ok], sh[ok], "o-", color=c, ms=4, label=r"$\tilde\sigma$=%.2f" % sg)
        x = np.array([0.005, 0.25])
        ax.plot(x, 4 * np.sqrt(x), "k:", lw=1, label=r"$\propto s^{1/2}$")
        ax.plot(x, 4 * x, "k-.", lw=1, label=r"$\propto s$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("s")
        ax.set_ylabel("%s(soft) - %s(Deserno)" % (key, key))
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "softshell_s0_convergence.png"), dpi=130)
    plt.close(fig)


def fig_collapse(rs):
    from scipy.integrate import quad
    W = lambda x, p: (0.5 * (1 + np.cos(np.pi * x))) ** p if abs(x) < 1 else 0.0
    I = {p: quad(lambda y: W(y * y, p), 0, 1)[0] for p in (1, 2, 4, 8)}
    fig, axs = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    mk = {1: "o", 2: "s", 4: "^", 8: "v"}
    for i, sg in enumerate((0.53, 1.0)):
        for j, key in enumerate(("S1", "E", "S2")):
            ax = axs[i, j]
            for p in (1, 2, 4, 8):
                pts = [(r["s"] * (I[p] / I[1]) ** 2, r[key] - r["deserno_" + key]) for r in rs
                       if r["set"] in ("main", "s0") and abs(r["sigma_t"] - sg) < 1e-6 and r["p"] == p and r["R"] == 10.0 and np.isfinite(r[key])]
                pts.sort()
                if pts:
                    ax.plot(*zip(*pts), mk[p] + "-", ms=4, lw=0.8, label="p = %d" % p)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(r"$\tilde\sigma$=%.2f  %s - Deserno" % (sg, key))
            if i == 1:
                ax.set_xlabel(r"$s_{eff} = s\,(I_p/I_1)^2$,  $I_p=\int_0^1 W(y^2)dy$")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("all shell powers collapse onto the p = 1 curve when s is replaced by s_eff")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "softshell_collapse_seff.png"), dpi=130)
    plt.close(fig)


def fig_vs_sim():
    p = os.path.join(DATA, "softshell_vs_sim_lines.csv")
    if not os.path.exists(p):
        return
    rs = read(p)
    sg = np.array([fl(r["sigma_t"]) for r in rs])
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for key, k2 in (("S1", "S1"), ("E", "E"), ("S2", "S2")):
        ax[0].plot(sg, [fl(r[key + "_cont"]) for r in rs], "-", color=COL[key], lw=2, label=key + " continuum soft shell (s = 0.25)")
        ax[0].plot(sg, [fl(r[key + "_deserno"]) for r in rs], "--", color=COL[key], lw=1, label=key + " Deserno")
        y = np.array([fl(r[key + "_sim"]) for r in rs])
        ax[0].plot(sg, y, MARK[key], color=COL[key], mfc="none", ms=7, label=key + " simulation")
    ax[0].set_xlabel(r"$\tilde\sigma$")
    ax[0].set_ylabel(r"$\tilde w$")
    ax[0].legend(fontsize=6, ncol=1)
    ax[0].set_title("lines: continuum soft shell vs simulation vs Deserno")
    for key in ("S1", "E", "S2"):
        ax[1].plot(sg, np.array([fl(r[key + "_sim"]) for r in rs]) - np.array([fl(r[key + "_cont"]) for r in rs]), MARK[key] + "-", color=COL[key], label="sim - continuum " + key)
        ax[1].plot(sg, np.array([fl(r[key + "_sim"]) for r in rs]) - np.array([fl(r[key + "_deserno"]) for r in rs]), MARK[key] + ":", color=COL[key], mfc="none", label="sim - Deserno " + key)
    ax[1].axhline(0, color="0.5", lw=0.8)
    ax[1].set_xlabel(r"$\tilde\sigma$")
    ax[1].set_ylabel(r"offset in $\tilde w$")
    ax[1].legend(fontsize=7)
    ax[1].set_title("offset left by the continuum shell (solid) vs by Deserno (dotted)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "softshell_vs_sim_lines.png"), dpi=130)
    plt.close(fig)


def fig_partial():
    p = os.path.join(DATA, "softshell_vs_sim_partial.csv")
    if not os.path.exists(p):
        return
    rs = read(p)
    sg = np.array([fl(r["sigma_t"]) for r in rs])
    w = np.array([fl(r["w"]) for r in rs])
    dE = np.array([fl(r["dE"]) for r in rs])
    dG = np.array([fl(r["dG"]) for r in rs])
    Es = np.array([fl(r["E_sim"]) for r in rs])
    Ec = np.array([fl(r["E_cont"]) for r in rs])
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    sc = ax[0].scatter(w, Es, c=sg, s=10, cmap="viridis", label="simulation")
    ax[0].scatter(w, Ec, c=sg, s=10, marker="x", cmap="viridis", label="continuum soft shell")
    ax[0].set_xlabel(r"$\tilde w$")
    ax[0].set_ylabel(r"$E/(\pi\kappa)$ of the partially wrapped state")
    ax[0].legend(fontsize=7)
    plt.colorbar(sc, ax=ax[0], label=r"$\tilde\sigma$")
    sc = ax[1].scatter(w, dE, c=sg, s=10, cmap="viridis")
    ax[1].axhline(0, color="0.5")
    ax[1].set_xlabel(r"$\tilde w$")
    ax[1].set_ylabel(r"$E_{sim} - E_{cont}$")
    plt.colorbar(sc, ax=ax[1], label=r"$\tilde\sigma$")
    sc = ax[2].scatter(w, dG, c=sg, s=10, cmap="viridis")
    ax[2].axhline(0, color="0.5")
    ax[2].set_xlabel(r"$\tilde w$")
    ax[2].set_ylabel(r"$z_{ad}^{sim} - G_{cont}$")
    plt.colorbar(sc, ax=ax[2], label=r"$\tilde\sigma$")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "softshell_vs_sim_partial.png"), dpi=130)
    plt.close(fig)


def fig_branches():
    p = os.path.join(DATA, "softshell_branches.npz")
    cache = os.path.join(DATA, "theory_cache.npz")
    if not (os.path.exists(p) and os.path.exists(cache)):
        return
    br = np.load(p)
    tc = np.load(cache)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for k, c in ((2, "tab:blue"), (4, "tab:orange"), (8, "tab:green")):
        sg = float(tc["s%02d_sigma" % k])
        z = tc["z"]
        wD = (1 - tc["s%02d_pd" % k]) ** 2
        cands = [f[2:] for f in br.files if f.startswith("G_") and abs(float(f[2:]) - sg) < 2e-3]
        if not cands:
            continue
        key = cands[0]
        ax[0].plot(z, wD, "--", color=c, label=r"Deserno $\tilde\sigma$=%.2f" % sg)
        ax[0].plot(br["G_" + key], br["wt_" + key], "-", color=c, label=r"soft shell s=0.25")
        eD = 4 * z + sg * z ** 2 + tc["s%02d_F" % k]
        ax[1].plot(z, eD, "--", color=c)
        ax[1].plot(br["G_" + key], br["e_" + key], "-", color=c)
    ax[0].set_xlim(0, 2.0)
    ax[0].set_ylim(0, 12)
    ax[0].set_xlabel("coverage z (Deserno) / G (soft shell)")
    ax[0].set_ylabel(r"$\tilde w(z) = de/dz$")
    ax[0].legend(fontsize=7)
    ax[1].set_xlabel("coverage z / G")
    ax[1].set_ylabel(r"$e = E_{el}/(\pi\kappa)$ (bending + tension)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "softshell_branch_w_of_G.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    rs = rows_main()
    fig_vs_s(rs)
    fig_vs_p(rs)
    fig_s0_scaling(rs)
    fig_collapse(rs)
    fig_vs_sim()
    fig_partial()
    fig_branches()
    print("figures written to", FIG)
