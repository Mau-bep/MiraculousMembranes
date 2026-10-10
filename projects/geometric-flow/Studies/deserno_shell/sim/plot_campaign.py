#!/usr/bin/env python
"""Figures of the shell-steepness campaign: simulated energies against the continuum soft-shell branches, and the lines vs the effective shell width.

    micromamba run -n mir_mem python plot_campaign.py   (reads ../data/campaign_s025.csv, campaign_lines.csv, campaign_theory_cache.json)
"""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FIG = os.path.join(HERE, "..", "figures")


def rd(f):
    return list(csv.DictReader(open(os.path.join(DATA, f))))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def seff_factor(p):
    y = np.linspace(0, 1, 20001)
    W = lambda x, pp: (0.5 * (1 + np.cos(np.pi * x))) ** pp
    return (np.trapezoid(W(y ** 2, p), y) / np.trapezoid(W(y ** 2, 1.0), y)) ** 2


def main():
    sims = rd("campaign_s025.csv")
    lines = rd("campaign_lines.csv")
    cache = json.load(open(os.path.join(DATA, "campaign_theory_cache.json")))
    tags = sorted({r["tag"] for r in sims}, key=lambda t: (-float(next(r["p"] for r in sims if r["tag"] == t)), t))
    fig, axs = plt.subplots(1, len(tags), figsize=(3.6 * len(tags), 3.6), sharey=False)
    axs = np.atleast_1d(axs)
    for ax, tag in zip(axs, tags):
        rr = [r for r in sims if r["tag"] == tag]
        s, p, sg = num(rr[0]["s"]), num(rr[0]["p"]), num(rr[0]["sigma_t"])
        key = "%.4f_%.4f_%.3f" % (sg, s, p)
        c = cache[key]
        G, e, wt = map(np.array, (c["G"], c["e"], c["wt"]))
        E = e - wt * G
        ax.plot(wt, E, "-", color="0.55", lw=1.2, label="continuum branches")
        for init, mk, col in (("flat", "o", "C3"), ("full", "s", "C0")):
            q = [r for r in rr if r["init"] == init]
            ax.plot([num(r["w_t"]) for r in q], [num(r["E_sim"]) for r in q], mk, color=col, ms=5, label=init + " start")
        wlo = min(num(r["w_t"]) for r in rr) - 0.3
        whi = max(num(r["w_t"]) for r in rr) + 0.3
        ax.set_xlim(wlo, whi)
        m = (wt > wlo) & (wt < whi)
        ax.set_ylim(E[m].min() - 0.3, E[m].max() + 0.3)
        ax.set_title("p = %g, s = %.2f" % (p, s))
        ax.set_xlabel(r"$\tilde w$")
    axs[0].set_ylabel(r"$E/(\pi\kappa)$, $\tilde\sigma$ = %.2f" % sg)
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "campaign_energies.png"), dpi=140)

    fig, ax = plt.subplots(1, 1, figsize=(6.2, 4.2))
    for L in lines:
        p, s = num(L["p"]), num(L["s"])
        se = s * seff_factor(p)
        for key, col, name in (("S1", "C3", "S1"), ("E", "k", "E"), ("S2", "C0", "S2")):
            ax.plot(se, num(L[key + "_cont"]), "_", color=col, ms=14, mew=2)
        ax.plot(se, num(L["E_sim"]), "o", color="k", mfc="none", ms=7)
        for key, col in (("S1", "C3"), ("S2", "C0")):
            lo, hi = num(L[key + "_sim_lo"]), num(L[key + "_sim_hi"])
            if np.isfinite(lo) and np.isfinite(hi):
                ax.plot([se, se], [lo, hi], "-", color=col, lw=3, alpha=0.5)
    sg = num(lines[0]["sigma_t"])
    for key, col in (("S1", "C3"), ("E", "k"), ("S2", "C0")):
        ax.axhline(num(lines[0][key + "_deserno"]), color=col, ls=":", lw=1)
    ax.plot([], [], "_", color="0.3", ms=14, mew=2, label="continuum soft shell (dash)")
    ax.plot([], [], "o", color="k", mfc="none", label="simulated E (flat/full crossing)")
    ax.plot([], [], "-", color="0.5", lw=3, alpha=0.5, label="simulated S1 (red), S2 (blue): bracket")
    ax.plot([], [], ":", color="0.3", label="Deserno (zero range)")
    ax.set_xlabel(r"effective shell width $s_{\rm eff} = s\,(I_p/I_1)^2$")
    ax.set_ylabel(r"$\tilde w$")
    ax.set_title(r"$\tilde\sigma$ = %.2f: S1 (red), E (black), S2 (blue)" % sg)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "campaign_lines_vs_seff.png"), dpi=140)


if __name__ == "__main__":
    main()
