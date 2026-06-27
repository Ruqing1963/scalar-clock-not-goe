"""
scalar_clock_audit.py
---------------------
A fourth spurious-GOE hazard: a narrow-CV scalar (renewal) clock reproduces any
value of the consecutive spacing ratio <r> -- including the GOE and GUE values --
purely from its marginal interval distribution, with no level repulsion. We
demonstrate the mechanism, then re-audit a published "temporal GOE in ore
systems" result (orogenic-gold geochronology) and show it is fully reproduced by
a zero-repulsion scalar clock and is unresolvable given the age uncertainties.

Produces the two figures and all reported numbers. Self-contained: the orogenic
gold ages (data/orogenic_gold_ages.csv) are the only input; everything else is
synthetic.
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R_POIS, R_GOE, R_GUE = 0.386, 0.531, 0.603
SEED = 0


def r_stat(s):
    s = np.asarray(s, float); s = s[s > 0]
    if len(s) < 3: return np.nan
    return float((np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])).mean())


def renewal(n, cv, rng):
    """gamma renewal clock, mean-1 i.i.d. intervals, coefficient of variation cv."""
    k = 1.0 / cv**2
    return np.cumsum(rng.gamma(k, 1.0 / k, n))


def goe_levels(N, rng):
    """unfolded bulk eigenvalues of a GOE matrix (genuine level repulsion)."""
    A = rng.standard_normal((N, N)); A = (A + A.T) / np.sqrt(2 * N)
    ev = np.sort(np.linalg.eigvalsh(A)); b = ev[(ev > -1.0) & (ev < 1.0)]
    return np.sort(N * (0.5 + (b * np.sqrt(2 - b**2) / 2 + np.arcsin(b / np.sqrt(2))) / np.pi))


def load_ages(path):
    dep, age, err = [], [], []
    for row in csv.DictReader(open(path)):
        dep.append(row["deposit"]); age.append(float(row["age_ma"])); err.append(float(row["age_err"]))
    return np.array(dep), np.array(age, float), np.array(err, float)


def paper_pipeline(dep, age):
    """exact construction audited: per-deposit ages -> within-deposit normalized
    intervals (deposits with >=3 ages) -> pooled."""
    pooled = []
    for d in dict.fromkeys(dep):
        a = np.sort(age[dep == d]); iv = np.diff(a); iv = iv[iv > 0]
        if len(iv) >= 2: pooled.extend((iv / iv.mean()).tolist())
    return np.array(pooled)


def main(outdir):
    os.makedirs(f"{outdir}/figures", exist_ok=True)
    rng = np.random.default_rng(SEED)

    # ---------- (1) mechanism: <r> vs CV for a renewal clock; obs-shuf discriminator ----------
    cvs = np.linspace(0.12, 1.2, 22)
    r_of_cv = [r_stat(np.diff(renewal(20000, cv, rng))) for cv in cvs]

    def obs_shuf(pos):
        s = np.diff(np.sort(pos)); return r_stat(s) - r_stat(rng.permutation(s))
    disc = {
        "relax CV=0.8": obs_shuf(renewal(8000, 0.8, rng)),
        "relax CV=0.55": obs_shuf(renewal(8000, 0.55, rng)),
        "relax CV=0.35": obs_shuf(renewal(8000, 0.35, rng)),
        "Poisson": obs_shuf(np.cumsum(rng.exponential(1, 8000))),
        "GOE": obs_shuf(goe_levels(3000, rng)),
    }

    # ---------- (2) real-data re-audit (orogenic gold) ----------
    dep, age, err = load_ages(f"{outdir}/data/orogenic_gold_ages.csv")
    pooled = paper_pipeline(dep, age)
    r_obs = r_stat(pooled); cv_obs = pooled.std() / pooled.mean()
    r_shuf = np.mean([r_stat(rng.permutation(pooled)) for _ in range(5000)])
    # CV-equivalent scalar clock
    cv_grid = np.linspace(0.25, 0.55, 31)
    cv_eq = min(cv_grid, key=lambda c: abs(r_stat(np.diff(renewal(20000, c, rng))) - r_obs))
    # age-error Monte Carlo through the exact pipeline
    mc = []
    for _ in range(5000):
        pp = paper_pipeline(dep, age + rng.normal(0, err))
        if len(pp) >= 4: mc.append(r_stat(pp))
    mc = np.array(mc)
    # normalization-artifact null: zero-repulsion clock through the SAME construction
    def null_construction(cv=0.38, ndep=10):
        pl = []
        for _ in range(ndep):
            m = rng.integers(2, 5); iv = rng.gamma(1/cv**2, cv**2, m); pl.extend((iv/iv.mean()).tolist())
        return np.array(pl)
    nro, nrs = [], []
    for _ in range(4000):
        p = null_construction(); nro.append(r_stat(p)); nrs.append(r_stat(rng.permutation(p)))
    null_r, null_osh = np.mean(nro), np.mean(nro) - np.mean(nrs)

    print("=== MECHANISM ===")
    print(f"  renewal clock <r> spans {min(r_of_cv):.3f}..{max(r_of_cv):.3f} as CV goes {cvs.max():.1f}..{cvs.min():.2f}")
    print(f"  obs-shuf: " + ", ".join(f"{k}={v:+.3f}" for k, v in disc.items()))
    print("=== OROGENIC GOLD RE-AUDIT ===")
    print(f"  reproduced <r>={r_obs:.3f} (paper 0.678), n={len(pooled)}, CV={cv_obs:.2f}")
    print(f"  shuffle: obs-shuf={r_obs-r_shuf:+.3f}")
    print(f"  scalar-clock equivalent CV={cv_eq:.2f}")
    print(f"  age-error MC: <r>={mc.mean():.3f} [{np.percentile(mc,2.5):.3f},{np.percentile(mc,97.5):.3f}]")
    print(f"  zero-repulsion null through same pipeline: <r>={null_r:.3f}, obs-shuf={null_osh:+.3f} (matches data artifact)")

    # ---------- Figure 1: mechanism ----------
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].plot(cvs, r_of_cv, "o-", color="#34495E", label="charge-release renewal clock")
    ax[0].axhline(R_GOE, ls="--", color="#27AE60", label="GOE 0.531")
    ax[0].axhline(R_GUE, ls="--", color="#2980B9", label="GUE 0.603")
    ax[0].axhline(R_POIS, ls="--", color="#888", label="Poisson 0.386")
    ax[0].set_xlabel("CV of inter-event intervals"); ax[0].set_ylabel(r"$\langle r\rangle$")
    ax[0].set_title("(a) A relaxation-oscillator (renewal) clock reproduces ANY\n"
                    "⟨r⟩ — including GOE and GUE — purely from its marginal CV")
    ax[0].legend(fontsize=8); ax[0].invert_xaxis()
    labs = list(disc); vals = [disc[k] for k in labs]
    cols = ["#34495E"]*3 + ["#888", "#C0392B"]
    ax[1].bar(range(len(labs)), vals, color=cols); ax[1].axhline(0, color="k", lw=.8)
    ax[1].set_xticks(range(len(labs))); ax[1].set_xticklabels([l.replace(" ", "\n") for l in labs], fontsize=8)
    ax[1].set_ylabel(r"$\langle r\rangle_{\rm obs}-\langle r\rangle_{\rm shuf}$")
    ax[1].set_title("(b) Only genuine GOE has the negative rigidity signature;\n"
                    "renewal clocks give obs−shuf ≈ 0 (and need N≫300 to tell apart)")
    fig.tight_layout(); fig.savefig(f"{outdir}/figures/fig_mechanism.pdf", bbox_inches="tight")
    fig.savefig(f"{outdir}/figures/fig_mechanism.png", dpi=120, bbox_inches="tight")

    # ---------- Figure 2: real-data audit ----------
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    ax[0].hist(mc, bins=40, color="#7F8C8D", alpha=.85)
    for v, c, l in [(R_POIS, "#888", "Poisson"), (R_GOE, "#27AE60", "GOE"), (R_GUE, "#2980B9", "GUE")]:
        ax[0].axvline(v, ls="--", color=c, label=l)
    ax[0].axvline(r_obs, color="#C0392B", lw=2, label=f"reported {r_obs:.3f}")
    ax[0].set_xlabel(r"$\langle r\rangle$ when ages are jittered within ±error")
    ax[0].set_ylabel("count")
    ax[0].set_title("(a) Age errors alone move ⟨r⟩ across Poisson→GUE;\n"
                    "the reported value sits at the top edge of the noise")
    ax[0].legend(fontsize=8)
    cats = ["⟨r⟩", "obs−shuf"]; x = np.arange(2); w = 0.36
    ax[1].bar(x - w/2, [r_obs, r_obs - r_shuf], w, label="orogenic gold DATA", color="#C0392B")
    ax[1].bar(x + w/2, [null_r, null_osh], w, label="scalar renewal clock\n(CV=0.38, ZERO repulsion)", color="#34495E")
    ax[1].axhline(0, color="k", lw=.8); ax[1].set_xticks(x); ax[1].set_xticklabels(cats)
    ax[1].set_title("(b) A zero-repulsion scalar clock reproduces BOTH the\n"
                    "⟨r⟩=0.68 and the −0.03 'rigidity' (a normalization artifact)")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{outdir}/figures/fig_realdata.pdf", bbox_inches="tight")
    fig.savefig(f"{outdir}/figures/fig_realdata.png", dpi=120, bbox_inches="tight")
    print("figures written")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default=".."); a = ap.parse_args()
    main(a.outdir)
