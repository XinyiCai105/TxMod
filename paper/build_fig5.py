
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import data, out
from figstyle import apply_figure_style, set_frame, panel_letter
import os, numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu
import statsmodels.api as sm
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_lib import load_recurrence, draw_recurrence

DEP_T, CGC_T, RBP_T = data("depletion_full_threshold_scan.tsv"), data("cgc_gene_level_table_v2.tsv"), data("rbp_binding_alterations_v2.tsv")

# ── panel B: composition-corrected mutation density
DEP = pd.read_csv(DEP_T, sep="\t")
# Every threshold at which the site definition is populated enough to estimate a
# ratio: n >= 100 sites. This keeps 0.1-0.8 and drops only 0.9, which has 3 sites.
dep = DEP[DEP.n_sites >= 100]

# ── gene sets: RBP-altered is nested inside m6A-altering; background is the complement
GL = pd.read_csv(CGC_T, sep="\t")
RB = pd.read_csv(RBP_T, sep="\t")
univ = set(GL.gene_name)
m6g = set(GL[GL.m6A == 1].gene_name)
rbg = set(RB.record_key.str.split("|").str[2].unique()) & univ
bgg = univ - m6g

# ── panel D: tumour/normal expression, one value per gene (largest |log2FC|)
V3 = pd.read_csv(data("TCGA_pan_cancer_DEG_table_v3.tsv"), sep="\t")
V3g = V3.sort_values("abs_log2FC").drop_duplicates("gene_name", keep="last")
EX = {k: V3g[V3g.gene_name.isin(S)].abs_log2FC.dropna().values
      for k, S in (("background", bgg), ("m6A", m6g), ("rbp", rbg))}
PE = {k: mannwhitneyu(EX[k], EX["background"], alternative="two-sided")[1] for k in ("m6A", "rbp")}

# ── panel C: Census enrichment, Fisher then logistic adjustment
a_ = int(((GL.m6A == 1) & (GL.is_CGC == 1)).sum()); b_ = int(((GL.m6A == 1) & (GL.is_CGC == 0)).sum())
c_ = int(((GL.m6A == 0) & (GL.is_CGC == 1)).sum()); d_ = int(((GL.m6A == 0) & (GL.is_CGC == 0)).sum())
or_raw, p_raw = fisher_exact([[a_, b_], [c_, d_]])
g0 = GL.dropna(subset=["total_utr_length", "n_mutations", "n_transcripts"])
X = sm.add_constant(pd.DataFrame({"m6A": g0.m6A.astype(float),
                                  "lu": np.log10(g0.total_utr_length + 1),
                                  "lm": np.log10(1 + g0.n_mutations),
                                  "nt": g0.n_transcripts.astype(float)}, index=g0.index))
fit = sm.Logit((g0.is_CGC == 1).astype(float), X).fit(disp=0)
or_adj = float(np.exp(fit.params["m6A"])); lo_, hi_ = np.exp(fit.conf_int().loc["m6A"])
p_adj = float(fit.pvalues["m6A"])

# ── panel A: recurrence across tumour samples
R = load_recurrence()

GRY="#9AA5B1"; PUR="#7B6AA0"; RED="#C1483A"; BLU="#1F5C8B"; TEA="#2E7D6E"; ORA="#C2622A"; INK="#1F2933"
rng = np.random.default_rng(0)
plt.close("all")
apply_figure_style(sizes=(9, 8, 7))
fig5 = plt.figure(figsize=(6.35, 5.80))
gs = fig5.add_gridspec(2, 2, hspace=.58, wspace=.66, left=.095, right=.975, top=.94, bottom=.085)
# Panel letters follow first citation in Results, which reaches the Census correction
# before recurrence, so the Census forest plot takes the top-left cell and recurrence
# the lower-left one. a1 still holds recurrence and a3 the Census panel.
a3, a2 = fig5.add_subplot(gs[0, 0]), fig5.add_subplot(gs[0, 1])
a1, a4 = fig5.add_subplot(gs[1, 0]), fig5.add_subplot(gs[1, 1])

draw_recurrence(a1, R, frame=set_frame)

a2.axhline(1, ls="--", lw=.8, color=GRY, zorder=1)
a2.plot(dep.threshold, dep.raw_center_flank_ratio, "-s", ms=3.4, lw=1.0, color=ORA,
        zorder=3, label="Uncorrected")
a2.errorbar(dep.threshold, dep.corrected_ror,
            yerr=[dep.corrected_ror - dep.ci_low, dep.ci_high - dep.corrected_ror],
            fmt="-o", ms=3.6, lw=1.0, elinewidth=.9, capsize=1.8, color=TEA,
            zorder=4, label="Composition-corrected")
for k_, (th, chi, n) in enumerate(zip(dep.threshold, dep.ci_high, dep.n_sites)):
    # Above the interval, alternating two heights so neighbouring counts cannot touch
    # at 7 pt; first and last aligned inward so the label stays inside the axes.
    ha_ = "left" if th == dep.threshold.min() else ("right" if th == dep.threshold.max() else "center")
    y_ = chi + (.014 if k_ % 2 == 0 else .034)
    if abs(y_ - 1.0) < .014: y_ = 1.006        # clear of the unity reference line
    lab_ = ("%.0fk" % (n/1000)) if n >= 10000 else (("%.1fk" % (n/1000)) if n >= 1000 else "%d" % n)
    a2.text(th, y_, lab_, ha=ha_,
            va="bottom", fontsize=7.0, color=GRY)

a2.set_xlabel("iM6A probability threshold"); a2.set_ylabel("Centre / flank density ratio")
a2.set_ylim(.76, 1.17); a2.set_xlim(.055, .845); a2.set_xticks(list(dep.threshold))
a2.set_title("Depletion around m\u2076A sites is a\ncomposition artefact", loc="left")
a2.legend(loc="lower left", frameon=False, fontsize=7, handlelength=1.4, labelspacing=.3, borderpad=.1)
set_frame(a2)

a3.axvline(1, ls=":", lw=.8, color=GRY, zorder=1)
rows = [("Uncorrected (Fisher)", or_raw, None, None, GRY),
        ("Logistic (3\u2032 UTR length,\nmutation count, isoforms)", or_adj, lo_, hi_, TEA)]
for k, (lab, v, l_, h_, col) in enumerate(rows):
    yy = len(rows) - 1 - k
    if l_ is not None:
        a3.plot([l_, h_], [yy, yy], lw=1.4, color=col, solid_capstyle="butt", zorder=3)
    a3.plot([v], [yy], "o", ms=5.0, color=col, zorder=4, markeredgecolor="white", markeredgewidth=.6)
    _m, _e = ("%.1e" % p_raw).split("e")
    txt = ("%.2f  P = %s \u00d7 10%s" % (v, _m, str(int(_e)).translate(str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b"))) if l_ is None
           else "%.2f (%.2f\u2013%.2f)  P = %.2f" % (v, l_, h_, p_adj))
    a3.text(v, yy + .20, txt, ha="center", fontsize=7, color=INK)
a3.set_yticks([1, 0]); a3.set_yticklabels([r[0] for r in rows], fontsize=7)
a3.set_xlabel("Odds ratio, Cancer Gene Census membership")
a3.set_xlim(.6, 2.4); a3.set_ylim(-.55, 1.62)
a3.set_title("Cancer-gene enrichment is explained\nby gene size", loc="left")
set_frame(a3)

# ── panel D: back-transformed medians with 95% bootstrap CIs (poster design)
rngB = np.random.default_rng(0)
def med_ci(v, B=4000):
    s = [np.median(rngB.choice(v, len(v), True)) for _ in range(B)]
    return 2 ** np.median(v), 2 ** np.percentile(s, 2.5), 2 ** np.percentile(s, 97.5)
FOLD = {k: med_ci(EX[k]) for k in ("background", "m6A", "rbp")}
rows4 = [("Any gene with a 3\u2032 UTR\n(n = %s)" % format(len(EX["background"]), ","), "background", BLU),
         ("m\u2076A-altering\n(n = %s)" % format(len(EX["m6A"]), ","), "m6A", TEA),
         ("m\u2076A-altering plus an\nRBP-binding change\n(n = %d)" % len(EX["rbp"]), "rbp", PUR)]
bg_m, bg_lo, bg_hi = FOLD["background"]
a4.axvspan(bg_lo, bg_hi, color=BLU, alpha=.10, zorder=1)
a4.axvline(bg_m, lw=.9, ls="--", color=BLU, zorder=2)
for i_, (lab, key, col) in enumerate(rows4):
    m_, lo2, hi2 = FOLD[key]
    y_ = 2 - i_
    a4.plot([lo2, hi2], [y_, y_], lw=1.8, color=col, solid_capstyle="butt", zorder=3)
    a4.plot([m_], [y_], marker="o", ms=5.2, color=col, zorder=4,
            markeredgecolor="white", markeredgewidth=.8)
    a4.text(m_, y_ + .20, "%.2f\u00d7" % m_, ha="center", va="bottom",
            fontsize=7, color=col, fontweight="bold")
for k_, key in ((1, "m6A"), (2, "rbp")):
    a4.text(FOLD[key][2] + .004, 2 - k_, "P = %.2f" % PE[key], ha="left", va="center",
            fontsize=7, color=GRY)
a4.set_yticks([2, 1, 0])
a4.set_yticklabels([r_[0] for r_ in rows4], fontsize=7, linespacing=1.30)
a4.set_xlim(1.42, 1.64); a4.set_ylim(-.72, 2.62)
a4.set_xticks([1.45, 1.50, 1.55, 1.60])
a4.set_xticklabels(["1.45\u00d7", "1.50\u00d7", "1.55\u00d7", "1.60\u00d7"])
a4.set_xlabel("Median tumour/normal\nexpression ratio", labelpad=2.5)
a4.text(bg_m + .004, -.66, "background level", fontsize=7, color=BLU, ha="left", va="bottom")
a4.set_title("Genes hit by these mutations are\nno more dysregulated", loc="left")
set_frame(a4)
a4.spines["left"].set_visible(False)
a4.tick_params(axis="y", length=0)

for ax_, L_ in zip((a3, a2, a1, a4), "ABCD"):
    panel_letter(ax_, L_, case="upper", dx=-0.20, dy=1.10)
fig5.savefig(out("Figure5_no_selective_signal.pdf"))
fig5.savefig(out("Figure5_no_selective_signal.png"), dpi=300)
