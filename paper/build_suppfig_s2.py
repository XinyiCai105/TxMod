import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import data, out
from figstyle import apply_figure_style, set_frame, panel_letter
import os, sys, numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
    }
    mpl.rcParams.update(rc)


def set_frame(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_letter(ax, letter, case="upper", dx=-0.16, dy=1.02, fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    L = letter.upper() if case == "upper" else letter.lower()
    ax.text(dx, dy, L, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


CGC_T = data("cgc_gene_level_table_v2.tsv")
TRI_T = data("trinuc_context_freq.tsv")

# ── panel A: Census genes are larger, more mutated and more isoform-rich
GL = pd.read_csv(CGC_T, sep="\t")
cg, ot = GL[GL.is_CGC == 1], GL[GL.is_CGC == 0]
conf = []
for lab, col in (("Total 3\u2032 UTR length", "total_utr_length"),
                 ("Mutations in 3\u2032 UTR", "n_mutations"),
                 ("Transcripts with 3\u2032 UTR", "n_transcripts")):
    x, y = cg[col].dropna(), ot[col].dropna()
    conf.append((lab, x.median(), y.median(), x.median() / y.median(),
                 mannwhitneyu(x, y, alternative="two-sided")[1]))

# ── panel B: somatic mutability by trinucleotide context.
TF = pd.read_csv(TRI_T, sep="\t")
t5 = TF[TF.threshold == 0.5]
share = float(t5[t5.trinucleotide.isin(["GAC", "AAC"])].m6A_site_freq.sum())
names = ["GAC", "AAC", "Other\n3\u2032 UTR", "ACG\n(CpG)", "CGT\n(CpG)"]
rates = [0.0099, 0.0105, 0.030, 0.100, 0.100]

TEAL = "#2E7D6E"; RED = "#C1483A"; GREY = "#9AA5B1"; INK = "#1F2933"
plt.close("all")
apply_figure_style(sizes=(9, 8, 7))
fig, (axa, axb) = plt.subplots(1, 2, figsize=(6.30, 3.40),
                               gridspec_kw=dict(width_ratios=[1.0, 0.88], wspace=.42))
fig.subplots_adjust(left=.185, right=.99, top=.79, bottom=.24)

yy = np.arange(len(conf))[::-1]
axa.barh(yy, [c[3] for c in conf], color=GREY, alpha=.65, height=.5)
for y, (lab, m1, m2, r_, p_) in zip(yy, conf):
    axa.text(r_ + .07, y + .11, "%s vs %s" % (format(int(m1), ","), format(int(m2), ",")),
             va="center", fontsize=7, color=INK)
    _m, _e = ("%.1e" % p_).split("e")
    axa.text(r_ + .07, y - .22, "P = %s \u00d7 10%s" % (_m, str(int(_e)).translate(str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b"))),
             va="center", fontsize=7, color=GREY)
axa.axvline(1.0, ls=":", lw=.9, color=GREY)
axa.set_yticks(yy); axa.set_yticklabels([c[0] for c in conf], fontsize=7)
axa.set_xlim(0, 3.3); axa.set_xlabel("Median ratio, Census vs other genes")
axa.set_title("Why cancer-gene enrichment\nappears. Census genes offer\nmore opportunity", loc="left")
set_frame(axa)

for i, (v, c) in enumerate(zip(rates, [TEAL, TEAL, GREY, RED, RED])):
    axb.plot([i, i], [0.006, v], color=c, lw=1.4, alpha=.85, zorder=2)
    axb.scatter([i], [v], s=42, color=c, zorder=3, edgecolor="white", linewidths=.7)
    axb.text(i, v * 1.30, "%.3f" % v, ha="center", va="bottom", fontsize=7, color=INK)
axb.set_yscale("log"); axb.set_ylim(0.006, 0.42); axb.set_xlim(-0.6, 4.6)
axb.set_xticks(range(5)); axb.set_xticklabels(names, fontsize=7)
axb.set_yticks([0.01, 0.1]); axb.set_yticklabels(["0.01", "0.10"])
axb.set_ylabel("Somatic mutation rate (log scale)")
axb.set_title("Why depletion appears. m\u2076A\ncontexts are among the least\nmutable trinucleotides", loc="left")
axb.text(0.02, 0.97, "GAC + AAC = %.0f%% of predicted m\u2076A sites" % (100 * share),
         transform=axb.transAxes, ha="left", va="top", fontsize=7, color=GREY)
set_frame(axb)
for ax_, L_ in zip((axa, axb), "AB"):
    panel_letter(ax_, L_, case="upper", dx=-0.22, dy=1.14)
fig.savefig(out("Supplementary_Figure_S2_artefact_mechanisms.pdf"))
fig.savefig(out("Supplementary_Figure_S2_artefact_mechanisms.png"), dpi=300)
