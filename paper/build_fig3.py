
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import data, out
from figstyle import apply_figure_style, set_frame, panel_letter
import os, numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_lib import draw_divergence_distribution

SM = pd.read_csv(data("Figure3C_inter_isoform_divergence_summary.tsv"), sep="\t")
T4 = pd.read_csv(data("Figure3AB_representative_cases.tsv"), sep="\t")
DV = pd.read_csv(data("Figure4D_divergence_v1_vs_v2.tsv"), sep="\t")

# ── panel D: every isoform record of the high-divergence mutations (no subsample)
hi = set(SM[SM.divergence_score >= 0.10].mutation_key)
use = ["ID", "locus", "ref_base", "alt_base", "delta_prob", "offset_1based", "utr_len"]
rec = []
for c in pd.read_csv(data("im6a_result_v2.tsv"), sep="\t", usecols=use, chunksize=400_000, low_memory=False):
    c["mutation_key"] = (c.ID.astype(str) + "|" + c.locus.astype(str) + "|"
                         + c.ref_base.astype(str) + ">" + c.alt_base.astype(str))
    rec.append(c[c.mutation_key.isin(hi)])
S = pd.concat(rec, ignore_index=True); del rec
S["d"] = pd.to_numeric(S.delta_prob, errors="coerce")
S["rel_pos"] = 100 * pd.to_numeric(S.offset_1based, errors="coerce") / pd.to_numeric(S.utr_len, errors="coerce")
S = S.dropna(subset=["d", "rel_pos"])
bins = [(0, 10), (10, 30), (30, 50), (50, 70), (70, 90), (90, 100)]
BB = pd.DataFrame([dict(bin="%d\u2013%d" % (lo, hi_), n=len(b),
                        frac=100 * (b.d.abs() >= 0.10).mean())
                   for lo, hi_ in bins
                   for b in [S[(S.rel_pos >= lo) & ((S.rel_pos < hi_) if hi_ < 100 else (S.rel_pos <= 100))]]])
mid = S[(S.rel_pos >= 10) & (S.rel_pos < 90)]; ends = S[(S.rel_pos < 10) | (S.rel_pos >= 90)]
orv, pv = fisher_exact([[int((mid.d.abs() >= 0.10).sum()), len(mid) - int((mid.d.abs() >= 0.10).sum())],
                        [int((ends.d.abs() >= 0.10).sum()), len(ends) - int((ends.d.abs() >= 0.10).sum())]])

# ── panels A, B: the two named cases, every isoform labelled
KEEP = ["CFTR", "PIK3R1"]
order = (T4.drop_duplicates("mutation_key")
           .merge(SM[["mutation_key", "divergence_score"]], on="mutation_key", how="left")
           .set_index("gene_name").loc[KEEP].reset_index())

BL = "#1F5C8B"; OR_ = "#C2622A"; GY = "#9AA5B1"; PUR = "#7B6AA0"; INK = "#1F2933"
plt.close("all")
apply_figure_style(sizes=(9, 8, 7))
fig3 = plt.figure(figsize=(6.22, 6.80))
gs = fig3.add_gridspec(2, 2, height_ratios=[1.62, 1.0], hspace=.56, wspace=.42,
                       left=.10, right=.975, top=.93, bottom=.075)
axs = [fig3.add_subplot(gs[0, 0]), fig3.add_subplot(gs[0, 1])]
axC = fig3.add_subplot(gs[1, 0]); axD = fig3.add_subplot(gs[1, 1])

for ax, row in zip(axs, order.itertuples(index=False)):
    g = T4[T4.mutation_key == row.mutation_key].copy()
    g["d"] = pd.to_numeric(g.delta_prob, errors="coerce")
    g = g.sort_values("d").reset_index(drop=True)
    n = len(g)
    col = [BL if v < 0 else OR_ for v in g.d]
    ax.barh(range(n), g.d, height=.62, color=col, edgecolor="none", zorder=3)
    ax.axvline(0, lw=.8, color=INK, zorder=2)
    for t_ in (-0.10, 0.10):
        ax.axvline(t_, lw=.7, ls=":", color=GY, zorder=2)
    lab = ["%s  %s/%s" % (t.split(".")[0][-7:], format(int(o), ","), format(int(u), ","))
           for t, o, u in zip(g.transcript_id, g.offset_1based, g.utr_len)]
    ax.set_yticks(range(n)); ax.set_yticklabels(lab, fontsize=7.6 if n <= 6 else 7.0)
    ax.set_ylim(-.8, n - .2); ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("\u0394 predicted m\u2076A probability")
    ax.set_title("%s  %s %s>%s\n%d isoforms, divergence %.2f"
                 % (row.gene_name, row.locus, row.ref_base, row.alt_base, n, row.divergence_score),
                 loc="left", linespacing=1.45)
    set_frame(ax); ax.tick_params(axis="y", length=0)

# ── panel C (shared with the Figure 1 overview)
draw_divergence_distribution(axC, SM, DV, frame=set_frame)

# ── panel D: effect size against position within each isoform's own 3' UTR
axD.bar(range(len(BB)), BB.frac, width=.62,
        color=[GY if i in (0, len(BB) - 1) else BL for i in range(len(BB))], edgecolor="none")
# Percentage and sample size both sit above the bar, so the axis carries only the
# position bin and the n = labels cannot run into one another at 7 pt.
for i, r_ in enumerate(BB.itertuples(index=False)):
    axD.text(i, r_.frac + 6.2, "%.0f%%" % r_.frac, ha="center", fontsize=8, fontweight="bold")
axD.set_xticks(range(len(BB)))
axD.set_xticklabels([r_.bin for r_ in BB.itertuples(index=False)], fontsize=7.0)
axD.set_xlabel("Mutation position within the\ntranscript-specific 3\u2032 UTR (%)", linespacing=1.3)
axD.set_ylabel("Isoform records with\n|\u0394| \u2265 0.10 (%)")
axD.set_ylim(0, 100); axD.set_xlim(-.7, len(BB) - .3)
axD.set_title("High-divergence mutations only (n = 1,539)", loc="left")
set_frame(axD)

for ax_, L_ in zip(axs + [axC, axD], "ABCD"):
    panel_letter(ax_, L_, case="upper", dx=-0.235, dy=1.07)
fig3.savefig(out("Figure3_isoform_divergence.pdf"))
fig3.savefig(out("Figure3_isoform_divergence.png"), dpi=300)
print("bins", [(r.bin, int(r.n), round(r.frac, 1)) for r in BB.itertuples(index=False)])
print("interior %.1f%% (n=%d) vs ends %.1f%% (n=%d) | OR %.1f" %
      (100 * (mid.d.abs() >= 0.10).mean(), len(mid), 100 * (ends.d.abs() >= 0.10).mean(), len(ends), orv))
print("records", len(S), "| mutations", len(hi))
