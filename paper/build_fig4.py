
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import data, out
from figstyle import apply_figure_style, set_frame, panel_letter
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_lib import draw_eclip

RB = pd.read_csv(data("rbp_binding_alterations_v2.tsv"), sep="\t")
EC = pd.read_csv(data("eclip_validation_v3_all20.tsv"), sep="\t").sort_values("enrichment", ascending=False).reset_index(drop=True)

# House convention, threaded from Figures 2 and 3: gain/increase orange, loss/decrease blue.
INC = "#C2622A"; DEC = "#1F5C8B"; GY = "#9AA5B1"; TEA = "#2E7D6E"; INK = "#1F2933"

up = RB[RB.delta_norm_score > 0].delta_norm_score.abs().values
dn = RB[RB.delta_norm_score < 0].delta_norm_score.abs().values
rng = np.random.default_rng(0)

plt.close("all")
apply_figure_style(sizes=(9, 8, 7))
fig4 = plt.figure(figsize=(6.35, 4.30))
gs = fig4.add_gridspec(1, 2, width_ratios=[1.0, 1.40], wspace=.42,
                       left=.095, right=.985, top=.87, bottom=.155)
axA = fig4.add_subplot(gs[0, 0]); axB = fig4.add_subplot(gs[0, 1])

# ── A: magnitude of every predicted occupancy change
for k, (v, col, lab) in enumerate(((up, INC, "increased"), (dn, DEC, "decreased"))):
    axA.scatter(np.full(len(v), k) + rng.normal(0, .085, len(v)), v, s=2.2, color=col,
                alpha=.30, linewidths=0, zorder=2)
    axA.plot([k - .30, k + .30], [np.median(v), np.median(v)], lw=2.0, color=INK,
             zorder=4, solid_capstyle="butt")
    axA.text(k, np.median(v) + .008, "%.3f" % np.median(v), ha="center", fontsize=7.5,
             fontweight="bold", color=INK, zorder=5)
axA.axhline(0.10, ls="--", lw=.8, color=GY, zorder=1)
axA.text(-.55, 0.0965, "calling threshold", ha="left", va="top", fontsize=7, color=GY)
axA.set_xticks([0, 1])
axA.set_xticklabels(["increased\n(n = %s)" % format(len(up), ","),
                     "decreased\n(n = %s)" % format(len(dn), ",")], fontsize=7.5, linespacing=1.35)
axA.set_ylabel("|\u0394 normalised PWM score|")
axA.set_ylim(.088, .27); axA.set_xlim(-.6, 1.6)
axA.set_title("Changes cluster near the\ncalling threshold", loc="left")
set_frame(axA)

# ── B (shared with the Figure 1 overview)
draw_eclip(axB, EC, frame=set_frame)

for ax_, L_ in zip((axA, axB), "AB"):
    panel_letter(ax_, L_, case="upper", dx=-0.22 if L_ == "A" else -0.30, dy=1.05)
fig4.savefig(out("Figure4_rbp_occupancy.pdf"), pad_inches=0)
fig4.savefig(out("Figure4_rbp_occupancy.png"), pad_inches=0, dpi=300)

# ── supplementary: the per-RBP record counts moved out of the main figure
plt.close("all")
apply_figure_style(sizes=(9, 8, 7))
cnt = (RB.assign(dir=np.where(RB.delta_norm_score > 0, "inc", "dec"))
         .groupby(["RBP", "dir"]).size().unstack(fill_value=0))
cnt["tot"] = cnt.sum(axis=1)
cnt = cnt.sort_values("tot")
figS = plt.figure(figsize=(4.9, 8.60))
axS = figS.add_subplot(111)
yy = np.arange(len(cnt))
axS.barh(yy, cnt.get("inc", 0), height=.66, color=INC, label="increased", zorder=3)
axS.barh(yy, -cnt.get("dec", 0), height=.66, color=DEC, label="decreased", zorder=3)
axS.axvline(0, lw=.8, color=INK, zorder=4)
axS.set_yticks(yy); axS.set_yticklabels(cnt.index, fontsize=7.0)
axS.set_ylim(-.8, len(cnt) - .2)
axS.set_xlabel("Binding-alteration records")
axS.set_title("All %d RBPs carrying at least one record" % len(cnt), loc="left")
axS.legend(loc="lower right", frameon=False, fontsize=7.5, handlelength=1.0,
           handleheight=.9, handletextpad=.4, labelspacing=.3)
set_frame(axS); axS.tick_params(axis="y", length=0)
figS.savefig(out("Supplementary_Figure_S1_rbp_panel_records.pdf"), pad_inches=0, bbox_inches="tight")
figS.savefig(out("Supplementary_Figure_S1_rbp_panel_records.png"), pad_inches=0, dpi=300, bbox_inches="tight")
print("A median increased %.3f (n=%d) / decreased %.3f (n=%d)" % (np.median(up), len(up), np.median(dn), len(dn)))
print("B nominally significant %s | BH q<0.05 %d | strongest %s q=%.3f"
      % (list(EC.RBP[EC.p_value < 0.05]), int((EC.p_adj_BH < 0.05).sum()),
         EC.RBP[EC.p_adj_BH.idxmin()], EC.p_adj_BH.min()))
print("supplementary figure %d RBPs, records in total %d" % (len(cnt), int(cnt.tot.sum())))
