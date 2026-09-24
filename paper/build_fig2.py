
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import data, out
from figstyle import apply_figure_style, set_frame, panel_letter
import gc
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panel_lib import unicode_log_labels, draw_two_reference
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Patch

BL="#1F5C8B"; OR_="#C2622A"; GY="#9AA5B1"; DOT="#8A939E"
DATA=data("im6a_result_v2.tsv")
VAL=data("validation_gain_loss_two_reference.tsv")

# ── data: panel A distribution
chunks=[]
RRACH_COUNTS = data("Figure2BC_rrach_counts.tsv")
for c in pd.read_csv(DATA,sep="\t",usecols=["delta_prob"],chunksize=500_000,low_memory=False):
    chunks.append(pd.to_numeric(c.delta_prob,errors="coerce"))
d=pd.concat(chunks,ignore_index=True).dropna(); del chunks; gc.collect()
n_loss=int((d<=-0.10).sum()); n_gain=int((d>=0.10).sum())

F=pd.read_csv(VAL,sep="\t")

plt.close("all")
apply_figure_style(sizes=(9,8,7))
fig=plt.figure(figsize=(6.50,5.44))
gs=fig.add_gridspec(2,2,hspace=.52,wspace=.34,left=.095,right=.985,top=.945,bottom=.085)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1])
axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])

# A ── distribution
axA.hist(d,bins=np.linspace(-1,1,201),color=GY,log=True)
axA.axvline(-0.10,ls="--",lw=.8,color=BL); axA.axvline(0.10,ls="--",lw=.8,color=OR_)
axA.text(-0.14,3e4,"loss\n%s"%format(n_loss,","),ha="right",va="center",fontsize=8,color=BL)
axA.text(0.14,3e4,"gain\n%s"%format(n_gain,","),ha="left",va="center",fontsize=8,color=OR_)
axA.set_xlabel("\u0394 predicted m\u2076A probability")
axA.set_ylabel("Mutation\u2013transcript pairs")
unicode_log_labels(axA)
axA.set_xlim(-1,1); axA.set_ylim(1,2e6)
axA.set_yticks([1,1e2,1e4,1e6])
axA.set_title("Most mutations leave m\u2076A potential unchanged",loc="left")
set_frame(axA)

# B, C ── RRACH grammar
# Panels B and C read their counts from the motif table rather than carrying
# literals: percentages, odds ratios and intervals are derived here.
RR = pd.read_csv(RRACH_COUNTS, sep="\t")
_rr = []
for (lab, title, col), r_ in zip((("m\u2076A loss", "Losses disrupt an existing RRACH", BL),
                                  ("m\u2076A gain", "Gains create a new RRACH", OR_)),
                                 RR.itertuples(index=False)):
    _rr.append((100 * r_.n_events_altering_motif / r_.n_events,
                100 * r_.n_background_altering_motif / r_.n_background,
                r_.odds_ratio, lab, title, col))
for ax,(obs,bg,orv,lab,title,col) in zip((axB,axC), _rr):
    ax.bar([0,1],[obs,bg],width=.5,color=[col,GY],edgecolor="none")
    for x_,v_ in zip([0,1],[obs,bg]):
        ax.text(x_,v_+2.5,"%.1f%%"%v_,ha="center",fontsize=8)
    ax.set_xticks([0,1]); ax.set_xticklabels([lab,"Neutral\nmutations"])
    ax.set_ylim(0,100); ax.set_yticks([0,20,40,60,80,100])
    ax.set_xlim(-.7,1.7)
    ax.text(.5,92,"OR = %.1f"%orv,ha="center",fontsize=8,fontweight="bold",color=col)
    ax.set_title(title,loc="left"); set_frame(ax)
axB.set_ylabel("Mutations altering\nan RRACH motif (%)")
axC.set_ylabel("Mutations altering\nan RRACH motif (%)")

# D ── measured-reference validation (shared with the Figure 1 overview)
draw_two_reference(axD, F, frame=set_frame)

for ax_,L_ in zip((axA,axB,axC,axD),"ABCD"):
    panel_letter(ax_,L_,case="upper",dx=-0.20,dy=1.10)

fig.savefig(out("Figure2_m6A_both_directions.pdf"), pad_inches=0)
fig.savefig(out("Figure2_m6A_both_directions.png"), pad_inches=0,dpi=300)
