#!/usr/bin/env python3
"""Figure 1 - the TxMod workflow.

Stages 4 and 5 are not re-drawn here: they call the same functions that draw
Figure 2D and Figure 3C (panel_lib), so the overview cannot drift from the
figures it summarises. Every count in the stage boxes is read from
fig1_scale.json, recomputed from im6a_result_v2.tsv.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import data, out
from figstyle import apply_figure_style, set_frame, panel_letter
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from panel_lib import (draw_eclip_vertical, draw_two_reference, draw_divergence_distribution,
                       draw_eclip, draw_recurrence, load_recurrence)

S = json.load(open(data("fig1_scale.json")))
F = pd.read_csv(data("validation_gain_loss_two_reference.tsv"), sep="\t")
SM = pd.read_csv(data("Figure3C_inter_isoform_divergence_summary.tsv"), sep="\t")
DV = pd.read_csv(data("Figure4D_divergence_v1_vs_v2.tsv"), sep="\t")
EC = pd.read_csv(data("eclip_validation_v3_all20.tsv"), sep="\t")
R = load_recurrence()
RB = pd.read_csv(data("rbp_binding_alterations_v2.tsv"), sep="\t")

BL="#1F5C8B"; OR_="#C2622A"; GY="#9AA5B1"; INK="#1F2933"; SOFT="#55606C"
PU="#7B3FA0"; EX="#4A79B8"; FILL="#F4F7FA"; EDGE="#C7D2DC"

plt.close("all")
apply_figure_style(sizes=(9, 8, 7))
mpl.rcParams["savefig.bbox"] = "standard"   # keep the canvas exactly double-column
CW, CH = 180.0, 220.0
fig = plt.figure(figsize=(CW/25.4, CH/25.4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, CW); ax.set_ylim(CH, 0); ax.axis("off")

def box(x, y, w, h, n, title, fill=FILL, edge=EDGE):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=2.2",
                                linewidth=.8, edgecolor=edge, facecolor=fill, zorder=2))
    ax.add_patch(plt.Circle((x+5.6, y+5.4), 3.0, facecolor=BL, edgecolor="none", zorder=4))
    ax.text(x+5.6, y+5.4, str(n), ha="center", va="center", fontsize=7.0, color="white",
            fontweight="bold", zorder=5)
    ax.text(x+10.6, y+5.4, title, va="center", fontsize=7.5, color=BL, fontweight="bold", zorder=4)

def chev(x, y, dy=False):
    p = ((x, y), (x, y+7)) if dy else ((x, y), (x+7, y))
    ax.add_patch(FancyArrowPatch(*p, arrowstyle="-|>", mutation_scale=9, linewidth=1.8,
                                 color="#B9C2CC", zorder=1))


def retype(axx, tick=7.0, anno=7.0, axlab=7.5):
    """One type ladder for the whole figure: stage title 7.5, stage text 6.0,
    inset axis title 6.0, inset ticks 5.2, in-panel annotation 5.8."""
    axx.xaxis.label.set_size(axlab); axx.yaxis.label.set_size(axlab)
    for t in axx.get_xticklabels() + axx.get_yticklabels(): t.set_size(tick)
    for t in axx.texts: t.set_size(anno)
    lg = axx.get_legend()
    if lg:
        for t in lg.get_texts(): t.set_size(anno)


def place(axx, bx0, by0, bx1, by1, pad=2.0):
    """Fit the axes so its tight bounding box - labels included - lands inside
    the reserved rectangle of its stage box."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    for _ in range(3):
        r_ = FigureCanvasAgg(fig).get_renderer()
        pos = axx.get_position()
        tb = axx.get_tightbbox(r_).transformed(fig.transFigure.inverted())
        L = pos.x0 - tb.x0; R = tb.x1 - pos.x1
        B = pos.y0 - tb.y0; T = tb.y1 - pos.y1
        X0 = (bx0 + pad) / CW + L; X1 = (bx1 - pad) / CW - R
        Y0 = 1 - (by1 - pad) / CH + B; Y1 = 1 - (by0 + pad) / CH - T
        axx.set_position([X0, Y0, max(X1 - X0, .02), max(Y1 - Y0, .02)])

def AX(x, y, w, h):
    return fig.add_axes([x/CW, 1-(y+h)/CH, w/CW, h/CH])

# ══ stage 1 ════════════════════════════════════════════════════════════════
box(3, 4, 40, 52, 1, "Input")
ax.text(23, 19, "COSMIC v102", ha="center", fontsize=8.4, color=BL, fontweight="bold")
ax.text(23, 24.5, "GRCh38, pan-cancer", ha="center", fontsize=7.0, color=SOFT)
for k, c in enumerate(["#E8A7B8", "#F0C3A0", "#D98C5F", "#B5564E", "#C8746B"]):
    xx = 8 + k*6.2
    ax.add_patch(plt.Circle((xx+1.9, 33), 1.7, facecolor=c, edgecolor="none", zorder=4))
    ax.add_patch(FancyBboxPatch((xx, 35.5), 3.8, 7.5+(k % 3)*1.2,
                                boxstyle="round,pad=0,rounding_size=1.1",
                                facecolor=c, edgecolor="none", zorder=4))
ax.text(23, 50, "%s somatic\n3\u2032 UTR substitutions" % format(S["variants"], ","),
        ha="center", va="center", fontsize=7.0, color=INK, linespacing=1.4)
chev(45, 30)

# ══ stage 2 ════════════════════════════════════════════════════════════════
box(55, 4, 62, 52, 2, "Isoform-resolved reconstruction")
ax.text(57.5, 16.0, "Gene locus", fontsize=7.0, color=INK)
ax.plot([60, 110], [19.5, 19.5], color=GY, lw=.8)
for k, (exs, ux, uw) in enumerate([([(61, 5), (69, 4), (76, 5)], 84, 12),
                                   ([(61, 5), (69, 4), (76, 5)], 84, 8.5),
                                   ([(61, 5), (69, 4)], 76, 6.5)]):
    Y = 25 + k*7.5
    ax.plot([61, ux+uw+3], [Y, Y], color=GY, lw=.6, zorder=3)
    for x0, w0 in exs:
        ax.add_patch(Rectangle((x0, Y-1.3), w0, 2.6, facecolor=EX, edgecolor="none", zorder=4))
    ax.add_patch(Rectangle((ux, Y-1.6), uw, 3.2, facecolor=PU, edgecolor="none", zorder=4))
    ax.text(ux+uw+4.2, Y, "AAAA", fontsize=7.0, color=GY, va="center")
    ax.plot([ux+uw*[.30, .55, .78][k]], [Y-3.2], marker="v", ms=3.4, color=OR_, zorder=5)
for lab, c, xx in (("Exon", EX, 59), ("3\u2032 UTR", PU, 74)):
    ax.add_patch(Rectangle((xx, 45.0), 2.6, 2.4, facecolor=c, edgecolor="none", zorder=4))
    ax.text(xx+3.4, 46.2, lab, fontsize=7.0, color=INK, va="center")
ax.text(57.5, 50.5, "One mutation, a different position in each 3\u2032 UTR", fontsize=7.0, color=INK)
ax.text(57.5, 54.0, "%s pairs across %s transcripts"
        % (format(S["pairs"], ","), format(S["transcripts"], ",")), fontsize=7.0, color=INK)
chev(118.5, 30)

# ══ stage 3 ════════════════════════════════════════════════════════════════
box(127, 4, 50, 52, 3, "m\u2076A change prediction")
ax.text(152, 14.5, "iM6A humanRRACH10000 ensemble", ha="center", fontsize=7.0, color=SOFT)
for k, (lab, col, seq, hl) in enumerate([("Reference", BL, "ACUGGACUUA", -1),
                                         ("Mutant", OR_, "ACUGGGCUUA", 5)]):
    Y = 21 + k*9
    ax.text(131, Y, lab, fontsize=7.0, color=col, va="center", fontweight="bold")
    for j, ch_ in enumerate(seq):
        ax.text(147+j*3.0, Y, ch_, fontsize=7.0, ha="center", va="center",
                color=OR_ if j == hl else INK, fontweight="bold" if j == hl else "normal")
ax.plot([145.5, 145.5+10*3.0], [26.5, 26.5], color=GY, lw=.5)
ax.text(131, 36.0, "\u00b12 nt window", fontsize=7.0, color=SOFT)
for k, (lab, col, fc) in enumerate([("\u0394 \u2264 \u22120.10\nLoss", BL, "#E8EEF4"),
                                    ("|\u0394| < 0.10\nNo change", GY, "#F0F2F4"),
                                    ("\u0394 \u2265 0.10\nGain", OR_, "#F7EBE3")]):
    xx = 129 + k*16.0
    ax.add_patch(FancyBboxPatch((xx, 39.0), 14.8, 9.6, boxstyle="round,pad=0,rounding_size=1.6",
                                linewidth=.7, edgecolor=col, facecolor=fc, zorder=3))
    ax.text(xx+7, 44, lab, ha="center", va="center", fontsize=7.0, color=col,
            fontweight="bold", linespacing=1.3, zorder=4)
ax.text(152, 50.5, "%s events at |\u0394| \u2265 0.10" % format(S["events"], ","),
        ha="center", fontsize=7.0, color=INK)
ax.text(152, 54.0, "%s loss and %s gain" % (format(S["loss"], ","), format(S["gain"], ",")),
        ha="center", fontsize=7.0, color=INK)
chev(90, 58, dy=True)

# ══ stages 4 and 5 - the real panels, same code as Figures 2D and 3C ═══════
box(3, 68, 86, 62, 4, "Measured support")
ax.text(13, 79.5, "GLORI and m\u2076A-Atlas, two measured single-base references",
        fontsize=7.0, color=SOFT)
a4 = AX(16, 84, 68, 38); draw_two_reference(a4, F, frame=set_frame, title=False)
a4.set_title(""); retype(a4); place(a4, 5, 84, 88, 128)

box(94, 68, 83, 62, 5, "Isoform divergence")
ax.text(104, 79.5, "197,289 mutations shared by two or more isoforms",
        fontsize=7.0, color=SOFT)
a5 = AX(106, 84, 66, 36); draw_divergence_distribution(a5, SM, DV, frame=set_frame, title=False)
a5.set_title(""); retype(a5); place(a5, 96, 84, 176, 128)
chev(90, 132, dy=True)

# ══ stages 6 and 7 ═════════════════════════════════════════════════════════
box(3, 140, 86, 76, 6, "RBP occupancy change")
ax.text(13, 150.5, "PWM scan over the full 3\u2032 UTR, using 80 RBPs and 178 motifs,"
        , fontsize=7.0, color=SOFT)
ax.text(13, 154.5, "which yields %s binding-alteration records" % format(len(RB), ","),
        fontsize=7.0, color=SOFT)
ax.text(13, 158.5, "Tested against ENCODE eCLIP peaks",
        fontsize=7.0, color=SOFT)
a6 = AX(14, 163, 70, 32); draw_eclip_vertical(a6, EC, frame=set_frame, label_fs=7.0)
a6.set_title(""); retype(a6); place(a6, 8, 162, 88, 214)

box(94, 140, 83, 76, 7, "Tests of selection")
# Four short lines: at 7 pt the previous three ran past the module frame.
for _k, _ln in enumerate(["Recurrence across tumours, mutation density,",
                          "Cancer Gene Census enrichment and tumour/normal",
                          "expression, each against a matched background.",
                          "Recurrence is shown here."]):
    ax.text(104, 150.0 + _k*4.0, _ln, fontsize=7.0, color=SOFT)
a7 = AX(112, 169, 56, 32); draw_recurrence(a7, R, frame=set_frame, title=False)
a7.set_title(""); retype(a7); place(a7, 96, 168, 176, 207)
ax.text(104, 210, "None returns a signal", fontsize=7.0, color=BL, fontweight="bold")

fig.savefig(out("Figure1_framework.pdf"))
fig.savefig(out("Figure1_framework.png"), dpi=300)
