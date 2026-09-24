"""Panels drawn once and used by both their home figure and the Figure 1 overview.

Both callers execute the same function on the same tables, so the overview
cannot diverge from the figures it summarises.
"""
import numpy as np
from matplotlib.patches import Patch

BL = "#1F5C8B"; OR_ = "#C2622A"; GY = "#9AA5B1"; DOT = "#8A939E"; PUR = "#7B6AA0"
REFS = ["GLORI (both cell lines)", "m6A-Atlas v2.0"]
LAB = ["GLORI", "Atlas"]


def draw_two_reference(ax, F, frame=None, title=True):
    """Figure 2D: enrichment of predicted events in two measured references."""
    lo = [F[(F.reference == r) & (F.event_class == "loss")].iloc[0] for r in REFS]
    ga = [F[(F.reference == r) & (F.event_class == "gain (scored A pre-exists)")].iloc[0] for r in REFS]
    x = np.arange(len(REFS)); w = .34
    ax.bar(x - w / 2, [v.enrichment for v in lo], w, color=BL, zorder=3)
    ax.bar(x + w / 2, [v.enrichment for v in ga], w, color=OR_, zorder=3)
    for i, (a_, b_) in enumerate(zip(lo, ga)):
        ax.text(i - w / 2, a_.enrichment + 1.0, "%.0f\u00d7" % a_.enrichment, ha="center",
                fontsize=7, fontweight="bold", color=BL)
        ax.text(i + w / 2, b_.enrichment + 1.0, "%.0f\u00d7" % b_.enrichment, ha="center",
                fontsize=7, fontweight="bold", color=OR_)
    ax.axhline(1, color=DOT, lw=.8, ls=":", zorder=2)
    ax.set_xticks(x); ax.set_xticklabels(LAB)
    ax.set_ylabel("Enrichment over matched\nneutral background (fold)")
    ax.set_ylim(0, 40); ax.set_xlim(-.6, len(REFS) - .4)
    if title:
        ax.set_title("Affected adenosines carry measured methylation", loc="left")
    ax.legend(handles=[Patch(facecolor=BL, edgecolor="none", label="loss"),
                       Patch(facecolor=OR_, edgecolor="none", label="gain")],
              loc="upper right", frameon=False, fontsize=7, handlelength=1.1,
              handleheight=.9, handletextpad=.5, labelspacing=.4, borderpad=.1)
    if frame: frame(ax)
    return lo, ga


def draw_divergence_distribution(ax, SM, DV, frame=None, title=True):
    """Figure 3C: inter-isoform divergence over every mutation shared by two or more isoforms."""
    ax.hist(SM.divergence_score.dropna(), bins=np.linspace(0, 0.8, 161), color=PUR, log=True)
    for th, lab in ((0.10, "\u22650.10"), (0.50, "\u22650.50")):
        n_ = int(DV.loc[DV.threshold == th, "n_events_v2"].iloc[0])
        ax.axvline(th, lw=.7, ls="--", color=GY, zorder=3)
        ax.text(th + .012, 3.2e4, "%s\n%s" % (lab, format(n_, ",")), fontsize=7, color=GY,
                va="top", ha="left", linespacing=1.25)
    ax.set_xlabel("Inter-isoform divergence, max(\u0394prob) \u2212 min(\u0394prob)")
    ax.set_ylabel("Shared mutations")
    ax.set_xlim(0, .8); ax.set_ylim(1, 3e5); ax.set_yticks([1, 1e2, 1e4])
    unicode_log_labels(ax)
    if title:
        ax.set_title("Shared mutations mostly agree", loc="left")
    if frame: frame(ax)


def draw_eclip(ax, EC, frame=None, title=True, label_fs=7.0):
    """Figure 4B: predicted alteration sites against measured ENCODE eCLIP peaks."""
    TEA = "#2E7D6E"; INK = "#1F2933"
    EC = EC.sort_values("enrichment", ascending=False).reset_index(drop=True)
    y = np.arange(len(EC))[::-1]
    sig = EC.p_value < 0.05
    undef = (EC.n_in_peak == 0) & (EC.bg_rate == 0)
    for i, r_ in EC.iterrows():
        col = TEA if sig[i] else GY
        ax.plot([r_.enrichment], [y[i]], marker="o", ms=4.4,
                color="white" if undef[i] else col, markeredgecolor=col,
                markeredgewidth=1.0, zorder=3)
    ax.axvline(1, ls=":", lw=.8, color=INK, zorder=1)
    for i, r_ in EC.iterrows():
        if r_.p_adj_BH < 0.2:
            ax.text(r_.enrichment * 1.25, y[i], "q = %.2f" % r_.p_adj_BH, va="center",
                    fontsize=7, color=TEA, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(["%s (%d)" % (r_.RBP, r_.n_events) for _, r_ in EC.iterrows()], fontsize=label_fs)
    ax.set_xscale("symlog", linthresh=0.3)
    ax.set_xticks([0, 0.3, 1, 3, 10]); ax.set_xticklabels(["0", "0.3", "1", "3", "10"])
    ax.set_xlim(-0.05, 26); ax.set_ylim(-.8, len(EC) - .2)
    ax.set_xlabel("Enrichment in measured eCLIP peaks\n(vs events predicted for other RBPs)", linespacing=1.3)
    if title:
        ax.set_title("No protein survives correction across the 20 tests", loc="left")
    if frame: frame(ax)
    ax.tick_params(axis="y", length=0)



def draw_eclip_vertical(ax, EC, frame=None, label_fs=7.0):
    """Same content as draw_eclip, proteins on the x axis with rotated labels.

    Used for the Figure 1 module inset, where 20 rows at 7 pt do not fit
    vertically but do fit horizontally.
    """
    TEA = "#2E7D6E"; INK = "#1F2933"
    EC = EC.sort_values("enrichment", ascending=False).reset_index(drop=True)
    x = np.arange(len(EC))
    sig = EC.p_value < 0.05
    undef = (EC.n_in_peak == 0) & (EC.bg_rate == 0)
    for i, r_ in EC.iterrows():
        col = TEA if sig[i] else GY
        ax.plot([x[i]], [r_.enrichment], marker="o", ms=4.0,
                color="white" if undef[i] else col, markeredgecolor=col,
                markeredgewidth=1.0, zorder=3)
    ax.axhline(1, ls=":", lw=.8, color=INK, zorder=1)
    for i, r_ in EC.iterrows():
        if r_.p_adj_BH < 0.2:
            # Beside the point rather than above it: stacked above, the two labels of
            # adjacent proteins collide with each other and with the y-axis ticks.
            ax.text(x[i] + .45, r_.enrichment * 1.9, "q = %.2f" % r_.p_adj_BH, ha="left", va="bottom",
                    fontsize=label_fs, color=TEA, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(["%s (%d)" % (r_.RBP, r_.n_events) for _, r_ in EC.iterrows()],
                       fontsize=label_fs, rotation=90, ha="center", va="top")
    ax.set_yscale("symlog", linthresh=0.3)
    ax.set_yticks([0, 0.3, 1, 3, 10]); ax.set_yticklabels(["0", "0.3", "1", "3", "10"])
    ax.set_ylim(-0.05, 60); ax.set_xlim(-.8, len(EC) - .2)
    ax.set_ylabel("Enrichment in\nmeasured eCLIP peaks")
    if frame is not None: frame(ax)
    return ax

def load_recurrence(cache=None, im6a=None, counts=None):
    """Sample-recurrence of 3' UTR mutations and of the high-confidence subset.

    Reads the cached summary (recurrence_counts.json) when present. Otherwise
    recomputes it from the per-record prediction table and the COSMIC
    sample-count table.
    """
    import json, os
    import pandas as pd
    from paths import data
    cache = cache or data("recurrence_counts.json", must_exist=False)
    if os.path.exists(cache):
        return json.load(open(cache))
    im6a = im6a or data("im6a_result_v2.tsv")
    counts = counts or data("cosv_sample_counts.tsv")
    if os.path.exists(cache):
        return json.load(open(cache))
    sc_all = pd.read_csv(counts, sep="\t")
    utrv, evv = set(), set()
    for c in pd.read_csv(im6a, sep="\t", usecols=["record_key", "delta_prob"], chunksize=500_000):
        utrv |= set(c.record_key.str.split("|").str[0])
        h = c[c.delta_prob.abs() >= 0.10]
        if len(h):
            evv |= set(h.record_key.str.split("|").str[0])
    sc = sc_all[sc_all.COSV.isin(utrv)]; sub = sc[sc.COSV.isin(evv)]
    R = {k: dict(n=int(g.COSV.nunique()),
                 pct=[float(100 * (g.n_samples == 1).mean()),
                      float(100 * (g.n_samples == 2).mean()),
                      float(100 * (g.n_samples >= 3).mean())])
         for k, g in (("all", sc), ("events", sub))}
    os.makedirs(os.path.dirname(cache) or ".", exist_ok=True)
    json.dump(R, open(cache, "w"))
    return R


def draw_recurrence(ax, R, frame=None, title=True):
    """Figure 5A: how many tumours each mutation is seen in."""
    PURP = "#7B6AA0"; RED = "#C1483A"
    # Three short lines rather than two long ones: at 7 pt the two-line form ran the
    # neighbouring tick labels together once the panel narrowed.
    grp = [("All 3\u2032 UTR\nmutations\n(n = %s)" % format(R["all"]["n"], ","), R["all"]["pct"]),
           ("High-confidence\nm\u2076A-altering\n(n = %s)" % format(R["events"]["n"], ","), R["events"]["pct"])]
    xs = np.arange(2); bot = np.zeros(2)
    for k, (lab, col) in enumerate([("1 sample", GY), ("2 samples", PURP), ("\u22653 samples", RED)]):
        v = np.array([g[1][k] for g in grp])
        ax.bar(xs, v, .5, bottom=bot, color=col, label=lab); bot += v
    for k, (lab, p) in enumerate(grp):
        ax.text(k, p[0] / 2, "%.1f%%" % p[0], ha="center", va="center", fontsize=7,
                color="white", fontweight="bold")
    ax.set_xticks(xs); ax.set_xticklabels([g[0] for g in grp], fontsize=7)
    ax.set_ylabel("Mutations (%)"); ax.set_xlim(-.6, 1.6)
    # The bars reach 100, so the key goes in the empty band above them rather than
    # on top of the data.
    ax.set_ylim(0, 124); ax.set_yticks([0, 20, 40, 60, 80, 100])
    if title:
        ax.set_title("Almost all occur in a single tumour", loc="left")
    ax.legend(loc="upper center", ncol=3, frameon=False, fontsize=7.0, handlelength=.9,
              handleheight=.9, handletextpad=.35, columnspacing=.9, borderpad=0)
    if frame: frame(ax)

SUP_TR = str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b")


def unicode_log_labels(ax, axis="y"):
    """Log tick labels as Unicode superscripts rather than mathtext.

    Mathtext renders the exponent at 0.7x the base size, which puts a 7 pt axis
    at 4.9 pt in the PDF and below the publisher's 7 pt floor. Unicode
    superscript glyphs are drawn at the full tick size.
    """
    a = ax.yaxis if axis == "y" else ax.xaxis
    ticks = [t for t in a.get_ticklocs() if t > 0]
    labs = ["10" + str(int(round(np.log10(t)))).translate(SUP_TR) for t in ticks]
    a.set_ticks(ticks)
    a.set_ticklabels(labs)
