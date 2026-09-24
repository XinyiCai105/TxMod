# Reproducing the manuscript figures

The scripts in this directory render every figure in the manuscript (Figures 1–5,
Supplementary Figures S1 and S2, and the graphical abstract). They read only the
processed intermediate tables deposited at Zenodo (DOI: [PLACEHOLDER]), so the
figures can be regenerated without the licensed COSMIC records and without
rerunning m⁶A prediction.

## Steps

1. Download and unpack the Zenodo deposit. The `.gz` files can stay compressed.
2. Install the plotting dependencies:

       pip install pandas numpy scipy statsmodels matplotlib

3. Render the figures:

       TXMOD_DATA=/path/to/TxMod_zenodo python paper/make_figures.py

   Output goes to `paper/output/`, or to the directory named by `TXMOD_OUT`.
   The whole run takes under a minute.

## Files

| File | Renders |
|---|---|
| `build_fig1.py` | Figure 1, the workflow with its inset panels |
| `build_fig2.py` | Figure 2 |
| `build_fig3.py` | Figure 3 |
| `build_fig4.py` | Figure 4 and Supplementary Figure S1 |
| `build_fig5.py` | Figure 5 |
| `build_suppfig_s2.py` | Supplementary Figure S2 |
| `build_graphical_abstract.py` | Graphical abstract (drawn, reads no data) |
| `panel_lib.py` | Panels shared between Figure 1 and the figures it summarises |
| `figstyle.py` | Font sizes, spines and PDF font embedding |
| `paths.py` | Where input tables are read from and figures written to |
| `make_figures.py` | Runs all of the above |

Every value printed on a figure is computed from the input tables at run time.
None is typed into the scripts. Figures are drawn at their final print size,
so the font sizes in the scripts are the printed sizes.

## Check

Rendered from the deposit in a clean directory, all 16 output files match the
submitted figures: identical page size and text layer in every PDF, and identical
pixels in every PNG except anti-aliasing in the graphical abstract (0.003 % of
pixels).
