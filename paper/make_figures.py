"""Render every manuscript figure from the deposited tables.

    TXMOD_DATA=/path/to/zenodo_deposit python paper/make_figures.py

Each build script is run in a fresh namespace. Output goes to TXMOD_OUT
(default paper/output/).
"""
import os, runpy, sys, time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SCRIPTS = ["build_fig1.py", "build_fig2.py", "build_fig3.py", "build_fig4.py",
           "build_fig5.py", "build_suppfig_s2.py", "build_graphical_abstract.py"]

if __name__ == "__main__":
    failed = []
    for s in SCRIPTS:
        t0 = time.time()
        try:
            plt.close("all")
            runpy.run_path(os.path.join(HERE, s), run_name="__main__")
            print(f"[ok]   {s}  ({time.time() - t0:.0f} s)")
        except Exception as e:
            failed.append(s)
            print(f"[fail] {s}: {type(e).__name__}: {e}")
    sys.exit(1 if failed else 0)
