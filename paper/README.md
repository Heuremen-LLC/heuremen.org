# The paper, built from source

`main.tex`, `refs.bib` and `figures/` build the PDF published at https://heuremen.org/a-readout-is-not-an-idle.pdf. `./build.sh` runs tectonic with the build date pinned to the publication date, because the title block prints `\today`. `python3 verify_paper.py` compares the rebuilt PDF with the published one: same page count and the same text, ignoring spacing and ligature glyphs. CI does both on every push (.github/workflows/paper.yml).

`figures/fig*.py` draw the figures from the run data. Results text and code are under `/experiments/` in this repository. The large raw arrays (`.npz`) are in the raw bundle of the Zenodo record (concept DOI 10.5281/zenodo.23169755): unpack it and run

```
READOUT_EXPERIMENTS=/path/to/unpacked/experiments python3 figures/fig2.py
```

`figstyle.py` maps each script's short run name (control_run, deciding_runs, ...) to the dated folder under `experiments/`. Figures 1 to 4 regenerate pixel for pixel identical to the ones in the paper (checked, maximum pixel difference 0). Figure 1 needs only this repository.

Heuremen LLC. The paper, run data and code are open with attribution (CC BY 4.0). See LICENSE.txt.
