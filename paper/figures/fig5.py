"""Fig. 5 (DRAFT, not finalized: the reach result is still being written up; run only when told).
Kick size against exact coupling distance from the nearest measured qubit (39 measured qubits packed in one region, three layouts),
reach run, median over untouched units per bin with 95 % bootstrap intervals, per processor, as committed by reach_run/analyze_reach.py
(RR1 'kick size at d=...'; size = kick x sign of the overall median, sign carried and printed by the analysis as -1 on every chip).
Inputs (read only): reach_run/results_reach.txt (committed output of analyze_reach.py on raw_reach_<chip>.npz / job_reach_<chip>.json).
Usage: python fig5.py [--out DIR]   (default DIR = this folder; writes fig5.pdf/png there)."""
import os, re, sys
import numpy as np
import matplotlib.pyplot as plt
import figstyle as fs

fs.setup()
OUT = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else fs.HERE
num = r'([-+]?[\d.]+(?:e[-+]?\d+)?)'
res, chip = {}, None
for line in open(os.path.join(fs.ROOT, 'reach_run', 'results_reach.txt')):
    m = re.match(r'^reach (ibm_\w+) job', line)
    if m: chip = m.group(1)
    if line.startswith('ACROSS CHIPS'): chip = None
    m = re.search(r'RR1 kick size at d(=|>=)(\d+)\s+%s \[%s, %s\] deg n=(\d+) -> (\S+)' % (num, num, num), line)
    if m and chip: res.setdefault(chip, []).append((int(m.group(2)), float(m.group(3)), float(m.group(4)), float(m.group(5)), int(m.group(6)), m.group(7)))
fig, ax = plt.subplots(figsize=(fs.COL1, 2.3))
cols = {'ibm_fez': fs.BLUE, 'ibm_kingston': fs.VERM, 'ibm_marrakesh': fs.GREEN}; mks = {'ibm_fez': 'o', 'ibm_kingston': 's', 'ibm_marrakesh': '^'}
for dx, c in zip((-0.15, 0, 0.15), fs.CHIPS):
    r = np.array([x[:4] for x in res[c]]); assert len(r) == 6, c
    ax.errorbar(r[:, 0] + dx, r[:, 1], yerr=[r[:, 1] - r[:, 2], r[:, 3] - r[:, 1]], fmt=mks[c], color=cols[c], ms=3.5, elinewidth=0.8, label=fs.SHORT[c])
    for x in res[c]: print('  %-14s d%s%d %+.4f [%+.4f, %+.4f] n=%d %s' % (c, '>=' if x[0] == 8 else '=', x[0], *x[1:]))
ax.axhline(0, color=fs.GREY, lw=0.5, zorder=0)
ax.set_xticks(range(3, 9)); ax.set_xticklabels(['3', '4', '5', '6', '7', r'$\geq$8'])
ax.set_xlabel('couplings to the nearest measured qubit'); ax.set_ylabel('kick size per readout (deg)')
ax.legend(loc='upper right')
for ext in ('pdf', 'png'): fig.savefig(os.path.join(OUT, 'fig5.%s' % ext), bbox_inches='tight', pad_inches=0.02)
print('wrote', os.path.join(OUT, 'fig5.pdf'))
