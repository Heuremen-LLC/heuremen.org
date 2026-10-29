"""Fig. 4. Cost of the readouts to untouched qubits, idle minus readout (positive = less after readouts than after the equal idle),
median over units with 95 % bootstrap intervals (2000 resamples), both copies pooled, per processor.
(a-c) coherence |<X> + i<Y>| of untouched qubits wired to a measured qubit (86 per processor), raw-readout run, N = 2, 4, 8,
      all shots (analyze_iq.py R2) and clean shots (R1: every measured neighbour verified in its |0> cluster on every readout).
(d-f) concurrence of untouched Bell pairs (wired, idle concurrence > 0), pairs run, N = 2, 4, all shots (analyze_pairs_iq.py R3)
      and clean shots (R1).
Inputs (read only; committed outputs parsed verbatim): control_run/results_iq.txt (analyze_iq.py), control_run/results_pairs.txt
(analyze_pairs_iq.py). Open markers: the committed analysis printed FINAL 'CANNOT_TELL (copies disagree)'; pooled value shown.
Sign: as printed by the analyses (idle minus readout), no folding."""
import os, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as fs

fs.setup()
CR = os.path.join(fs.ROOT, 'control_run'); num = r'([-+]?[\d.]+)'
def sections(path):
    out, chip = {}, None
    for line in open(path):
        m = re.match(r'^(ibm_\w+)', line)
        if m: chip = m.group(1)
        if chip: out.setdefault(chip, []).append(line)
    return {c: ''.join(v) for c, v in out.items()}
iq, pr = sections(os.path.join(CR, 'results_iq.txt')), sections(os.path.join(CR, 'results_pairs.txt'))
def get(txt, name):
    m = re.search(r'%s\s+pooled\s+%s \[%s, %s\] n=(\d+)' % (re.escape(name), num, num, num), txt)
    f = re.search(r'%s\s+FINAL -> (.*)' % re.escape(name), txt).group(1)
    return float(m.group(1)), float(m.group(2)), float(m.group(3)), 'disagree' in f, int(m.group(4)), f.strip()

SEL = (('all shots', fs.BLACK, 'o', -0.18), ('clean shots', fs.GREEN, 's', 0.18))
fig, axs = plt.subplots(2, 3, figsize=(fs.COL2, 3.6), sharey='row')
rows = []
for j, chip in enumerate(fs.CHIPS):
    for r, (Ns, names, src, ylab) in enumerate((
            ((2, 4, 8), lambda N: ('N=%d R2 all-shot loss' % N, 'N=%d R1 clean-shot loss' % N), iq, 'coherence cost\n(idle $-$ readout)'),
            ((2, 4), lambda N: ('R3 N=%d all-shot concurrence loss (K=39)' % N, 'R1 N=%d clean-shot concurrence loss (K=39)' % N), pr,
             'concurrence cost\n(idle $-$ readout)'))):
        ax = axs[r, j]
        for N in Ns:
            for (lab, col, mk, dx), name in zip(SEL, names(N)):
                m, lo, hi, dis, n, fin = get(src[chip], name)
                ax.errorbar(np.log2(N) + dx, m, yerr=[[m - lo], [hi - m]], fmt=mk, color=col, mfc='white' if dis else col, mew=0.8, ms=3.5, elinewidth=0.8)
                rows.append((chip, name, m, lo, hi, n, fin))
        ax.axhline(0, color=fs.GREY, lw=0.5, zorder=0)
        ax.set_xticks([np.log2(N) for N in Ns]); ax.set_xticklabels([str(N) for N in Ns]); ax.set_xlim(0.4, np.log2(Ns[-1]) + 0.6)
        ax.set_xlabel('number of readouts $N$')
        ax.text(0.04, 0.95, fs.SHORT[chip], transform=ax.transAxes, ha='left', va='top', fontsize=7.5)
        if j == 0: ax.set_ylabel(ylab)
        fs.panel(ax, 'abcdef'[3 * r + j])
h = [Line2D([], [], marker=mk, color=col, ls='none', ms=3.5) for _, col, mk, _ in SEL] + [Line2D([], [], marker='o', color='k', mfc='white', ls='none', ms=3.5)]
axs[0, 0].legend(h, [s[0] for s in SEL] + ['copies disagree'], loc='upper right', fontsize=6.5, handletextpad=0.2)
fig.tight_layout(h_pad=0.8, w_pad=0.8)
for r_ in rows: print('  %-14s %-44s %+.4f [%+.4f, %+.4f] n=%d %s' % r_)
fs.save(fig, 4)
