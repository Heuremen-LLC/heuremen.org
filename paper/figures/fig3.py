"""Fig. 3. Controls on the kick, per processor, median over untouched qubits with 95 % bootstrap intervals, as committed.
(a) kick per readout with 39, 10 and 1 of the 39 qubits read out (control run, N = 8, analyze.py R1 and R2, pooled copies), and the kick
    on the minority of untouched qubits whose idle phase turns the other way (control run and raw-readout run, check L 'forward' group,
    slope of the readout-minus-idle phase over N = 2, 4, 8).
(b) kick with the measured qubits in |1> minus kick with them in |0> (each against its own idle in the same state), deciding run
    (excitation) and the kingston repeat, side by side, never pooled.
Inputs (read only; the committed outputs of the analyses, parsed verbatim):
  control_run/results_prereg.txt   (analyze.py)                 R1 kick pooled; R2 kick with K=10 / K=1 measured pooled; 'copies disagree' flags
  control_run/results_check_l.txt  (check_l_backward.py)        forward group, control run and raw-readout run
  deciding_runs/results_excitation.txt (analyze_excitation.py)  RJ.3 per chip
  excitation_kingston/results_exk.txt  (analyze_excitation.py on the repeat)  RJ.3 kingston
SIGN CONVENTION: carried exactly as the analyses print it. All four analyses read the Y axis with the same convention, under which the
kick is NEGATIVE (the paper text quotes kick sizes with the first run's positive sign). In (b) a POSITIVE value therefore means the kick
with the measured qubits in |1> is SMALLER in size than in |0>. Open markers: the committed analysis printed 'copies disagree'
(CANNOT TELL); the pooled value and interval are shown."""
import os, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as fs

fs.setup()
R = fs.ROOT; num = r'([-+]?[\d.]+(?:e[-+]?\d+)?)'

def sections(path, head=r'^(?:\S+\s+)?(ibm_\w+)'):
    out, chip = {}, None
    for line in open(path):
        m = re.match(r'^(?:excitation )?(ibm_\w+)', line)
        if m: chip = m.group(1)
        if chip: out.setdefault(chip, []).append(line)
    return out

# control run, analyze.py
pre = sections(os.path.join(R, 'control_run/results_prereg.txt'))
def prereg(chip, name):
    txt = ''.join(pre[chip]); m = re.search(r'\s%s\s+pooled\s+%s deg \[%s, %s\]' % (re.escape(name), num, num, num), txt)
    dis = re.search(r'\s%s\s+copies disagree' % re.escape(name), txt) is not None
    return float(m.group(1)), float(m.group(2)), float(m.group(3)), dis
# check L
chl = sections(os.path.join(R, 'control_run/results_check_l.txt'))
def checkl(chip, run):
    txt = ''.join(chl[chip]); m = re.search(r'%s\s+forward \(n=(\d+)\) %s \[%s, %s\]' % (re.escape(run), num, num, num), txt)
    return float(m.group(2)), float(m.group(3)), float(m.group(4)), False, int(m.group(1))
# excitation RJ.3
def rj3(path):
    out = {}
    for chip, lines in sections(path).items():
        for line in lines:
            m = re.search(r'RJ\.3 G = kick in 1 minus kick in 0\s+%s \[%s, %s\] deg/step' % (num, num, num), line)
            if m: out[chip] = tuple(float(v) for v in m.groups()) + (False,)
    return out
exc = rj3(os.path.join(R, 'deciding_runs/results_excitation.txt')); exk = rj3(os.path.join(R, 'excitation_kingston/results_exk.txt'))

series_a = [('39 read out', fs.VERM, 'o', lambda c: prereg(c, 'R1 kick')),
            ('10 read out', fs.ORANGE, 's', lambda c: prereg(c, 'R2 kick with K=10 measured')),
            ('1 read out', fs.SKY, 'D', lambda c: prereg(c, 'R2 kick with K=1 measured')),
            ('minority, control run', fs.GREEN, '^', lambda c: checkl(c, 'control run')),
            ('minority, raw-readout run', fs.PURPLE, 'v', lambda c: checkl(c, 'raw-readout run'))]

fig, (ax, bx) = plt.subplots(2, 1, figsize=(fs.COL1, 3.7), gridspec_kw=dict(height_ratios=[1.55, 1]))
xs = np.arange(3); off = np.linspace(-0.3, 0.3, len(series_a)); printed = []
for (lab, col, mk, get), dx in zip(series_a, off):
    for i, chip in enumerate(fs.CHIPS):
        v = get(chip); m, lo, hi, dis = v[:4]
        ax.errorbar(xs[i] + dx, m, yerr=[[m - lo], [hi - m]], fmt=mk, color=col, mfc='white' if dis else col, mew=0.8, ms=3.5, elinewidth=0.8)
        printed.append(('a', lab, chip, m, lo, hi, 'copies disagree' if dis else ''))
ax.axhline(0, color=fs.GREY, lw=0.5, zorder=0)
ax.set_xticks(xs); ax.set_xticklabels([fs.SHORT[c] for c in fs.CHIPS]); ax.set_xlim(-0.5, 2.5)
ax.set_ylabel('kick per readout (deg)')
h = [Line2D([], [], marker=mk, color=col, ls='none', ms=3.5) for _, col, mk, _ in series_a]
h.append(Line2D([], [], marker='o', color='k', mfc='white', ls='none', ms=3.5))
ax.legend(h, [s[0] for s in series_a] + ['open: copies disagree'], loc='lower left', ncol=2, fontsize=6.5, columnspacing=0.8, handletextpad=0.2)
ax.set_ylim(-6.6, 0.6)
fs.panel(ax, 'a')

for i, chip in enumerate(fs.CHIPS):
    for dx, src, col, mk, lab in ((-0.12, exc, fs.BLUE, 'o', 'deciding run'), (0.12, exk, fs.BLACK, 's', 'kingston repeat')):
        if chip not in src: continue
        m, lo, hi, dis = src[chip]
        bx.errorbar(xs[i] + (dx if chip in exk else 0), m, yerr=[[m - lo], [hi - m]], fmt=mk, color=col, mfc=col, ms=3.5, elinewidth=0.8)
        printed.append(('b', lab, chip, m, lo, hi, ''))
bx.axhline(0, color=fs.GREY, lw=0.5, zorder=0)
bx.set_xticks(xs); bx.set_xticklabels([fs.SHORT[c] for c in fs.CHIPS]); bx.set_xlim(-0.5, 2.5)
bx.set_ylabel('kick in $|1\\rangle$ $-$ kick in $|0\\rangle$\n(deg per readout)')
bx.legend([Line2D([], [], marker='o', color=fs.BLUE, ls='none', ms=3.5), Line2D([], [], marker='s', color=fs.BLACK, ls='none', ms=3.5)],
          ['deciding run', 'kingston repeat'], loc='upper left', fontsize=6.5)
fs.panel(bx, 'b')
fig.tight_layout(h_pad=0.6)
for p in printed: print('  %s %-26s %-14s %+.4f [%+.4f, %+.4f] %s' % p)
fs.save(fig, 3)
