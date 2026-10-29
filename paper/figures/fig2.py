"""Fig. 2. Phase of the 117 untouched qubits against the number of steps N, idle arm and readout arm on the same axes,
median over qubits with 95 % bootstrap bands (2000 resamples), one panel per processor; control run (timestamped).
Inputs (read only): control_run/job_<chip>.json, control_run/raw_<chip>.npz, control_run/results_prereg.txt.
Quantities exactly as control_run/analyze.py: <X>, <Y> per untouched qubit = mean of (1 - 2 bit) of the final readout in the
XX / YY settings of arm S, both copies pooled; phase = arctan2(<Y>, <X>). The plotted phase of each qubit is taken relative to its
N = 0 phase (the N = 0 circuit is shared by both arms) and unwrapped along N. The kick printed in each panel is analyze.py's R1
(readout minus idle phase at N = 8, wrapped, divided by 8, median over qubits, pooled copies), recomputed here with analyze.py's
bootstrap stream (seed 20261021, copy 0, copy 1, pooled) and asserted equal to the committed results_prereg.txt.
SIGN CONVENTION: carried as in analyze.py, which reads the Y axis with the opposite convention to the first run, so the kick is
negative here; the paper text quotes its size (1.62, 2.86, 2.39 deg) with the first run's positive sign."""
import json, os, re
import numpy as np
import matplotlib.pyplot as plt
import figstyle as fs

fs.setup()
CR = os.path.join(fs.ROOT, 'control_run'); B = 2000
def bits(arr, k): nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi
def boot_median(x, rng):          # analyze.py's boot_median, unchanged
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    b = np.median(rng.choice(x, (B, len(x))), axis=1)
    return float(np.median(x)), tuple(np.percentile(b, [2.5, 97.5])), len(x)
committed = {}
chip = None
for line in open(os.path.join(CR, 'results_prereg.txt')):
    m = re.match(r'(ibm_\w+)\s+job', line)
    if m: chip = m.group(1)
    m = re.match(r'\s+R1 kick\s+pooled\s+([-+.\d]+) deg \[([-+.\d]+), ([-+.\d]+)\]', line)
    if m and chip: committed[chip] = tuple(float(v) for v in m.groups())

fig, axes = plt.subplots(1, 3, figsize=(fs.COL2, 2.15), sharey=True)
Ns = [0, 2, 4, 8]
for ax, chip, lab in zip(axes, fs.CHIPS, 'abc'):
    rec = json.load(open(os.path.join(CR, 'job_%s.json' % chip))); z = np.load(os.path.join(CR, 'raw_%s.npz' % chip))
    meta = rec['meta']; nS = len(rec['S'])
    def expv(mode, N, st, copies):
        ii = [i for i, m in enumerate(meta) if m['arm'] == 'S' and m['mode'] == mode and m['N'] == N and m['setting'] == st and m['copy'] in copies]
        a = np.concatenate([z['fin_%d' % i] for i in ii]); return np.array([(1 - 2 * bits(a, k)).mean() for k in range(nS)])
    def phase(mode, N, copies=(0, 1)): return np.arctan2(expv(mode, N, 'YY', copies), expv(mode, N, 'XX', copies))
    # R1 exactly as analyze.py (same rng stream: copy 0, copy 1, pooled)
    rng = np.random.default_rng(20261021); r1 = None
    for cps in ((0,), (1,), (0, 1)):
        k = np.degrees(wrap(phase('M', 8, cps) - phase('D', 8, cps))) / 8; r1 = boot_median(k, rng)
    m1, (lo1, hi1), _ = r1
    c = committed[chip]; assert np.allclose([m1, lo1, hi1], c, atol=6e-5), (chip, (m1, lo1, hi1), c)
    # curves
    p0 = phase('D', 0)
    rng = np.random.default_rng(1)
    for mode, col, name in (('D', fs.IDLE, 'idle arm'), ('M', fs.READ, 'readout arm')):
        P = np.stack([np.zeros(nS)] + [wrap(phase(mode, N) - p0) for N in Ns[1:]])
        P = np.degrees(np.unwrap(P, axis=0))
        med = np.median(P, axis=1); bs = np.array([np.median(P[:, rng.integers(0, nS, nS)], axis=1) for _ in range(B)])
        lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
        ax.fill_between(Ns, lo, hi, color=col, alpha=0.22, lw=0)
        ax.plot(Ns, med, 'o-', color=col, ms=3, label=name)
    ax.set_xticks(Ns); ax.set_xlabel('number of steps $N$')
    f = lambda v: ('%+.2f' % v).replace('-', '\u2212')
    ax.text(0.97, 0.96, '%s: kick %s$^\\circ$ [%s, %s]\nper readout' % (fs.SHORT[chip], f(m1), f(lo1), f(hi1)),
            transform=ax.transAxes, ha='right', va='top', fontsize=6.5)
    ax.axhline(0, color=fs.GREY, lw=0.4, zorder=0)
    fs.panel(ax, lab)
axes[0].set_ylim(-84, 27)
axes[0].set_ylabel(r'phase relative to $N=0$ (deg)')
axes[0].legend(loc='lower left')
fig.tight_layout(w_pad=1.0)
fs.save(fig, 2)
