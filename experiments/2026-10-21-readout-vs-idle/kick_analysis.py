"""Additive-kick analysis of the 21 October readout-vs-idle run (written AFTER the data were seen, at an outside referee's request).
Per spectator: excess phase per step = s_M0 - s_D. Fit excess = a * s_D + c. Extra elapsed time predicts a ~ +0.19, c = 0;
a fixed phase kick per readout predicts a = 0, c != 0. Also: median excess on spectators precessing backwards (s_D < 0),
and R0 with the |rate| >= 4 deg/us conditioning cut used in the October campaign. Bootstrap over spectators."""
import json, glob, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(1)
def p1(arr, nbits):
    nby = arr.shape[1]; return np.array([((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).mean() for k in range(nbits)])
def slope(ns, ph): return np.polyfit(ns, np.unwrap(ph), 1)[0]
for jf in sorted(glob.glob(HERE + '/job_*.json')):
    rec = json.load(open(jf)); chip = rec['chip']; raw = HERE + '/raw_%s.npz' % chip
    if not os.path.exists(raw): print(chip, 'NOTHING_RETURNED'); continue
    z = np.load(raw); S = rec['S']; meta = rec['meta']; L = rec['L_dt'] * 4e-3
    P = {(m['mode'], m['N'], m['quad']): p1(z['fin_%d' % i], len(S)) for i, m in enumerate(meta)}
    NS = sorted({m['N'] for m in meta}); s = {}
    for mode in ('D', 'M0', 'M1'):
        ph = np.array([np.arctan2(2 * P[(mode, N, 'Y')] - 1, 2 * P[(mode, N, 'X')] - 1) for N in NS])
        s[mode] = np.degrees(np.array([slope(np.array(NS, float), ph[:, k]) for k in range(len(S))]))
    x = s['D']; d = s['M0'] - x; d1 = s['M1'] - x
    a, c = np.polyfit(x, d, 1)
    bs = np.array([np.polyfit(x[i], d[i], 1) for i in [rng.integers(0, len(x), len(x)) for _ in range(2000)]])
    ca, cc = np.percentile(bs[:, 0], [2.5, 97.5]), np.percentile(bs[:, 1], [2.5, 97.5])
    kick0 = np.median(d); kb = np.percentile([np.median(rng.choice(d, len(d))) for _ in range(2000)], [2.5, 97.5])
    kick1 = np.median(d1); kb1 = np.percentile([np.median(rng.choice(d1, len(d1))) for _ in range(2000)], [2.5, 97.5])
    neg = x < 0; keep = np.abs(x) / L >= 4
    print('%-14s L=%.3f us  n=%d' % (chip, L, len(x)))
    print('   excess = a*s_D + c :  a = %+.3f [%+.3f, %+.3f]   c = %+.2f [%+.2f, %+.2f] deg/readout' % (a, *ca, c, *cc))
    print('   median kick per readout, measured qubits in 0: %+.2f [%+.2f, %+.2f] deg;  in 1: %+.2f [%+.2f, %+.2f] deg' % (kick0, *kb, kick1, *kb1))
    print('   spectators precessing backwards (n=%d): median excess %+.2f deg/readout' % (neg.sum(), np.median(d[neg])))
    print('   R0 with |rate| >= 4 deg/us (n=%d): %.3f' % (keep.sum(), np.median(s['M0'][keep] / x[keep])))
