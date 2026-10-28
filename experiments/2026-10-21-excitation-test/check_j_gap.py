"""Check J exactly as CHECK_J_GAP.md.text."""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(20261033)
def bits(arr, k): nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def p1(arr, n): return np.array([bits(arr, k).mean() for k in range(n)])
def boot(x):
    b = np.median(rng.choice(x, (2000, len(x))), 1); lo, hi = np.percentile(b, [2.5, 97.5]); return np.median(x), lo, hi
def state(lo, hi): return 'ABOVE (bigger in 1)' if lo > 0 else 'BELOW (smaller in 1)' if hi < 0 else 'CONSISTENT (no difference resolved)'
rows = []
for chip in ['ibm_kingston', 'ibm_marrakesh']:
    rec = json.load(open(os.path.join(HERE, 'job_%s.json' % chip))); z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); S = rec['S']; meta = rec['meta']
    P = {(m['mode'], m['N'], m['quad']): p1(z['fin_%d' % i], len(S)) for i, m in enumerate(meta)}; s = {}
    for mo in ('D0', 'M0', 'D1', 'M1'):
        ph = np.array([np.arctan2(2 * P[(mo, N, 'Y')] - 1, 2 * P[(mo, N, 'X')] - 1) for N in (0, 4, 8)])
        s[mo] = np.degrees(np.array([np.polyfit([0, 4, 8], np.unwrap(ph[:, k]), 1)[0] for k in range(len(S))]))
    G = (s['M1'] - s['D1']) - (s['M0'] - s['D0']); m, lo, hi = boot(G); k0 = np.median(s['M0'] - s['D0'])
    print('%-14s G = kick(in 1) - kick(in 0): %+.3f deg/step [%+.3f, %+.3f] n=%d -> %s | kick in 0 %+.3f; G / kick(0) = %+.2f [%+.2f, %+.2f]'
          % (chip, m, lo, hi, len(G), state(lo, hi), k0, m / k0, lo / k0, hi / k0)); rows.append((m, lo, hi))
v = np.array([r[0] for r in rows]); sd = np.array([(r[2] - r[1]) / 3.92 for r in rows]); w = 1 / sd ** 2
mc = np.sum(w * v) / np.sum(w); se = 1 / np.sqrt(np.sum(w)); chi = np.sum(((v - mc) / sd) ** 2) / (len(v) - 1)
print('combined G %+.3f [%+.3f, %+.3f] -> %s | chip-to-chip spread chi2/dof %.2f' % (mc, mc - 1.96 * se, mc + 1.96 * se, state(mc - 1.96 * se, mc + 1.96 * se), chi))
