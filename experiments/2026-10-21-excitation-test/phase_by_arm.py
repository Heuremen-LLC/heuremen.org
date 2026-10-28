"""Phase per step by arm (analysis after the data, at an outside referee's request): least-squares slope of each untouched
qubit's phase over N = 0, 4, 8, in degrees per step; median of arm differences over qubits, 95 % bootstrap."""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(5)
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def p1(arr, n): return np.array([bits(arr, k).mean() for k in range(n)])
def med(x): b = [np.median(rng.choice(x, len(x))) for _ in range(2000)]; return np.median(x), np.percentile(b, [2.5, 97.5])
for chip in ['ibm_kingston', 'ibm_marrakesh']:
    rec = json.load(open(os.path.join(HERE, 'job_%s.json' % chip))); z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); S = rec['S']; meta = rec['meta']
    P = {(m['mode'], m['N'], m['quad']): p1(z['fin_%d' % i], len(S)) for i, m in enumerate(meta)}; s = {}
    for mo in ('D0', 'M0', 'D1', 'M1'):
        ph = np.array([np.arctan2(2 * P[(mo, N, 'Y')] - 1, 2 * P[(mo, N, 'X')] - 1) for N in (0, 4, 8)])
        s[mo] = np.degrees(np.array([np.polyfit([0, 4, 8], np.unwrap(ph[:, k]), 1)[0] for k in range(len(S))]))
    for lab, a, b in (('readout kick, measured qubits in 0 (M0-D0)', 'M0', 'D0'), ('readout kick, measured qubits in 1 (M1-D1)', 'M1', 'D1'),
                      ('measured qubits in 1, no readout (D1-D0)', 'D1', 'D0'), ('M1-D0 (mixes the two)', 'M1', 'D0')):
        m, ci = med(s[a] - s[b]); print('%-14s %-44s %+.2f deg/step [%+.2f, %+.2f]' % (chip, lab, m, *ci))
