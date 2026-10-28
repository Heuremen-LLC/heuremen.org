"""Plain concurrence difference (the same ruler fix as coherence_difference.py): for each Bell pair, concurrence after N readouts
minus after N equal idles; median over pairs, 95 % bootstrap over pairs. Reads job_*.json and raw_*.npz. Written after the data."""
import json, os, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(11)
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def conc(c):
    p = c['ZZ']; XX, YY, XY, YX = c['XX'], c['YY'], c['XY'], c['YX']
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))
def med(x): b = [np.median(rng.choice(x, len(x))) for _ in range(2000)]; return np.median(x), np.percentile(b, [2.5, 97.5])
for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip))
    nP = len(rec['pairs']); pc = {}
    for i, m in enumerate(rec['meta']):
        arr = z['fin_%d' % i]; st = m['setting']
        for j in range(nP):
            a = bits(arr, 2 * j); c = bits(arr, 2 * j + 1); cell = pc.setdefault((m['mode'], m['N'], j), {})
            if st == 'ZZ': cell['ZZ'] = [np.mean((a == x) & (c == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))]
            else: cell[st] = np.mean((1 - 2 * a) * (1 - 2 * c))
    C = lambda mo, N: np.array([conc(pc[(mo, N, j)]) for j in range(nP)])
    out = []
    for N in (8, 16):
        d = C('M0', N) - C('D', N); m, ci = med(d)
        out.append('N=%2d %+.4f [%+.4f, %+.4f]  idle median %.3f' % (N, m, *ci, np.median(C('D', N))))
    print('%-14s readout minus idle concurrence (median over %d pairs, 95%%): ' % (chip, nP) + ' | '.join(out))
