"""Analysis exactly as PREREG.md (hashed before submission). Prints the declared cut and no cut."""
import json, glob, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(7)
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def boot(x, f=np.median, n=2000):
    x = np.asarray(x); b = [f(rng.choice(x, len(x))) for _ in range(n)]; return f(x), np.percentile(b, [16, 84]), np.percentile(b, [2.5, 97.5])
def state(lo, hi, ref): return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'
def conc(cell):
    p = cell['ZZ']; XX, YY, XY, YX = (cell[s] for s in ('XX', 'YY', 'XY', 'YX'))
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))
for jf in sorted(glob.glob(HERE + '/job_*.json')):
    rec = json.load(open(jf)); chip = rec['chip']; raw = HERE + '/raw_%s.npz' % chip
    if not os.path.exists(raw): print(chip, 'NOTHING_RETURNED yet'); continue
    z = np.load(raw); P = rec['pairs']; E = set(rec['E']); meta = rec['meta']
    C = {}  # (mode,N) -> array over pairs
    cells = {}
    for i, m in enumerate(meta):
        arr = z['fin_%d' % i]
        for j in range(len(P)):
            a = bits(arr, 2 * j); c = bits(arr, 2 * j + 1); key = (m['mode'], m['N'], j)
            cell = cells.setdefault(key, {})
            if m['setting'] == 'ZZ': cell['ZZ'] = [np.mean((a == x) & (c == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))]
            else: cell[m['setting']] = np.mean((1 - 2 * a) * (1 - 2 * c))
    for (mode, N, j), cell in cells.items(): C.setdefault((mode, N), np.zeros(len(P)))[j] = conc(cell)
    C0 = C[('D', 0)]; kD = (C0 - C[('D', 16)]) / 16; kM = (C0 - C[('M0', 16)]) / 16
    ns = np.array([0, 8, 16.])
    sl = lambda y: np.polyfit(ns, y, 1)[0]
    sD = np.array([-sl([C0[j], C[('D', 8)][j], C[('D', 16)][j]]) for j in range(len(P))])
    sM = np.array([-sl([C0[j], C[('M0', 8)][j], C[('M0', 16)][j]]) for j in range(len(P))])
    L_us = rec['L_dt'] * 4e-3
    print('%s  L=%.3f us  pairs %d  shots %d' % (chip, L_us, len(P), rec['shots']))
    print('   median C: D0 %.3f  D8 %.3f  D16 %.3f  M8 %.3f  M16 %.3f' % tuple(np.median(C[k]) for k in [('D', 0), ('D', 8), ('D', 16), ('M0', 8), ('M0', 16)]))
    for label, keep in (('declared cut C0>=0.5, k_D>=0.002/step', (C0 >= 0.5) & (kD >= 0.002)), ('no cut (k_D>0 for ratio)', kD > 0)):
        for nm, num, den in (('endpoint', kM, kD), ('slope', sM, sD)):
            k2 = keep & (den > 0)
            if k2.sum() < 3: print('   [%s] %s: n=%d CANNOT_TELL' % (label, nm, k2.sum())); continue
            q, i68, i95 = boot(num[k2] / den[k2])
            print('   [%s] Q_%s = %.3f  68%% [%.3f, %.3f]  95%% [%.3f, %.3f]  n=%d  -> %s vs 1, %s vs 1.2' % (label, nm, q, *i68, *i95, k2.sum(), state(i95[0], i95[1], 1), state(i95[0], i95[1], 1.2)))
    d, j68, j95 = boot(C[('M0', 16)] - C[('D', 16)])
    print('   C_M0(16) - C_D(16) = %+.4f  95%% [%+.4f, %+.4f]  -> %s vs 0' % (d, *j95, state(j95[0], j95[1], 0)))
