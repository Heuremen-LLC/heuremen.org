"""Analysis exactly as PREREG.md (hashed before submission). Prints every chip, with the declared cut and with no cut."""
import json, glob, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(7)
def p1(arr, nbits):
    nby = arr.shape[1]; out = np.zeros(nbits)
    for k in range(nbits):
        out[k] = ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).mean()
    return out
def slope(ns, ph):
    ph = np.unwrap(ph); A = np.vstack([ns, np.ones_like(ns)]).T; return np.linalg.lstsq(A, ph, rcond=None)[0][0]
def boot(x, f=np.median, n=2000):
    x = np.asarray(x); b = [f(rng.choice(x, len(x))) for _ in range(n)]; return f(x), np.percentile(b, [16, 84]), np.percentile(b, [2.5, 97.5])
def state(lo95, hi95, ref):
    return 'ABOVE' if lo95 > ref else 'BELOW' if hi95 < ref else 'CONSISTENT'
rows = []
for jf in sorted(glob.glob(HERE + '/job_*.json')):
    rec = json.load(open(jf)); chip = rec['chip']; raw = HERE + '/raw_%s.npz' % chip
    if not os.path.exists(raw): print(chip, 'NOTHING_RETURNED yet'); continue
    z = np.load(raw); meta = rec['meta']; S = rec['S']; L_us = rec['L_dt'] * 4e-3 * 1.0
    P = {}
    for i, m in enumerate(meta): P[(m['mode'], m['N'], m['quad'])] = p1(z['fin_%d' % i], len(S))
    NS = sorted({m['N'] for m in meta}); s = {}
    for mode in ('D', 'M0', 'M1'):
        ph = np.array([np.arctan2(2 * P[(mode, N, 'Y')] - 1, 2 * P[(mode, N, 'X')] - 1) for N in NS])  # (len NS, spectators)
        s[mode] = np.array([slope(np.array(NS, float), ph[:, k]) for k in range(len(S))])
    sD = np.degrees(s['D'])
    for label, keep in (('declared cut |s_D|>=1 deg/step', np.abs(sD) >= 1.0), ('no cut', np.ones(len(S), bool))):
        R0 = s['M0'][keep] / s['D'][keep]; R1 = s['M1'][keep] / s['D'][keep]
        m0, i68, i95 = boot(R0); m1, j68, j95 = boot(R1); d, k68, k95 = boot(R1 - R0)
        print('%-14s L=%.3f us  [%s] n=%d' % (chip, L_us, label, keep.sum()))
        print('   R0 = %.4f  68%% [%.4f, %.4f]  95%% [%.4f, %.4f]  -> %s vs 1   extra time per measurement dT = %+.1f ns' % (m0, *i68, *i95, state(i95[0], i95[1], 1), (m0 - 1) * L_us * 1e3))
        print('   R1 = %.4f  95%% [%.4f, %.4f]    R1-R0 = %+.4f  95%% [%+.4f, %+.4f] -> %s vs 0' % (m1, *j95, d, *k95, state(k95[0], k95[1], 0)))
        if label.startswith('declared'): rows.append((chip, L_us, m0, i95, (m0 - 1) * L_us * 1e3))
if len(rows) >= 2:
    print('\nAcross readout lengths (reading 3): chip, L, R0, dT ns')
    for r in rows: print('   %-14s %.3f us  R0 %.4f  dT %+.1f ns' % (r[0], r[1], r[2], r[4]))
