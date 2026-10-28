"""Analysis exactly as PREREG.md (hashed before submission). Declared cut and no cut."""
import json, glob, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(7)
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def boot(x, f=np.median, n=2000):
    x = np.asarray(x); b = [f(rng.choice(x, len(x))) for _ in range(n)]; return f(x), np.percentile(b, [16, 84]), np.percentile(b, [2.5, 97.5])
def state(lo, hi, ref): return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'
def conc(c):
    p = c['ZZ']; XX, YY, XY, YX = c['XX'], c['YY'], c['XY'], c['YX']
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))
CONDS = [('D', 0), ('D', 8), ('D', 16), ('D', 32), ('M0', 8), ('M0', 16), ('M0', 32)]
summary = []
for jf in sorted(glob.glob(HERE + '/job_*.json')):
    rec = json.load(open(jf)); chip = rec['chip']; raw = HERE + '/raw_%s.npz' % chip
    if not os.path.exists(raw): print(chip, 'NOTHING_RETURNED yet'); continue
    z = np.load(raw); P = rec['pairs']; SG = rec['singles']; meta = rec['meta']; nP = len(P)
    pc, sc = {}, {}
    for i, m in enumerate(meta):
        arr = z['fin_%d' % i]; key = (m['mode'], m['N']); st = m['setting']
        for j in range(nP):
            a = bits(arr, 2 * j); c = bits(arr, 2 * j + 1); cell = pc.setdefault(key + (j,), {})
            if st == 'ZZ': cell['ZZ'] = [np.mean((a == x) & (c == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))]
            else: cell[st] = np.mean((1 - 2 * a) * (1 - 2 * c))
        for j in range(len(SG)):
            s = 1 - 2 * bits(arr, 2 * nP + j); cell = sc.setdefault(key + (j,), {'X': [], 'Y': []})
            if st[0] in 'XY': cell[st[0]].append(s.mean())
    C = {k: np.array([conc(pc[k + (j,)]) for j in range(nP)]) for k in CONDS}
    Co = {k: np.array([abs(complex(np.mean(sc[k + (j,)]['X']), np.mean(sc[k + (j,)]['Y']))) for j in range(len(SG))]) for k in CONDS}
    print('%s  L=%.3f us  pairs %d  singles %d  shots %d' % (chip, rec['L_dt'] * 4e-3, nP, len(SG), rec['shots']))
    for name, V in (('ENTANGLEMENT (concurrence)', C), ('SUPERPOSITION (coherence)', Co)):
        print('  %s  median: D0 %.3f | D8 %.3f D16 %.3f D32 %.3f | M8 %.3f M16 %.3f M32 %.3f' % ((name,) + tuple(np.median(V[k]) for k in CONDS)))
        v0 = V[('D', 0)]
        for N in (8, 16, 32):
            LD = v0 - V[('D', N)]; LM = v0 - V[('M0', N)]
            for label, keep in (('declared', (v0 >= 0.5) & (LD >= 0.02)), ('no cut', LD > 0)):
                if keep.sum() < 3: print('    N=%2d [%s] n=%d CANNOT_TELL' % (N, label, keep.sum())); continue
                q, i68, i95 = boot(LM[keep] / LD[keep])
                tag = '  <- PRIMARY' if (N == 32 and label == 'declared') else ''
                print('    N=%2d [%-8s] Q = %.3f  68%% [%.3f, %.3f]  95%% [%.3f, %.3f]  n=%d -> %s vs 1%s' % (N, label, q, *i68, *i95, keep.sum(), state(i95[0], i95[1], 1), tag))
                if tag: summary.append((chip, name.split()[0], q, i95, state(i95[0], i95[1], 1)))
        d, _, d95 = boot(V[('M0', 32)] - V[('D', 32)])
        print('    value_M(32) - value_D(32) = %+.4f  95%% [%+.4f, %+.4f] -> %s vs 0' % (d, *d95, state(d95[0], d95[1], 0)))
print('\nPRIMARY (N=32, declared cut):')
for r in summary: print('  %-14s %-14s Q=%.3f 95%% [%.3f, %.3f] %s' % (r[0], r[1], r[2], *r[3], r[4]))
