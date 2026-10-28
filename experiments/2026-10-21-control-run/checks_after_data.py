"""Checks A and B exactly as CHECKS_AFTER_DATA.md (defined after the control-run results were seen, before these were computed)."""
import json, glob, os, sys
import numpy as np
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
B = 2000
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi
def bm(x, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 3: return np.nan, (np.nan, np.nan), len(x)
    b = np.median(rng.choice(x, (B, len(x))), axis=1); return float(np.median(x)), tuple(np.percentile(b, [2.5, 97.5])), len(x)
def state(ci, ref=0):
    lo, hi = ci
    if not np.isfinite(lo): return 'CANNOT_TELL'
    return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'
def conc(c):
    p = c['ZZ']; XX, YY, XY, YX = c['XX'], c['YY'], c['XY'], c['YX']
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))

for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); meta = rec['meta']
    props = json.load(open(os.path.join(HERE, 'backend_properties_%s.json' % chip)))
    edges = {tuple(sorted(g['qubits'])) for g in props['gates'] if g['gate'] == 'cz'}
    adj = {}
    for a, b in edges: adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    E = rec['E']; S = rec['S']; P = rec['pairs']; SG = rec['singles']; nS = len(S); rng = np.random.default_rng(20261022)
    idx = lambda **kw: [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    print('=' * 90); print(chip)
    # ---------- Check A ----------
    def stuck(arm, N, K=39):
        """per measured qubit (order of the measured set): f = shots with two consecutive 1s, q = persistent part"""
        ii = idx(arm=arm, mode='M', N=N) if arm == 'P' else idx(arm=arm, mode='M', N=N, K=K)
        f, q = [], []
        regs = [[z['m%d_%d' % (k, i)] for k in range(N)] for i in ii]
        for e in range(len(E)):
            rows = [np.stack([bits(r[k], e) for k in range(N)], 1) for r in regs]   # shots x N per circuit
            R = np.concatenate(rows)
            f.append(np.mean(np.any(R[:, :-1] & R[:, 1:], axis=1)))
            q.append(np.mean([max(0.0, np.mean(R[:, k] & R[:, k + 1]) - np.mean(R[:, k]) * np.mean(R[:, k + 1])) for k in range(N - 1)]))
        return dict(zip(E, f)), dict(zip(E, q))
    def coh(mode, N):
        ex = lambda st: np.mean(np.concatenate([np.stack([1 - 2 * bits(z['fin_%d' % i], k) for k in range(nS)], 1) for i in idx(arm='S', mode=mode, N=N, setting=st)]), 0)
        return np.hypot(ex('XX'), ex('YY'))
    print('Check A, superposition (arm S):')
    for N in (2, 4, 8):
        f, q = stuck('S', N); loss = coh('D', N) - coh('M', N)
        nbE = [[e for e in adj.get(s, ()) if e in f] for s in S]
        Bc = np.array([2 * sum(f[e] for e in nb) for nb in nbE]); Be = np.array([2 * sum(q[e] for e in nb) for nb in nbE])
        has = np.array([len(nb) > 0 for nb in nbE])
        m1, c1, n1 = bm((loss - Bc)[has], rng); m2, c2, _ = bm((loss - Be)[has], rng); m0, c0, n0 = bm(loss[~has], rng)
        out = 'EXCLUDED' if state(c1) == 'ABOVE' else 'CANNOT TELL (excluded on the persistent estimate only)' if state(c2) == 'ABOVE' else 'NOT EXCLUDED'
        print('  N=%d  median stuck f %.4f q %.4f | wired qubits n=%d: loss-B %+.4f [%+.4f,%+.4f], loss-B_est %+.4f [%+.4f,%+.4f] -> %s | no measured neighbour n=%d: loss %+.4f [%+.4f,%+.4f] %s'
              % (N, np.median(list(f.values())), np.median(list(q.values())), n1, m1, *c1, m2, *c2, out, n0, m0, *c0, state(c0)))
    print('Check A, entanglement (arm P):')
    nP = len(P)
    def concs(mode, N):
        cells = [dict() for _ in range(nP)]
        for st in ('ZZ', 'XX', 'YY', 'XY', 'YX'):
            a = np.concatenate([z['fin_%d' % i] for i in idx(arm='P', mode=mode, N=N, setting=st)])
            for j in range(nP):
                b1 = bits(a, 2 * j); b2 = bits(a, 2 * j + 1)
                if st == 'ZZ': cells[j]['ZZ'] = [np.mean((b1 == x) & (b2 == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))]
                else: cells[j][st] = np.mean((1 - 2 * b1) * (1 - 2 * b2))
        return np.array([conc(c) for c in cells])
    for N in (2, 4):
        f, q = stuck('P', N); cD = concs('D', N); loss = cD - concs('M', N); keep = cD > 0
        nbE = [[e for e in (adj.get(a, set()) | adj.get(b, set())) if e in f] for a, b in P]
        Bc = np.array([2 * sum(f[e] for e in nb) for nb in nbE]); Be = np.array([2 * sum(q[e] for e in nb) for nb in nbE])
        has = np.array([len(nb) > 0 for nb in nbE]) & keep
        m1, c1, n1 = bm((loss - Bc)[has], rng); m2, c2, _ = bm((loss - Be)[has], rng); m0, c0, n0 = bm(loss[keep & ~has], rng)
        out = 'EXCLUDED' if state(c1) == 'ABOVE' else 'CANNOT TELL (excluded on the persistent estimate only)' if state(c2) == 'ABOVE' else 'NOT EXCLUDED'
        print('  N=%d  wired pairs n=%d: loss-B %+.4f [%+.4f,%+.4f], loss-B_est %+.4f [%+.4f,%+.4f] -> %s | pairs with no measured neighbour n=%d: loss %+.4f %s'
              % (N, n1, m1, *c1, m2, *c2, out, n0, m0, state(c0)))
    # ---------- Check B ----------
    print('Check B, where the kick lands (N = 8, pooled copies):')
    def phase(arm, mode, K):
        kw = dict(arm=arm, mode=mode, N=8) if arm == 'S' else dict(arm=arm, mode=mode, N=8, K=K)
        ex = lambda st: np.mean(np.concatenate([np.stack([1 - 2 * bits(z['fin_%d' % i], k) for k in range(nS)], 1) for i in idx(setting=st, **kw)]), 0)
        return np.arctan2(ex('YY'), ex('XX'))
    nbmed = {}
    for K, arm, Kset in ((1, 'C', rec['K1']), (10, 'C', rec['K10']), (39, 'S', E)):
        kick = np.degrees(wrap(phase(arm, 'M', K) - phase(arm, 'D', K))) / 8
        nb = np.array([len(adj.get(s, set()) & set(Kset)) > 0 for s in S])
        mN, cN, nN = bm(kick[nb], rng); mF, cF, nF = bm(kick[~nb], rng); nbmed[K] = mN
        out = 'NEIGHBOUR-LOCAL' if state(cN) != 'CONSISTENT' and state(cF) == 'CONSISTENT' else 'CHIP-WIDE' if state(cF) != 'CONSISTENT' else 'CANNOT TELL'
        print('  K=%2d  neighbours n=%3d: %+.3f [%+.3f,%+.3f] %s | others n=%3d: %+.3f [%+.3f,%+.3f] %s -> %s'
              % (K, nN, mN, *cN, state(cN), nF, mF, *cF, state(cF), out))
    print('  neighbour median kick, K=10 vs K=39: %+.3f vs %+.3f' % (nbmed[10], nbmed[39]))
