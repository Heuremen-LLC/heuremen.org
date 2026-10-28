"""Check C exactly as CHECK_C_CLEAN_SHOTS.md."""
import json, glob, os, sys
import numpy as np
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)); B = 2000
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def bm(x, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 3: return np.nan, (np.nan, np.nan), len(x)
    b = np.median(rng.choice(x, (B, len(x))), axis=1); return float(np.median(x)), tuple(np.percentile(b, [2.5, 97.5])), len(x)
def state(ci):
    lo, hi = ci
    return 'CANNOT_TELL' if not np.isfinite(lo) else 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
def conc(c):
    p = c['ZZ']; XX, YY, XY, YX = c['XX'], c['YY'], c['XY'], c['YX']
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))
for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); meta = rec['meta']
    props = json.load(open(os.path.join(HERE, 'backend_properties_%s.json' % chip)))
    adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    E = rec['E']; Epos = {e: i for i, e in enumerate(E)}; S = rec['S']; P = rec['pairs']; nS, nP = len(S), len(P); rng = np.random.default_rng(20261023)
    idx = lambda **kw: [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    def any1(i, N):   # shots x measured qubits: 1 if that measured qubit read 1 on any of the N readouts
        return np.any(np.stack([np.stack([bits(z['m%d_%d' % (k, i)], e) for e in range(len(E))], 1) for k in range(N)], 0), axis=0)
    print('=' * 90); print(chip)
    for N in (2, 4, 8):
        nb = [[Epos[e] for e in adj.get(s, ()) if e in Epos] for s in S]
        def coh(mode, clean):
            ex = {}
            for st in ('XX', 'YY'):
                num = np.zeros(nS); den = np.zeros(nS)
                for i in idx(arm='S', mode=mode, N=N, setting=st):
                    a = z['fin_%d' % i]; vals = np.stack([1 - 2 * bits(a, k) for k in range(nS)], 1)
                    if clean:
                        A = any1(i, N); keep = np.stack([~A[:, n].any(axis=1) if n else np.ones(len(a), bool) for n in nb], 1)
                    else: keep = np.ones_like(vals, bool)
                    num += (vals * keep).sum(0); den += keep.sum(0)
                ex[st] = num / den
            return np.hypot(ex['XX'], ex['YY']), den
        cD, _ = coh('D', False); cMa, nall = coh('M', False); cMc, nkeep = coh('M', True)
        has = np.array([len(n) > 0 for n in nb])
        mA, cA, _ = bm((cD - cMa)[has], rng); mC, cC, n = bm((cD - cMc)[has], rng)
        print('  coherence N=%d: all shots loss %+.4f [%+.4f,%+.4f] %s | clean shots loss %+.4f [%+.4f,%+.4f] %s  n=%d, median shots kept %.1f%%'
              % (N, mA, *cA, state(cA), mC, *cC, state(cC), n, 100 * np.median((nkeep / nall)[has])))
    for N in (2, 4):
        nb = [[Epos[e] for e in (adj.get(a, set()) | adj.get(b, set())) if e in Epos] for a, b in P]
        def cc(mode, clean):
            cells = [dict() for _ in range(nP)]; kept = np.zeros(nP); tot = np.zeros(nP)
            for st in ('ZZ', 'XX', 'YY', 'XY', 'YX'):
                acc = [[] for _ in range(nP)]
                for i in idx(arm='P', mode=mode, N=N, setting=st):
                    a = z['fin_%d' % i]; A = any1(i, N) if clean else None
                    for j in range(nP):
                        keep = ~A[:, nb[j]].any(axis=1) if (clean and nb[j]) else np.ones(len(a), bool)
                        acc[j].append((bits(a, 2 * j)[keep], bits(a, 2 * j + 1)[keep])); kept[j] += keep.sum(); tot[j] += len(a)
                for j in range(nP):
                    b1 = np.concatenate([x for x, _ in acc[j]]); b2 = np.concatenate([y for _, y in acc[j]])
                    if st == 'ZZ': cells[j]['ZZ'] = [np.mean((b1 == x) & (b2 == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))]
                    else: cells[j][st] = np.mean((1 - 2 * b1) * (1 - 2 * b2))
            return np.array([conc(c) for c in cells]), kept / tot
        cD, _ = cc('D', False); cMa, _ = cc('M', False); cMc, frac = cc('M', True)
        keepP = (cD > 0) & np.array([len(n) > 0 for n in nb])
        mA, cA, _ = bm((cD - cMa)[keepP], rng); mC, cC, n = bm((cD - cMc)[keepP], rng)
        print('  concurrence N=%d: all shots loss %+.4f [%+.4f,%+.4f] %s | clean shots loss %+.4f [%+.4f,%+.4f] %s  n=%d, median shots kept %.1f%%'
              % (N, mA, *cA, state(cA), mC, *cC, state(cC), n, 100 * np.median(frac[keepP])))
