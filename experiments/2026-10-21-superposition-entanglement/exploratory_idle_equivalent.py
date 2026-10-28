"""EXPLORATORY (written after the data was seen): idle-equivalent steps per readout.
For each qubit (superposition) or pair (entanglement, unclipped witness that keeps its sign), the idle curve D(0,8,16,32)
is interpolated to find how many idle steps lose as much as N readouts did. Ratio = idle steps / N.
N=32 is censored (values beyond the idle curve's range are dropped) and is printed but not read."""
import json, numpy as np
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
rng = np.random.default_rng(7)
def boot(x, n=2000):
    x = np.asarray(x); b = [np.median(rng.choice(x, len(x))) for _ in range(n)]; return np.median(x), np.percentile(b, [2.5, 97.5])
def neff(dN, dv, v):
    dv = np.maximum.accumulate(dv[::-1])[::-1]
    if v > dv[0] or v < dv[-1]: return np.nan
    for i in range(len(dN) - 1):
        a, b = dv[i], dv[i + 1]
        if a >= v >= b: return dN[i] if a == b else dN[i] + (dN[i + 1] - dN[i]) * (a - v) / (a - b)
    return np.nan
CONDS = [('D', 0), ('D', 8), ('D', 16), ('D', 32), ('M0', 8), ('M0', 16), ('M0', 32)]
for chip in ['ibm_fez', 'ibm_kingston', 'ibm_marrakesh']:
    rec = json.load(open('job_%s.json' % chip)); z = np.load('raw_%s.npz' % chip); P = rec['pairs']; SG = rec['singles']; nP = len(P); meta = rec['meta']
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
    def W(c):
        p = c['ZZ']; r03 = abs(complex(c['XX'] - c['YY'], -(c['XY'] + c['YX']))) / 4; return 2 * (r03 - np.sqrt(p[1] * p[2]))
    Ent = {k: np.array([W(pc[k + (j,)]) for j in range(nP)]) for k in CONDS}
    Coh = {k: np.array([abs(complex(np.mean(sc[k + (j,)]['X']), np.mean(sc[k + (j,)]['Y']))) for j in range(len(SG))]) for k in CONDS}
    dN = np.array([0, 8, 16, 32.])
    for name, V, nu in (('superposition', Coh, len(SG)), ('entanglement, unclipped', Ent, nP)):
        out = []
        for N in (8, 16, 32):
            r = [x for x in (neff(dN, np.array([V[('D', n)][j] for n in (0, 8, 16, 32)]), V[('M0', N)][j]) / N for j in range(nu)) if np.isfinite(x)]
            if len(r) >= 5: m, ci = boot(r); out.append('N=%d: %.3f [%.3f, %.3f] n=%d' % (N, m, ci[0], ci[1], len(r)))
            else: out.append('N=%d: CANNOT_TELL n=%d' % (N, len(r)))
        print('%-14s %-24s ' % (chip, name) + ' | '.join(out))
