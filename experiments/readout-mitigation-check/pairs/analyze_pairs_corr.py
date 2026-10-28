"""Analysis exactly as PREREG_PAIRS.md. Reads pairsjob_<chip>_copy<c>.json, pairsraw_<chip>_copy<c>.npz and pairs_backend_properties_<chip>.json
from the folder given as the first argument (default: this folder)."""
import json, glob, os, sys
import numpy as np
HERE = sys.argv[1]; MULT = float(sys.argv[2]); B = 2000; D2 = 3.5 ** 2
SETS = ('ZZ', 'XX', 'YY', 'XY', 'YX')
def bm(x, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 3: return np.nan, (np.nan, np.nan), len(x)
    b = np.median(rng.choice(x, (B, len(x))), axis=1); return float(np.median(x)), tuple(np.percentile(b, [2.5, 97.5])), len(x)
def state(ci):
    lo, hi = ci
    return 'CANNOT_TELL' if not np.isfinite(lo) else 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
def maha2(x, mu, cov):
    d = np.stack([x.real - mu.real, x.imag - mu.imag], -1); return np.einsum('...i,ij,...j->...', d, np.linalg.inv(cov), d)
def cov2(x): return np.cov(np.stack([x.real, x.imag])) + 1e-12 * np.eye(2)
def conc(c):
    p = c['ZZ']; XX, YY, XY, YX = c['XX'], c['YY'], c['XY'], c['YX']
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))
chips = sorted({os.path.basename(f).split('_copy')[0][len('pairsjob_'):] for f in glob.glob(os.path.join(HERE, 'pairsjob_*_copy*.json'))})
for chip in chips:
    copies = {}
    for c in (0, 1):
        jf = os.path.join(HERE, 'pairsjob_%s_copy%d.json' % (chip, c)); raw = os.path.join(HERE, 'pairsraw_%s_copy%d.npz' % (chip, c))
        if os.path.exists(jf) and os.path.exists(raw): copies[c] = (json.load(open(jf)), np.load(raw))
    if not copies: print(chip, 'NOTHING_RETURNED'); continue
    rec0 = next(iter(copies.values()))[0]; E = rec0['E']; P = rec0['pairs']; nP = len(P); n = len(E) + len(rec0['S']); K1 = rec0['K1']
    props = json.load(open(os.path.join(HERE, 'pairs_backend_properties_%s.json' % chip)))
    adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    Epos = {e: i for i, e in enumerate(E)}
    nb39 = [[Epos[e] for e in (adj.get(a, set()) | adj.get(b, set())) if e in Epos] for a, b in P]
    rng = np.random.default_rng(20261028)
    cal = {st: np.concatenate([z['c_%d' % i] for (r, z) in copies.values() for i, m in enumerate(r['meta']) if m['arm'] == 'CAL' and m['state'] == st]) for st in '01'}
    mu0 = cal['0'].mean(0); mu1 = cal['1'].mean(0); u = mu1 - mu0
    proj = lambda x, q: np.real((x - mu0[q]) * np.conj(u[q])) / np.abs(u[q]) ** 2
    cov0 = [cov2(cal['0'][:, q][proj(cal['0'][:, q], q) < 0.5]) for q in range(n)]; cov1 = [cov2(cal['1'][:, q][proj(cal['1'][:, q], q) > 0.5]) for q in range(n)]
    E01 = MULT * np.array([(proj(cal['0'][:, q], q) > 0.5).mean() for q in range(n)]); E10 = MULT * np.array([(proj(cal['1'][:, q], q) < 0.5).mean() for q in range(n)])
    print('INJOB_CAL', chip, 'median over pair qubits e01 %.4f e10 %.4f (x1)' % (np.median([E01[q] for p in P for q in p]) / max(MULT,1e-9), np.median([E10[q] for p in P for q in p]) / max(MULT,1e-9)))
    print('=' * 96); print(chip, 'copies returned:', sorted(copies))
    def concs(cps, mode, N, K, clean):
        cells = [dict() for _ in range(nP)]; kept = np.zeros(nP); tot = np.zeros(nP)
        for st in SETS:
            acc = [[] for _ in range(nP)]
            for c in cps:
                if c not in copies: continue
                r, z = copies[c]
                for i, m in enumerate(r['meta']):
                    if not (m['arm'] == 'P' and m['mode'] == mode and m['N'] == N and m['setting'] == st and (mode == 'D' or m['K'] == K)): continue
                    fin = z['fin_%d' % i]; nsh = fin.shape[0]
                    if clean:
                        Kset = E if K == 39 else K1; ok = np.ones((nsh, len(Kset)), bool)
                        for k in range(N):
                            mm = z['m%d_%d' % (k, i)]
                            for j, e in enumerate(Kset): ok[:, j] &= maha2(mm[:, j], mu0[e], cov0[e]) < D2
                        kpos = {e: j for j, e in enumerate(Kset)}
                    for j, (a, b) in enumerate(P):
                        if clean:
                            nn = [kpos[E[x]] for x in nb39[j] if E[x] in kpos]; keep = ok[:, nn].all(1) if nn else np.ones(nsh, bool)
                        else: keep = np.ones(nsh, bool)
                        b1 = (proj(fin[:, 2 * j], a) > 0.5).astype(int)[keep]; b2 = (proj(fin[:, 2 * j + 1], b) > 0.5).astype(int)[keep]
                        acc[j].append((b1, b2)); kept[j] += keep.sum(); tot[j] += nsh
            for j in range(nP):
                if not acc[j]: continue
                b1 = np.concatenate([x for x, _ in acc[j]]); b2 = np.concatenate([y for _, y in acc[j]])
                pm = np.array([np.mean((b1 == x) & (b2 == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))])
                (a, b) = P[j]
                A1 = np.array([[1 - E01[a], E10[a]], [E01[a], 1 - E10[a]]]); A2 = np.array([[1 - E01[b], E10[b]], [E01[b], 1 - E10[b]]])
                pc = np.linalg.solve(np.kron(A1, A2), pm)
                if st == 'ZZ': cells[j]['ZZ'] = list(pc)
                else: cells[j][st] = pc[0] - pc[1] - pc[2] + pc[3]
        return np.array([conc(c) if len(c) == 5 else np.nan for c in cells]), kept / np.maximum(tot, 1)
    CP = {'copy 0': [0], 'copy 1': [1], 'pooled': [0, 1]}
    def report(name, vals):
        sts = {lab: state(vals[lab][1]) for lab in CP}
        for lab in CP: m, ci, k = vals[lab]; print('  %-44s %-7s %+.4f [%+.4f, %+.4f] n=%d -> %s' % (name, lab, m, *ci, k, sts[lab]))
        print('  %-44s FINAL -> %s' % (name, sts['pooled'] if sts['copy 0'] == sts['copy 1'] else 'CANNOT_TELL (copies disagree)'))
    wired = np.array([len(x) > 0 for x in nb39])
    for N in (2, 4):
        r1, r3 = {}, {}
        for lab, cps in CP.items():
            cD, _ = concs(cps, 'D', N, 39, False); cMa, _ = concs(cps, 'M', N, 39, False); cMc, kept = concs(cps, 'M', N, 39, True)
            keepP = wired & (cD > 0)
            r1[lab] = bm((cD - cMc)[keepP], rng); r3[lab] = bm((cD - cMa)[keepP], rng)
            if lab == 'pooled': print('  N=%d: %d pairs enter (wired, idle concurrence > 0); R4 shots kept (median) %.1f%%' % (N, keepP.sum(), 100 * np.nanmedian(kept[keepP])))
        report('R1 N=%d clean-shot concurrence loss (K=39)' % N, r1); report('R3 N=%d all-shot concurrence loss (K=39)' % N, r3)
    r2 = {}
    for lab, cps in CP.items():
        cD, _ = concs(cps, 'D', 4, 39, False); c39, _ = concs(cps, 'M', 4, 39, False); c1, _ = concs(cps, 'M', 4, 1, False); keepP = cD > 0
        r2[lab] = bm(((cD - c39) - (cD - c1))[keepP], rng)
    report('R2 N=4 loss(K=39) - loss(K=1)', r2)
    out = []
    for c, (r, z) in copies.items():
        for i, m in enumerate(r['meta']):
            if m['arm'] == 'P' and m['mode'] == 'M' and m['K'] == 39:
                for k in range(m['N']):
                    mm = z['m%d_%d' % (k, i)]
                    for j, e in enumerate(E): out.append(((maha2(mm[:, j], mu0[e], cov0[e]) > D2) & (maha2(mm[:, j], mu1[e], cov1[e]) > D2)).mean())
    print('  R5 mid readouts outside both clusters: %.4f%%' % (100 * np.mean(out)))
