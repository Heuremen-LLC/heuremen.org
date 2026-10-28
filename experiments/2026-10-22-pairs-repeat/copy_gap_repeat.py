"""Copy gap (copy 0 minus copy 1, with interval) for every pair-run-repeat reading the committed analysis labelled by its outcome labels. Defined after the data. Raw IQ: arxiv_one_readout/pairs_repeat (on request)."""
"""Check G exactly as CHECK_G_COPY_GAP.md. Same classification and quantities as analyze_iq.py / analyze_pairs_iq.py."""
import json, glob, os
import numpy as np
import sys
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/bones-forge/arxiv_one_readout/pairs_repeat'); B = 2000; D2 = 3.5 ** 2; rng = np.random.default_rng(20261029)
def maha2(x, mu, cov):
    d = np.stack([x.real - mu.real, x.imag - mu.imag], -1); return np.einsum('...i,ij,...j->...', d, np.linalg.inv(cov), d)
def cov2(x): return np.cov(np.stack([x.real, x.imag])) + 1e-12 * np.eye(2)
def conc(c):
    p = c['ZZ']; XX, YY, XY, YX = c['XX'], c['YY'], c['XY'], c['YX']
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))
def gap(v0, v1, label):
    d = np.asarray(v0) - np.asarray(v1); d = d[np.isfinite(d)]
    b = np.median(rng.choice(d, (B, len(d))), axis=1); lo, hi = np.percentile(b, [2.5, 97.5])
    out = 'no shift resolved -> pooled reading stands' if lo <= 0 <= hi else 'shift between copies -> CANNOT TELL stands'
    print('  %-48s copy0-copy1 median %+.4f [%+.4f, %+.4f] n=%d -> %s' % (label, np.median(d), lo, hi, len(d), out))
def calib(cal0, cal1, n):
    mu0 = cal0.mean(0); mu1 = cal1.mean(0); u = mu1 - mu0
    proj = lambda x, q: np.real((x - mu0[q]) * np.conj(u[q])) / np.abs(u[q]) ** 2
    cov0 = [cov2(cal0[:, q][proj(cal0[:, q], q) < 0.5]) for q in range(n)]
    return mu0, proj, cov0
def adjacency(props):
    adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    return adj
print('PAIR RUN (concurrence)')
for chip, cells in (('ibm_fez', ('R1-2','R1-4', 'R2')), ('ibm_kingston', ('R1-2', 'R1-4','R2')), ('ibm_marrakesh', ('R1-2','R1-4','R2'))):
    copies = {c: (json.load(open(os.path.join(HERE, 'pairsjob_%s_copy%d.json' % (chip, c)))), np.load(os.path.join(HERE, 'pairsraw_%s_copy%d.npz' % (chip, c)))) for c in (0, 1)}
    r0 = copies[0][0]; E = r0['E']; P = r0['pairs']; nP = len(P); n = len(E) + len(r0['S']); K1 = r0['K1']
    adj = adjacency(json.load(open(os.path.join(HERE, 'pairs_backend_properties_%s.json' % chip)))); Epos = {e: i for i, e in enumerate(E)}
    nb39 = [[Epos[e] for e in (adj.get(a, set()) | adj.get(b, set())) if e in Epos] for a, b in P]; wired = np.array([len(x) > 0 for x in nb39])
    cal = {st: np.concatenate([z['c_%d' % i] for (r, z) in copies.values() for i, m in enumerate(r['meta']) if m['arm'] == 'CAL' and m['state'] == st]) for st in '01'}
    mu0, proj, cov0 = calib(cal['0'], cal['1'], n)
    def concs(c, mode, N, K, clean):
        r, z = copies[c]; cells_ = [dict() for _ in range(nP)]
        for st in ('ZZ', 'XX', 'YY', 'XY', 'YX'):
            acc = [[] for _ in range(nP)]
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
                    acc[j].append(((proj(fin[:, 2 * j], a) > 0.5).astype(int)[keep], (proj(fin[:, 2 * j + 1], b) > 0.5).astype(int)[keep]))
            for j in range(nP):
                b1 = np.concatenate([x for x, _ in acc[j]]); b2 = np.concatenate([y for _, y in acc[j]])
                if st == 'ZZ': cells_[j]['ZZ'] = [np.mean((b1 == x) & (b2 == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))]
                else: cells_[j][st] = np.mean((1 - 2 * b1) * (1 - 2 * b2))
        return np.array([conc(x) for x in cells_])
    for cell in cells:
        if cell.startswith('R1'):
            N = int(cell[-1]); v = []
            for c in (0, 1):
                cD = concs(c, 'D', N, 39, False); v.append(np.where(wired & (cD > 0), cD - concs(c, 'M', N, 39, True), np.nan))
            gap(v[0], v[1], '%s N=%d clean-shot concurrence loss' % (chip, N))
        else:
            v = []
            for c in (0, 1):
                cD = concs(c, 'D', 4, 39, False); v.append(np.where(cD > 0, concs(c, 'M', 4, 1, False) - concs(c, 'M', 4, 39, False), np.nan))
            gap(v[0], v[1], '%s N=4 loss(39)-loss(1)' % chip)
