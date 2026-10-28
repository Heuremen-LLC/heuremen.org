"""Check H exactly as CHECK_H_DISTANCE_AND_COUNT.md."""
import json, os
from collections import deque
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = 2000; D2 = 3.5 ** 2; rng = np.random.default_rng(20261030)
def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi
def state(lo, hi): return 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
def boot(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 5: return np.nan, np.nan, np.nan, len(x), 'CANNOT TELL'
    b = np.median(rng.choice(x, (B, len(x))), axis=1); lo, hi = np.percentile(b, [2.5, 97.5]); return np.median(x), lo, hi, len(x), state(lo, hi)
def boot_diff(a, b_):
    a = a[np.isfinite(a)]; b_ = b_[np.isfinite(b_)]
    if min(len(a), len(b_)) < 5: return np.nan, np.nan, np.nan, 'CANNOT TELL'
    d = np.median(rng.choice(a, (B, len(a))), 1) - np.median(rng.choice(b_, (B, len(b_))), 1); lo, hi = np.percentile(d, [2.5, 97.5])
    return np.median(a) - np.median(b_), lo, hi, state(lo, hi)
def graph(props):
    adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    return adj
def dist(adj, E, S):
    d = {e: 0 for e in E}; q = deque(E)
    while q:
        u = q.popleft()
        for v in adj.get(u, ()):
            if v not in d: d[v] = d[u] + 1; q.append(v)
    return np.array([d.get(s, 99) for s in S])
def bits_packed(arr, k): nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def maha2(x, mu, cov):
    d = np.stack([x.real - mu.real, x.imag - mu.imag], -1); return np.einsum('...i,ij,...j->...', d, np.linalg.inv(cov), d)
def cov2(x): return np.cov(np.stack([x.real, x.imag])) + 1e-12 * np.eye(2)
def run1(chip):
    rec = json.load(open(os.path.join(HERE, 'job_%s.json' % chip))); z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); meta = rec['meta']
    E, S = rec['E'], rec['S']; nS = len(S); adj = graph(json.load(open(os.path.join(HERE, 'backend_properties_%s.json' % chip))))
    def ex(mode, N, st):
        ii = [i for i, m in enumerate(meta) if m['arm'] == 'S' and m['mode'] == mode and m['N'] == N and m['setting'] == st]
        return np.mean(np.concatenate([np.stack([1 - 2 * bits_packed(z['fin_%d' % i], k) for k in range(nS)], 1) for i in ii]), 0)
    X = {(mo, N): (ex(mo, N, 'XX'), ex(mo, N, 'YY')) for mo in 'DM' for N in (2, 4, 8)}
    clean = None
    return E, S, adj, X, clean
def run2(chip):
    rec = json.load(open(os.path.join(HERE, 'iqjob_%s.json' % chip))); z = np.load(os.path.join(HERE, 'iqraw_%s.npz' % chip)); meta = rec['meta']
    E, S = rec['E'], rec['S']; nS = len(S); n = len(E) + nS; adj = graph(json.load(open(os.path.join(HERE, 'iq_backend_properties_%s.json' % chip))))
    idx = lambda **kw: [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    C0 = np.concatenate([z['c_%d' % i] for i in idx(arm='CAL', state='0')]); C1 = np.concatenate([z['c_%d' % i] for i in idx(arm='CAL', state='1')])
    mu0 = C0.mean(0); u = C1.mean(0) - mu0; proj = lambda x, q: np.real((x - mu0[q]) * np.conj(u[q])) / np.abs(u[q]) ** 2
    cov0 = [cov2(C0[:, q][proj(C0[:, q], q) < 0.5]) for q in range(n)]
    Epos = {e: i for i, e in enumerate(E)}; nb = [[Epos[e] for e in adj.get(s, ()) if e in Epos] for s in S]
    def ex(mode, N, st, clean):
        num = np.zeros(nS); den = np.zeros(nS)
        for i in idx(arm='S', mode=mode, N=N, setting=st):
            fin = z['fin_%d' % i]; v = 1 - 2 * np.stack([(proj(fin[:, k], S[k]) > 0.5).astype(int) for k in range(nS)], 1)
            if clean:
                ok = np.ones((fin.shape[0], len(E)), bool)
                for k in range(N):
                    m = z['m%d_%d' % (k, i)]
                    for e in range(len(E)): ok[:, e] &= maha2(m[:, e], mu0[E[e]], cov0[E[e]]) < D2
                keep = np.stack([ok[:, x].all(1) if x else np.zeros(len(v), bool) for x in nb], 1)
            else: keep = np.ones_like(v, bool)
            num += (v * keep).sum(0); den += keep.sum(0)
        with np.errstate(invalid='ignore', divide='ignore'): return np.where(den > 0, num / np.maximum(den, 1), np.nan)
    X = {(mo, N): (ex(mo, N, 'XX', False), ex(mo, N, 'YY', False)) for mo in 'DM' for N in (2, 4, 8)}
    Xc = {N: (ex('M', N, 'XX', True), ex('M', N, 'YY', True)) for N in (2, 8)}
    return E, S, adj, X, Xc
out = {}
for chip in ('ibm_fez', 'ibm_kingston', 'ibm_marrakesh'):
    print('=' * 100); print(chip)
    for name, fn in (('run 1', run1), ('run 2', run2)):
        E, S, adj, X, Xc = fn(chip); d = dist(adj, E, S)
        ph = lambda mo, N: np.arctan2(X[(mo, N)][1], X[(mo, N)][0]); co = lambda mo, N: np.hypot(*X[(mo, N)])
        tot = {N: np.degrees(wrap(ph('M', N) - ph('D', N))) for N in (2, 4, 8)}; kick = {N: tot[N] / N for N in tot}
        r = {}
        print('  %s  spectators at distance 1/2/3+: %d/%d/%d' % (name, (d == 1).sum(), (d == 2).sum(), (d >= 3).sum()))
        for lab, sel in (('d=1', d == 1), ('d=2', d == 2), ('d>=3', d >= 3)):
            m, lo, hi, k, s = boot(kick[8][sel]); print('    H1 kick N=8 %-5s %+.3f deg [%+.3f, %+.3f] n=%d -> %s' % (lab, m, lo, hi, k, s))
        m, lo, hi, s = boot_diff(kick[8][d == 1], kick[8][d >= 2]); r['H1'] = s
        print('    H1 kick(d=1) - kick(d>=2)      %+.3f deg [%+.3f, %+.3f] -> %s' % (m, lo, hi, s))
        m, lo, hi, k, s = boot(kick[8] - kick[2]); r['H2'] = s; print('    H2 kick(N=8) - kick(N=2)       %+.3f deg [%+.3f, %+.3f] n=%d -> %s' % (m, lo, hi, k, s))
        Ns = np.array([2, 4, 8.]); Y = np.stack([tot[N] for N in (2, 4, 8)], 1); A = np.stack([Ns, np.ones(3)], 1)
        coef = np.linalg.lstsq(A, Y.T, rcond=None)[0]
        m, lo, hi, k, s = boot(coef[0]); print('    H2 slope (deg per readout)      %+.3f [%+.3f, %+.3f] n=%d -> %s' % (m, lo, hi, k, s))
        m, lo, hi, k, s = boot(coef[1]); r['H2i'] = s; print('    H2 intercept (deg)              %+.3f [%+.3f, %+.3f] n=%d -> %s' % (m, lo, hi, k, s))
        if Xc is None:
            loss = {N: co('D', N) - co('M', N) for N in (2, 8)}; lab = 'all spectators, all shots'
        else:
            loss = {N: co('D', N) - np.hypot(*Xc[N]) for N in (2, 8)}; lab = 'wired spectators, clean shots'
        m, lo, hi, k, s = boot(loss[8] - loss[2]); r['H3'] = s; print('    H3 loss(8) - loss(2) [%s] %+.4f [%+.4f, %+.4f] n=%d -> %s' % (lab, m, lo, hi, k, s))
        out[(chip, name)] = r
    print('  REPORTED (run 1 and run 2 agree):', {key: (out[(chip, 'run 1')][key] if out[(chip, 'run 1')][key] == out[(chip, 'run 2')][key] else 'CANNOT TELL (runs differ: %s / %s)' % (out[(chip, 'run 1')][key], out[(chip, 'run 2')][key])) for key in ('H1', 'H2', 'H2i', 'H3')})
