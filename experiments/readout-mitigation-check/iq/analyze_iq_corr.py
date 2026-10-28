"""Analysis exactly as PREREG_IQ.md. Reads iqjob_<chip>.json, iqraw_<chip>.npz and iq_backend_properties_<chip>.json
from the folder given as the first argument (default: this folder). Standard numpy only."""
import json, glob, os, sys
import numpy as np
HERE = '/Users/heuremen/bones-forge/arxiv_one_readout/control_run/'; MULT = float(sys.argv[1]); B = 2000; CAL_STORE = {}; D2 = 3.5 ** 2

def bm(x, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 3: return np.nan, (np.nan, np.nan), len(x)
    b = np.median(rng.choice(x, (B, len(x))), axis=1); return float(np.median(x)), tuple(np.percentile(b, [2.5, 97.5])), len(x)
def state(ci):
    lo, hi = ci
    return 'CANNOT_TELL' if not np.isfinite(lo) else 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
def maha2(x, mu, cov):
    d = np.stack([x.real - mu.real, x.imag - mu.imag], -1); ic = np.linalg.inv(cov)
    return np.einsum('...i,ij,...j->...', d, ic, d)
def cov2(x):
    return np.cov(np.stack([x.real, x.imag])) + 1e-12 * np.eye(2)

for jf in sorted(glob.glob(os.path.join(HERE, 'iqjob_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; meta = rec['meta']; raw = os.path.join(HERE, 'iqraw_%s.npz' % chip)
    if not os.path.exists(raw): print(chip, 'NOTHING_RETURNED'); continue
    z = np.load(raw); E = rec['E']; S = rec['S']; nS = len(S); n = len(E) + nS
    props = json.load(open(os.path.join(HERE, 'iq_backend_properties_%s.json' % chip)))
    adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    Epos = {e: i for i, e in enumerate(E)}; nb = [[Epos[e] for e in adj.get(s, ()) if e in Epos] for s in S]; wired = np.array([len(x) > 0 for x in nb])
    idx = lambda **kw: [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    rng = np.random.default_rng(20261027)
    # calibration (both copies)
    C0 = np.concatenate([z['c_%d' % i] for i in idx(arm='CAL', state='0')]); C1 = np.concatenate([z['c_%d' % i] for i in idx(arm='CAL', state='1')])
    mu0 = C0.mean(0); mu1 = C1.mean(0); u = mu1 - mu0
    proj = lambda x, q: np.real((x - mu0[q]) * np.conj(u[q])) / np.abs(u[q]) ** 2
    cov0 = [cov2(C0[:, q][proj(C0[:, q], q) < 0.5]) for q in range(n)]; cov1 = [cov2(C1[:, q][proj(C1[:, q], q) > 0.5]) for q in range(n)]
    print('=' * 96); print('%s  job %s' % (chip, rec['job_id']))
    def spec_bits(i):
        fin = z['fin_%d' % i]; return np.stack([(proj(fin[:, k], S[k]) > 0.5).astype(int) for k in range(nS)], 1)
    def mid_ok(i, N):   # shots x 39: True if that measured qubit was verified 0 on every one of the N readouts
        ok = np.ones((z['fin_%d' % i].shape[0], len(E)), bool)
        for k in range(N):
            m = z['m%d_%d' % (k, i)]
            for e in range(len(E)): ok[:, e] &= maha2(m[:, e], mu0[E[e]], cov0[E[e]]) < D2
        return ok
    def coh(mode, N, copies, clean):
        ex = {}
        for st in ('XX', 'YY'):
            num = np.zeros(nS); den = np.zeros(nS)
            for i in [i for c in copies for i in idx(arm='S', mode=mode, N=N, setting=st, copy=c)]:
                v = 1 - 2 * spec_bits(i)
                if clean and N > 0:
                    ok = mid_ok(i, N); keep = np.stack([ok[:, x].all(1) if x else np.ones(len(v), bool) for x in nb], 1)
                else: keep = np.ones_like(v, bool)
                num += (v * keep).sum(0); den += keep.sum(0)
            v = num / np.maximum(den, 1)
            e01 = MULT * np.array([(proj(C0[:, S[k]], S[k]) > 0.5).mean() for k in range(nS)]); e10 = MULT * np.array([(proj(C1[:, S[k]], S[k]) < 0.5).mean() for k in range(nS)])
            CAL_STORE[chip] = (e01 / max(MULT, 1e-9), e10 / max(MULT, 1e-9))
            p1 = (1 - v) / 2; ex[st] = 1 - 2 * (p1 - e01) / (1 - e01 - e10)
        return np.hypot(ex['XX'], ex['YY']), den
    CP = {'copy 0': [0], 'copy 1': [1], 'pooled': [0, 1]}
    for N in (2, 4, 8):
        res = {}
        for lab, cps in CP.items():
            cD, tot = coh('D', N, cps, False); cMa, _ = coh('M', N, cps, False); cMc, kept = coh('M', N, cps, True)
            res[lab] = (bm((cD - cMc)[wired], rng), bm((cD - cMa)[wired], rng), np.median((kept / np.maximum(tot, 1))[wired]))
        for name, j in (('R1 clean-shot loss', 0), ('R2 all-shot loss', 1)):
            sts = {lab: state(res[lab][j][1]) for lab in CP}
            final = sts['pooled'] if sts['copy 0'] == sts['copy 1'] else 'CANNOT_TELL (copies disagree)'
            for lab in CP:
                m, ci, k = res[lab][j]
                print('  N=%d %-20s %-7s %+.4f [%+.4f, %+.4f] n=%d -> %s' % (N, name, lab, m, *ci, k, sts[lab]))
            print('  N=%d %-20s FINAL -> %s' % (N, name, final))
        print('  N=%d R4 shots kept by the selection (median over spectators, pooled): %.1f%%' % (N, 100 * res['pooled'][2]))
        out = []; tot = 0
        for i in idx(arm='S', mode='M', N=N):
            for k in range(N):
                m = z['m%d_%d' % (k, i)]
                for e in range(len(E)):
                    q = E[e]; out.append(((maha2(m[:, e], mu0[q], cov0[q]) > D2) & (maha2(m[:, e], mu1[q], cov1[q]) > D2)).mean())
        print('  N=%d R3 mid readouts outside both clusters: %.4f%%' % (N, 100 * np.mean(out)))

for c,(a,b) in CAL_STORE.items(): print('CAL %s median e01 %.4f e10 %.4f (spectators, in-job)' % (c, np.median(a), np.median(b)))
