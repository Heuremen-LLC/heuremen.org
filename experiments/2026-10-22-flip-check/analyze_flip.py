"""Analysis exactly as PLAN_FLIP.md (X-FLIP CHECK). Reads flipjob_<chip>.json, flipraw_<chip>.npz, flip_backend_properties_<chip>.json
from the folder given as the first argument (default: this folder). Standard numpy only. Ruler and classification as control_run/analyze_iq.py.
Heuremen LLC, Wayfinder and Bones, 22 October 2026."""
import json, glob, os, sys
import numpy as np
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)); B = 2000; D2 = 3.5 ** 2; NS = (4, 8)

def bm(x, rng, f=np.median):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 3: return np.nan, (np.nan, np.nan), len(x)
    b = f(rng.choice(x, (B, len(x))), axis=1); return float(f(x)), tuple(np.percentile(b, [2.5, 97.5])), len(x)
def state(ci):
    lo, hi = ci
    return 'CANNOT_TELL' if not np.isfinite(lo) else 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
def maha2(x, mu, cov):
    d = np.stack([x.real - mu.real, x.imag - mu.imag], -1); ic = np.linalg.inv(cov)
    return np.einsum('...i,ij,...j->...', d, ic, d)
def cov2(x):
    return np.cov(np.stack([x.real, x.imag])) + 1e-12 * np.eye(2)
# KNOWN-GOOD reference: 21 Oct raw-readout run (control_run/results_iq.txt), R1 clean-shot loss, pooled, where its FINAL was resolved (ABOVE)
KG_R1 = {('ibm_kingston', 4): (0.0056, 0.0013, 0.0105), ('ibm_kingston', 8): (0.0062, 0.0013, 0.0155),
         ('ibm_marrakesh', 4): (0.0087, 0.0032, 0.0160), ('ibm_marrakesh', 8): (0.0205, 0.0112, 0.0284)}
# Addendum before firing (PLAN_FLIP.md): marrakesh side reading without the untouched qubits wired to these measured qubits (fixed now)
SIDE_DROP = {'ibm_marrakesh': (82, 100, 119)}
SIDE_LABEL = 'marrakesh without neighbours of overlap/miss qubits (printed beside, does not replace)'
CP = {'copy 0': [0], 'copy 1': [1], 'pooled': [0, 1]}
def final(sts): return sts['pooled'] if sts['copy 0'] == sts['copy 1'] else 'CANNOT_TELL (copies disagree)'
VERDICT = {}

for jf in sorted(glob.glob(os.path.join(HERE, 'flipjob_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; meta = rec['meta']; raw = os.path.join(HERE, 'flipraw_%s.npz' % chip)
    if not os.path.exists(raw): print(chip, 'NOTHING_RETURNED'); continue
    z = np.load(raw); E = rec['E']; S = rec['S']; nS = len(S); n = len(E) + nS
    props = json.load(open(os.path.join(HERE, 'flip_backend_properties_%s.json' % chip)))
    adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    Epos = {e: i for i, e in enumerate(E)}; nb = [[Epos[e] for e in adj.get(s, ()) if e in Epos] for s in S]; wired = np.array([len(x) > 0 for x in nb])
    idx = lambda **kw: [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    rng = np.random.default_rng(20261028)
    # --- calibration (both copies), as analyze_iq.py
    C0 = np.concatenate([z['c_%d' % i] for i in idx(arm='CAL', state='0')]); C1 = np.concatenate([z['c_%d' % i] for i in idx(arm='CAL', state='1')])
    mu0 = C0.mean(0); mu1 = C1.mean(0); u = mu1 - mu0
    proj = lambda x, q: np.real((x - mu0[q]) * np.conj(u[q])) / np.abs(u[q]) ** 2
    cov0 = [cov2(C0[:, q][proj(C0[:, q], q) < 0.5]) for q in range(n)]; cov1 = [cov2(C1[:, q][proj(C1[:, q], q) > 0.5]) for q in range(n)]
    # --- NEW: post-flip cluster from CAL 0X: second read of shots whose first read was verified 0
    A = np.concatenate([z['c0_%d' % i] for i in idx(arm='CAL', state='0X')]); Bx = np.concatenate([z['c1_%d' % i] for i in idx(arm='CAL', state='0X')])
    mu1x = np.zeros(n, complex); cov1x = []
    for q in range(n):
        v = Bx[:, q][maha2(A[:, q], mu0[q], cov0[q]) < D2]; mu1x[q] = v.mean(); cov1x.append(cov2(v[proj(v, q) > 0.5]))
    # flip-readable guard (fixed now): a measured qubit whose CAL 0X second read (verified-0 shots) fails the flip rule in > 20 % of shots
    # cannot report a flip; it is left out of F1/F2, and every spectator wired to it is left out of L1/L2/L3
    cal_fail = np.zeros(len(E))
    for e, q in enumerate(E):
        v = Bx[:, q][maha2(A[:, q], mu0[q], cov0[q]) < D2]
        cal_fail[e] = np.mean((maha2(v, mu0[q], cov0[q]) < D2) | (maha2(v, mu1x[q], cov1x[q]) >= D2))
    good = cal_fail <= 0.20; useS = wired & np.array([all(good[x] for x in xs) for xs in nb])
    print('=' * 100); print('%s  job %s' % (chip, rec['job_id']))
    print('  CAL 0X flip-rule failure per measured qubit: median %.2f%%, max %.2f%%; flip-unreadable (> 20%%): %s; wired spectators %d, used for L1/L2/L3 %d'
          % (100 * np.median(cal_fail), 100 * cal_fail.max(), [E[e] for e in np.where(~good)[0]], wired.sum(), useS.sum()))
    _c = {}
    def cached(f):
        def g(*a):
            if (f.__name__,) + a not in _c: _c[(f.__name__,) + a] = f(*a)
            return _c[(f.__name__,) + a]
        return g
    @cached
    def spec_bits(i):
        fin = z['fin_%d' % i]; return np.stack([(proj(fin[:, k], S[k]) > 0.5).astype(int) for k in range(nS)], 1)
    @cached
    def mid_ok(i, N):   # shots x 39: True if that measured qubit was verified 0 on every one of its N mid readouts
        ok = np.ones((z['fin_%d' % i].shape[0], len(E)), bool)
        for k in range(N):
            m = z['m%d_%d' % (k, i)]
            for e in range(len(E)): ok[:, e] &= maha2(m[:, e], mu0[E[e]], cov0[E[e]]) < D2
        return ok
    @cached
    def flip(i, ref='0X'):   # shots x 39: (flipped, failed inside 0, failed outside both)
        ef = z['efin_%d' % i]; F = np.zeros(ef.shape, bool); I0 = np.zeros(ef.shape, bool); OB = np.zeros(ef.shape, bool)
        for e in range(len(E)):
            q = E[e]; m1, c1 = (mu1x[q], cov1x[q]) if ref == '0X' else (mu1[q], cov1[q])
            d0 = maha2(ef[:, e], mu0[q], cov0[q]) < D2; d1 = maha2(ef[:, e], m1, c1) < D2
            F[:, e] = d1 & ~d0; I0[:, e] = d0; OB[:, e] = ~d0 & ~d1
        return F, I0, OB
    def coh(mode, N, copies, use_mid, use_flip):
        ex = {}; kept = np.zeros(nS); tot = np.zeros(nS)
        for st in ('XX', 'YY'):
            num = np.zeros(nS); den = np.zeros(nS)
            for i in [i for c in copies for i in idx(arm='S', mode=mode, N=N, setting=st, copy=c)]:
                v = 1 - 2 * spec_bits(i); ok = np.ones((len(v), len(E)), bool)
                if use_mid and N > 0: ok &= mid_ok(i, N)
                if use_flip: ok &= flip(i)[0]
                keep = np.stack([ok[:, x].all(1) if x else np.ones(len(v), bool) for x in nb], 1)
                num += (v * keep).sum(0); den += keep.sum(0); kept += keep.sum(0); tot += len(v)
            ex[st] = num / np.maximum(den, 1)
        return np.hypot(ex['XX'], ex['YY']), kept / np.maximum(tot, 1)
    def frate(mode, N, copies, hidden, ref='0X'):   # per measured qubit: fraction of shots that failed to flip (and its two parts)
        f = np.zeros(len(E)); a = np.zeros(len(E)); o = np.zeros(len(E)); d = np.zeros(len(E))
        for i in [i for c in copies for i in idx(arm='S', mode=mode, N=N, copy=c)]:
            Fl, I0, OB = flip(i, ref); w = mid_ok(i, N) if (hidden and mode == 'M' and N > 0) else np.ones_like(Fl)
            f += (~Fl & w).sum(0); a += (I0 & w).sum(0); o += (OB & w).sum(0); d += w.sum(0)
        d = np.maximum(d, 1); return f / d, a / d, o / d
    # --- KNOWN-GOOD 1: the 0X calibration must flip (verified-0 shots land in the post-flip cluster) in nearly all shots: median >= 95 %
    print('  KNOWN-GOOD 0X: flip success, median over measured qubits %.2f%% (flip-readable %.2f%%) -> %s'
          % (100 * np.median(1 - cal_fail), 100 * np.median(1 - cal_fail[good]), 'PASS' if np.median(1 - cal_fail[good]) >= 0.95 else 'FAIL'))
    # --- F1: flip failure, per N, per arm (mean over the 39 measured qubits)
    for N in (0,) + NS:
        for mode in (('D',) if N == 0 else ('D', 'M')):
            for lab, cps in CP.items():
                f, a, o = frate(mode, N, cps, False); fc, _, _ = frate(mode, N, cps, False, ref='CAL1')
                f, a, o, fc = f[good], a[good], o[good], fc[good]
                line = '  F1 N=%d %-7s %-7s fail %.4f%% (inside 0 %.4f%%, outside both %.4f%%) | with CAL1 cluster %.4f%%' % (N, {'D': 'idle', 'M': 'readout'}[mode], lab, 100 * f.mean(), 100 * a.mean(), 100 * o.mean(), 100 * fc.mean())
                if mode == 'M':
                    fh = frate('M', N, cps, True)[0][good]; line += ' | among shots verified 0 on every mid read %.4f%%' % (100 * fh.mean())
                print(line)
    # --- KNOWN-GOOD 2: idle arm (no readouts) flip failure must sit at the calibration's own floor: per qubit F1(idle, N) - CAL 0X failure,
    # mean over flip-readable qubits, bootstrap over qubits; must not read ABOVE in the pooled FINAL
    for N in (0,) + NS:
        sts = {}; vals = {}
        for lab, cps in CP.items():
            vals[lab] = bm((frate('D', N, cps, False)[0] - cal_fail)[good], rng, np.mean); sts[lab] = state(vals[lab][1])
        m, ci, k = vals['pooled']; f = final(sts)
        print('  KNOWN-GOOD idle N=%d: flip failure minus CAL 0X floor %+.4f%% [%+.4f, %+.4f] -> %s -> %s'
              % (N, 100 * m, 100 * ci[0], 100 * ci[1], f, 'FAIL (idle above the calibration floor)' if f == 'ABOVE' else 'PASS'))
    for N in NS:
        res = {}; vec = {}
        for lab, cps in CP.items():
            fD, _, _ = frate('D', N, cps, False); fMh, _, _ = frate('M', N, cps, True); fM, _, _ = frate('M', N, cps, False)
            cD, _ = coh('D', N, cps, False, False); cDf, kDf = coh('D', N, cps, False, True)
            cM1, k1 = coh('M', N, cps, True, False); cM2, k2 = coh('M', N, cps, True, True)
            L1 = (cD - cM1)[useS]; L2 = (cDf - cM2)[useS]
            vec[lab] = (cD - cM1, cDf - cM2)
            res[lab] = dict(F2=bm((fMh - fD)[good], rng, np.mean), F2u=bm((fM - fD)[good], rng, np.mean), L1=bm(L1, rng), L2=bm(L2, rng), L3=bm(L2 - L1, rng),
                            L1w=bm((cD - cM1)[wired], rng), k1=np.median(k1[useS]), k2=np.median(k2[useS]), kD=np.median(kDf[useS]))
        fin = {}
        for name, desc in (('F2', 'hidden-leak rate (readout, verified-0 shots) - idle'), ('F2u', 'all-shot flip failure, readout - idle'),
                           ('L1w', 'loss, mid-verified-0, ALL wired (run_iq R1 repeat)'),
                           ('L1', 'loss, mid-verified-0 selection'), ('L2', 'loss, mid-verified-0 AND flipped'), ('L3', 'L2 - L1')):
            sts = {lab: state(res[lab][name][1]) for lab in CP}; fin[name] = final(sts)
            for lab in CP:
                m, ci, k = res[lab][name]; sc = 100 if name.startswith('F') else 1
                print('  N=%d %-4s %-48s %-7s %+.4f%s [%+.4f, %+.4f] n=%d -> %s' % (N, name, desc, lab, sc * m, '%' if sc == 100 else ' ', sc * ci[0], sc * ci[1], k, sts[lab]))
            print('  N=%d %-4s FINAL -> %s' % (N, name, fin[name]))
        if (chip, N) in KG_R1:
            r, lo, hi = KG_R1[(chip, N)]; m = res['pooled']['L1w'][0]
            print('  KNOWN-GOOD N=%d: L1w %+.4f vs 21 Oct clean-shot loss %+.4f [%+.4f, %+.4f] -> %s' % (N, m, r, lo, hi, 'PASS' if lo <= m <= hi else 'FAIL'))
        if chip in SIDE_DROP:   # side reading: its own rng, so every committed reading above is untouched
            srng = np.random.default_rng(20261029 + N); drop = set(SIDE_DROP[chip])
            notnb = np.array([not any(E[e] in drop for e in xs) for xs in nb]); uX = useS & notnb; wX = wired & notnb
            sres = {lab: dict(L1=bm(vec[lab][0][uX], srng), L2=bm(vec[lab][1][uX], srng), L3=bm((vec[lab][1] - vec[lab][0])[uX], srng),
                              L1w=bm(vec[lab][0][wX], srng)) for lab in CP}
            print('  N=%d --- %s: spectators used %d (L1/L2/L3), %d (L1w); dropped neighbours of %s' % (N, SIDE_LABEL, uX.sum(), wX.sum(), sorted(drop)))
            for name in ('L1w', 'L1', 'L2', 'L3'):
                sts = {lab: state(sres[lab][name][1]) for lab in CP}
                for lab in CP:
                    m, ci, k = sres[lab][name]
                    print('  N=%d SIDE %-4s %-7s %+.4f [%+.4f, %+.4f] n=%d -> %s' % (N, name, lab, m, ci[0], ci[1], k, sts[lab]))
                print('  N=%d SIDE %-4s FINAL -> %s' % (N, name, final(sts)))
            if (chip, N) in KG_R1:
                r, lo, hi = KG_R1[(chip, N)]; m = sres['pooled']['L1w'][0]
                print('  KNOWN-GOOD SIDE N=%d: L1w without those neighbours %+.4f vs 21 Oct clean-shot loss %+.4f [%+.4f, %+.4f] -> %s' % (N, m, r, lo, hi, 'PASS' if lo <= m <= hi else 'FAIL'))
        p = res['pooled']; print('  N=%d shots kept (median over wired spectators, pooled): L1 readout %.1f%%, L2 readout %.1f%%, L2 idle %.1f%%' % (N, 100 * p['k1'], 100 * p['k2'], 100 * p['kD']))
        if fin['L2'] == 'ABOVE': v = 'L2 ABOVE: a neighbour stuck in a level that reads as 0 is EXCLUDED as the cause (to the extent the X-flip catches it)'
        elif fin['L2'] == 'CONSISTENT' and fin['L1'] == 'ABOVE': v = 'L2 CONSISTENT while L1 ABOVE: the loss MAY come from neighbours stuck in a level that reads as 0'
        else: v = 'L1 %s, L2 %s: neither branch of the rule; reported as is' % (fin['L1'], fin['L2'])
        v += ' | L3 %s, F2 %s' % (fin['L3'], fin['F2'])
        VERDICT[(chip, N)] = (fin, v); print('  N=%d DECISION -> %s' % (N, v))
