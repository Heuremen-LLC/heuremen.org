"""Analysis exactly as PREREG.md (control run, 21 Oct 2026). Written before any result file was opened; its hash is committed first.
Reads job_<chip>.json and raw_<chip>.npz in this folder (or the folder given as the first argument). Prints every reading on copy 0,
copy 1 and pooled, each with one of four outcomes. Standard numpy only."""
import json, glob, os, sys
import numpy as np

HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
B = 2000

def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)

def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi

def boot_median(x, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 3: return np.nan, (np.nan, np.nan), len(x)
    b = np.median(rng.choice(x, (B, len(x))), axis=1)
    return float(np.median(x)), tuple(np.percentile(b, [2.5, 97.5])), len(x)

def state(ci, ref):
    lo, hi = ci
    if not np.isfinite(lo): return 'CANNOT_TELL'
    return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'

def conc(c):
    p = c['ZZ']; XX, YY, XY, YX = c['XX'], c['YY'], c['XY'], c['YX']
    r03 = abs(complex(XX - YY, -(XY + YX))) / 4; r12 = abs(complex(XX + YY, XY - YX)) / 4
    return 2 * max(0.0, r03 - np.sqrt(p[1] * p[2]), r12 - np.sqrt(p[0] * p[3]))

def analyse(rec, z, out):
    chip = rec['chip']; meta = rec['meta']; S = rec['S']; P = rec['pairs']; SG = rec['singles']; E = rec['E']
    nS, nP = len(S), len(P)
    rng = np.random.default_rng(20261021)
    def idx(**kw):
        return [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    def expv(arm, mode, N, setting, copies, K=None):
        """mean of (1 - 2 bit) per spectator over the given copies (shots pooled)"""
        kw = dict(arm=arm, mode=mode, N=N, setting=setting)
        if K is not None: kw['K'] = K
        arrs = [z['fin_%d' % i] for c in copies for i in idx(copy=c, **kw)]
        assert arrs, (arm, mode, N, setting, copies, K)
        a = np.concatenate(arrs)
        return np.array([(1 - 2 * bits(a, k)).mean() for k in range(nS)])
    def phase(arm, mode, N, copies, K=None):
        x = expv(arm, mode, N, 'XX', copies, K); y = expv(arm, mode, N, 'YY', copies, K); return np.arctan2(y, x), np.hypot(x, y)
    CP = {'copy 0': [0], 'copy 1': [1], 'pooled': [0, 1]}
    print('=' * 100, file=out); print('%s  job %s  E %d  S %d  pairs %d  L=%s dt' % (chip, rec['job_id'], len(E), nS, nP, rec['L_dt']), file=out)
    results = {}
    def report(name, per_copy, ref, unit=''):
        st = {}
        for lab in ('copy 0', 'copy 1', 'pooled'):
            m, ci, n = per_copy[lab]; st[lab] = state(ci, ref)
            print('  %-44s %-7s %+.4f%s [%+.4f, %+.4f] n=%d -> %s vs %g' % (name, lab, m, unit, ci[0], ci[1], n, st[lab], ref), file=out)
        final = st['pooled'] if st['copy 0'] == st['copy 1'] else 'CANNOT_TELL'
        if st['copy 0'] != st['copy 1']: print('  %-44s copies disagree -> CANNOT_TELL' % name, file=out)
        results[name] = (per_copy['pooled'][0], per_copy['pooled'][1], final); return final
    # R1 kick and R3 stretch (arm S)
    print('R1 KICK (deg per readout, arm S, N = 8)', file=out)
    kick, rate = {}, {}
    for lab, cps in CP.items():
        pM, _ = phase('S', 'M', 8, cps); pD8, _ = phase('S', 'D', 8, cps); pD0, _ = phase('S', 'D', 0, cps)
        kick[lab] = np.degrees(wrap(pM - pD8)) / 8; rate[lab] = np.degrees(wrap(pD8 - pD0)) / 8
    report('R1 kick', {lab: boot_median(kick[lab], rng) for lab in CP}, 0, ' deg')
    # R2 dose (arm C vs S)
    print('R2 DOSE (arm C, N = 8): kick with K measured vs the 39-measured kick', file=out)
    k39 = boot_median(kick['pooled'], rng)[0]
    for K in (1, 10):
        kk = {}
        for lab, cps in CP.items():
            pM, _ = phase('C', 'M', 8, cps, K); pD, _ = phase('C', 'D', 8, cps, K); kk[lab] = np.degrees(wrap(pM - pD)) / 8
        bm = {lab: boot_median(kk[lab], rng) for lab in CP}
        report('R2 kick with K=%d measured' % K, bm, 0, ' deg')
        lo, hi = bm['pooled'][1]; glob_ok = lo <= k39 <= hi; prop = k39 * K / 39; prop_ok = lo <= prop <= hi
        verdict = 'BOTH (cannot tell)' if glob_ok and prop_ok else 'GLOBAL' if glob_ok else 'PER-QUBIT' if prop_ok else 'NEITHER'
        print('  R2 K=%d: interval [%+.3f, %+.3f]; GLOBAL predicts %+.3f, PER-QUBIT predicts %+.3f -> %s' % (K, lo, hi, k39, prop, verdict), file=out)
        results['R2 dose K=%d' % K] = (bm['pooled'][0], (lo, hi), verdict)
    # R3 stretch: cross-copy regression
    print('R3 STRETCH: slope of kick on rate, rate and kick from different copies', file=out)
    for xa, ya in (('copy 0', 'copy 1'), ('copy 1', 'copy 0')):
        x, y = rate[xa], kick[ya]; a = np.polyfit(x, y, 1)[0]
        bs = [np.polyfit(x[i], y[i], 1)[0] for i in [rng.integers(0, len(x), len(x)) for _ in range(B)]]; ci = tuple(np.percentile(bs, [2.5, 97.5]))
        print('  rate from %s, kick from %s: slope a = %+.4f [%+.4f, %+.4f] -> %s vs 0; %s vs 0.0909; distance to 0 %.1f widths, to 0.0909 %.1f widths'
              % (xa, ya, a, *ci, state(ci, 0), state(ci, 0.0909), abs(a) / ((ci[1] - ci[0]) / 3.92), abs(a - 0.0909) / ((ci[1] - ci[0]) / 3.92)), file=out)
    back = rate['pooled'] < 0
    report('R3 kick on backward-precessing spectators', {lab: boot_median(kick[lab][back], rng) for lab in CP}, 0, ' deg')
    # R4 superposition cost (arm S)
    print('R4 SUPERPOSITION COST: coherence after readouts minus after equal idle (arm S)', file=out)
    diffs = {lab: {} for lab in CP}
    for lab, cps in CP.items():
        for N in (2, 4, 8):
            _, cM = phase('S', 'M', N, cps); _, cD = phase('S', 'D', N, cps); diffs[lab][N] = cM - cD
    for N in (2, 4, 8):
        report('R4 coherence M-D at N=%d' % N, {lab: boot_median(diffs[lab][N], rng) for lab in CP}, 0)
    slope = {lab: np.array([np.polyfit([2, 4, 8], [diffs[lab][N][j] for N in (2, 4, 8)], 1)[0] for j in range(nS)]) for lab in CP}
    report('R4 growth: slope of M-D against N', {lab: boot_median(slope[lab], rng) for lab in CP}, 0)
    # R5 entanglement cost (arm P)
    print('R5 ENTANGLEMENT COST: concurrence after readouts minus after equal idle (arm P)', file=out)
    def concs(mode, N, cps):
        cells = [dict() for _ in range(nP)]
        for st in ('ZZ', 'XX', 'YY', 'XY', 'YX'):
            a = np.concatenate([z['fin_%d' % i] for c in cps for i in idx(arm='P', mode=mode, N=N, setting=st, copy=c)])
            for j in range(nP):
                b1 = bits(a, 2 * j); b2 = bits(a, 2 * j + 1)
                if st == 'ZZ': cells[j]['ZZ'] = [np.mean((b1 == x) & (b2 == y)) for x, y in ((0, 0), (0, 1), (1, 0), (1, 1))]
                else: cells[j][st] = np.mean((1 - 2 * b1) * (1 - 2 * b2))
        return np.array([conc(c) for c in cells])
    for N in (2, 4):
        per = {}
        for lab, cps in CP.items():
            cD = concs('D', N, cps); cM = concs('M', N, cps); keep = cD > 0
            per[lab] = boot_median((cM - cD)[keep], rng)
            if lab == 'pooled': print('  N=%d: %d of %d pairs excluded (idle concurrence 0); idle median %.3f' % (N, (~keep).sum(), nP, np.median(cD)), file=out)
        report('R5 concurrence M-D at N=%d' % N, per, 0)
    # R6 leakage (arm L)
    print('R6 LEAKAGE: P(final=1 | last readout=1), leak arm minus reference, per measured qubit', file=out)
    def cond(N, e_init, cps):
        last = 'm%d' % (N - 1); ii = [i for c in cps for i in idx(arm='L', N=N, e_init=e_init, copy=c)]
        mid = np.concatenate([z['%s_%d' % (last, i)] for i in ii]); fin = np.concatenate([z['efin_%d' % i] for i in ii])
        r = []
        for k in range(len(E)):
            m1 = bits(mid, k) == 1
            r.append(bits(fin, k)[m1].mean() if m1.sum() >= 5 else np.nan)
        return np.array(r)
    report('R6 leakage', {lab: boot_median(cond(8, '0', cps) - cond(1, '1', cps), rng) for lab, cps in CP.items()}, 0)
    # R7 excitation rate (bound only)
    ii = idx(arm='S', mode='M', N=8)
    p = np.mean([bits(z['m%d_%d' % (k, i)], q).mean() for i in ii for k in range(8) for q in range(len(E))])
    print('R7 EXCITATION: fraction of 1 outcomes in the mid readouts (arm S, N = 8): %.4f' % p, file=out)
    results['R7 excitation p'] = (p, None, 'printed')
    return results

if __name__ == '__main__':
    out = open(os.path.join(HERE, 'results_prereg.txt'), 'w')
    class Tee:
        def write(self, s): out.write(s); sys.stdout.write(s)
        def flush(self): pass
    tee = Tee(); allres = {}
    for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
        rec = json.load(open(jf)); raw = os.path.join(HERE, 'raw_%s.npz' % rec['chip'])
        if not os.path.exists(raw): print(rec['chip'], 'NOTHING_RETURNED', file=tee); continue
        allres[rec['chip']] = analyse(rec, np.load(raw), tee)
    print('\nSUMMARY (pooled outcome; CANNOT_TELL where copies disagree)', file=tee)
    for chip, res in allres.items():
        for name, (m, ci, st) in res.items(): print('  %-14s %-46s %s' % (chip, name, st), file=tee)
