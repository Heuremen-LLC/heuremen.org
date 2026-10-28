"""dose_read.py -- R1 kick (arm S, 39 read) and R2 dose (arm C, K = 1 and 10 read) from a control-run job fired with --arms S,C.
The reading code below is copied unchanged from the committed control-run analyze.py (functions bits, wrap, boot_median, state, idx,
expv, phase, report and the R1/R2 blocks); only the arms this job holds are read. Usage: dose_read.py [FOLDER]  (job_*.json + raw_*.npz)"""
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
def analyse(rec, z, out):
    meta = rec['meta']; S = rec['S']; nS = len(S); rng = np.random.default_rng(20261021)
    def idx(**kw): return [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    def expv(arm, mode, N, setting, copies, K=None):
        kw = dict(arm=arm, mode=mode, N=N, setting=setting)
        if K is not None: kw['K'] = K
        arrs = [z['fin_%d' % i] for c in copies for i in idx(copy=c, **kw)]; assert arrs, (arm, mode, N, setting, copies, K)
        a = np.concatenate(arrs); return np.array([(1 - 2 * bits(a, k)).mean() for k in range(nS)])
    def phase(arm, mode, N, copies, K=None):
        x = expv(arm, mode, N, 'XX', copies, K); y = expv(arm, mode, N, 'YY', copies, K); return np.arctan2(y, x), np.hypot(x, y)
    CP = {'copy 0': [0], 'copy 1': [1], 'pooled': [0, 1]}; results = {}
    print('=' * 100, file=out); print('%s  job %s  S %d' % (rec['chip'], rec['job_id'], nS), file=out)
    def report(name, per_copy, ref, unit=''):
        st = {}
        for lab in ('copy 0', 'copy 1', 'pooled'):
            m, ci, n = per_copy[lab]; st[lab] = state(ci, ref)
            print('  %-44s %-7s %+.4f%s [%+.4f, %+.4f] n=%d -> %s vs %g' % (name, lab, m, unit, ci[0], ci[1], n, st[lab], ref), file=out)
        final = st['pooled'] if st['copy 0'] == st['copy 1'] else 'CANNOT_TELL'
        if st['copy 0'] != st['copy 1']: print('  %-44s copies disagree -> CANNOT_TELL' % name, file=out)
        results[name] = (per_copy['pooled'][0], per_copy['pooled'][1], final); return final
    kick = {}
    for lab, cps in CP.items():
        pM, _ = phase('S', 'M', 8, cps); pD8, _ = phase('S', 'D', 8, cps); kick[lab] = np.degrees(wrap(pM - pD8)) / 8
    report('R1 kick', {lab: boot_median(kick[lab], rng) for lab in CP}, 0, ' deg')
    k39 = boot_median(kick['pooled'], rng)[0]
    for K in (1, 10):
        kk = {}
        for lab, cps in CP.items():
            pM, _ = phase('C', 'M', 8, cps, K); pD, _ = phase('C', 'D', 8, cps, K); kk[lab] = np.degrees(wrap(pM - pD)) / 8
        bm = {lab: boot_median(kk[lab], rng) for lab in CP}
        final = report('R2 kick with K=%d measured' % K, bm, 0, ' deg')
        g = kk['copy 0'] - kk['copy 1']; gb = boot_median(g, rng)
        print('  R2 K=%d copy gap (copy 0 - copy 1, per qubit): %+.4f [%+.4f, %+.4f] -> %s' % (K, gb[0], gb[1][0], gb[1][1], 'copies shifted' if state(gb[1], 0) != 'CONSISTENT' else 'no shift'), file=out)
    return results
if __name__ == '__main__':
    for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
        rec = json.load(open(jf)); raw = os.path.join(HERE, 'raw_%s.npz' % rec['chip'])
        if os.path.exists(raw): analyse(rec, np.load(raw), sys.stdout)
        else: print(rec['chip'], 'NOTHING_RETURNED')
