"""Check I exactly as CHECK_I_RUN_GAP.md. Reuses check H's readers unchanged."""
import os, sys
import numpy as np
sys.argv = sys.argv[:1]; HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'check_h_distance_and_count.py')).read().split('out = {}')[0]; exec(src)
rng = np.random.default_rng(20261031)
def per_run(fn, chip):
    E, S, adj, X, Xc = fn(chip); d = dist(adj, E, S)
    ph = lambda mo, N: np.arctan2(X[(mo, N)][1], X[(mo, N)][0]); co = lambda mo, N: np.hypot(*X[(mo, N)])
    tot = {N: np.degrees(wrap(ph('M', N) - ph('D', N))) for N in (2, 4, 8)}
    A = np.stack([np.array([2, 4, 8.]), np.ones(3)], 1); coef = np.linalg.lstsq(A, np.stack([tot[N] for N in (2, 4, 8)], 1).T, rcond=None)[0]
    v = dict(k8=tot[8] / 8, k2=tot[2] / 2, slope=coef[0], icpt=coef[1], dl=(co('D', 8) - co('M', 8)) - (co('D', 2) - co('M', 2)))
    if Xc is not None: v['dl_clean'] = (co('D', 8) - np.hypot(*Xc[8])) - (co('D', 2) - np.hypot(*Xc[2]))
    return d, v
def gap_then_pool(a, b, label, unit):
    g = boot(a - b); p = boot((a + b) / 2)
    if g[4] == 'CONSISTENT':
        print('    %-34s run gap %+.4f [%+.4f, %+.4f] -> none resolved | POOLED %+.4f [%+.4f, %+.4f] %s n=%d -> %s' % (label, g[0], g[1], g[2], p[0], p[1], p[2], unit, p[3], p[4])); return p
    print('    %-34s run gap %+.4f [%+.4f, %+.4f] -> runs differ | run1 %+.4f run2 %+.4f -> CANNOT TELL' % (label, g[0], g[1], g[2], np.median(a), np.median(b))); return None
for chip in ('ibm_fez', 'ibm_kingston', 'ibm_marrakesh'):
    print('=' * 100); print(chip)
    d, r1 = per_run(run1, chip); _, r2 = per_run(run2, chip)
    neg = all(np.median(r[k][sel]) < 0 for r in (r1, r2) for k in ('k8',) for sel in (d == 1, d == 2))
    print('    kick below 0 in every group in both runs: %s' % neg)
    p = gap_then_pool(r1['k8'], r2['k8'], 'kick per readout, N=8', 'deg')
    pk = (r1['k8'] + r2['k8']) / 2
    if neg:
        m, lo, hi, s = boot_diff(-pk[d == 1], -pk[d == 2])
        print('    H1 size(d=1) - size(d=2), pooled     %+.3f deg [%+.3f, %+.3f] (n %d / %d) -> %s' % (m, lo, hi, (d == 1).sum(), (d == 2).sum(),
              {'ABOVE': 'ABOVE: larger next to a measured qubit', 'BELOW': 'BELOW: larger two steps away', 'CONSISTENT': 'CONSISTENT: neither resolved'}.get(s, s)))
        for lab, sel in (('d=1', d == 1), ('d=2', d == 2)): m, lo, hi, k, s = boot(-pk[sel]); print('       size %s pooled %.3f deg [%.3f, %.3f] n=%d' % (lab, m, lo, hi, k))
    else: print('    H1 CANNOT TELL (sign not the same in every group)')
    q = gap_then_pool(r1['k8'] - r1['k2'], r2['k8'] - r2['k2'], 'H2 kick(8) - kick(2)', 'deg')
    if q and p: print('       bound: |change| <= %.3f deg = %.0f %% of the kick size %.3f deg' % (max(abs(q[1]), abs(q[2])), 100 * max(abs(q[1]), abs(q[2])) / abs(p[0]), abs(p[0])))
    gap_then_pool(r1['slope'], r2['slope'], 'H2 slope', 'deg/readout')
    gap_then_pool(r1['icpt'], r2['icpt'], 'H2 intercept', 'deg')
    gap_then_pool(r1['dl'], r2['dl'], 'H3 loss(8)-loss(2) all shots', '')
    m, lo, hi, k, s = boot(r2['dl_clean']); print('    H3 run 2 clean shots, wired (beside, not mixed)  %+.4f [%+.4f, %+.4f] n=%d -> %s' % (m, lo, hi, k, s))
