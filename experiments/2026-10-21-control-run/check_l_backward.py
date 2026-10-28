"""Check L exactly as CHECK_L_BACKWARD.md. Reuses check H's two readers unchanged."""
import os, sys
import numpy as np
sys.argv = sys.argv[:1]; HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'check_h_distance_and_count.py')).read().split('out = {}')[0]; exec(src)
rng = np.random.default_rng(20261034)
def b(x):
    x = x[np.isfinite(x)]; bb = np.median(rng.choice(x, (2000, len(x))), 1); lo, hi = np.percentile(bb, [2.5, 97.5])
    return np.median(x), lo, hi, len(x), ('ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT')
for chip in ('ibm_fez', 'ibm_kingston', 'ibm_marrakesh'):
    print('=' * 100); print(chip)
    for name, fn in (('control run', run1), ('raw-readout run', run2)):
        E, S, adj, X, Xc = fn(chip); Ns = np.array([2., 4., 8.])
        ph = lambda mo: np.unwrap(np.stack([np.arctan2(X[(mo, N)][1], X[(mo, N)][0]) for N in (2, 4, 8)]), axis=0)
        pD, pM = ph('D'), ph('M')
        sD = np.degrees(np.polyfit(Ns, pD, 1)[0]); kick = np.degrees(np.polyfit(Ns, np.unwrap(pM - pD, axis=0), 1)[0])
        fw, bw = sD >= 1, sD <= -1
        f = b(kick[fw]); g = b(kick[bw]); sg = np.sign(f[0]) or 1.0
        d = np.median(rng.choice(kick[bw], (2000, bw.sum())), 1) - np.median(rng.choice(kick[fw], (2000, fw.sum())), 1); lo, hi = np.percentile(d, [2.5, 97.5])
        print('  %-16s forward (n=%d) %+.3f [%+.3f, %+.3f] -> %s | backward (n=%d) %+.3f [%+.3f, %+.3f] -> %s | backward - forward %+.3f [%+.3f, %+.3f] | backward same sign as forward: %s'
              % (name, f[3], f[0], f[1], f[2], f[4], g[3], g[0], g[1], g[2], g[4], g[0] - f[0], lo, hi,
                 'RESOLVED' if (sg * g[1] > 0 if sg > 0 else sg * g[2] > 0) else 'NOT RESOLVED'))
