"""Check F exactly as CHECK_F_DOSE_COHERENCE.md."""
import json, glob, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = 2000
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def bm(x, rng):
    b = np.median(rng.choice(x, (B, len(x))), axis=1); return float(np.median(x)), tuple(np.percentile(b, [2.5, 97.5]))
st = lambda ci: 'ABOVE' if ci[0] > 0 else 'BELOW' if ci[1] < 0 else 'CONSISTENT'
for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); meta = rec['meta']
    nS = len(rec['S']); rng = np.random.default_rng(20261026)
    idx = lambda **kw: [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    def coh(arm, mode, K):
        kw = dict(arm=arm, mode=mode, N=8) if arm == 'S' else dict(arm=arm, mode=mode, N=8, K=K)
        ex = {s: np.mean(np.concatenate([np.stack([1 - 2 * bits(z['fin_%d' % i], k) for k in range(nS)], 1) for i in idx(setting=s, **kw)]), 0) for s in ('XX', 'YY')}
        return np.hypot(ex['XX'], ex['YY'])
    loss = {39: coh('S', 'D', 39) - coh('S', 'M', 39), 10: coh('C', 'D', 10) - coh('C', 'M', 10), 1: coh('C', 'D', 1) - coh('C', 'M', 1)}
    print(chip)
    for K in (1, 10, 39):
        m, ci = bm(loss[K], rng); print('   loss with %2d read out: %+.4f [%+.4f, %+.4f] %s' % (K, m, *ci, st(ci)))
    m, ci = bm(loss[39] - loss[1], rng)
    print('   loss(39) - loss(1): %+.4f [%+.4f, %+.4f] %s -> %s' % (m, *ci, st(ci), 'extra time per step EXCLUDED' if ci[0] > 0 else 'cannot tell'))
