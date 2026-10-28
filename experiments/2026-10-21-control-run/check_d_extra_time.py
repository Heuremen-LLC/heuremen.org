"""Check D exactly as CHECK_D_EXTRA_TIME.md. a_hi per chip is read from results_prereg.txt (the committed analysis output).
Readout length L in microseconds = L_dt * 0.004 (matches IBM's stated 1.66 / 2.18 / 2.684 us)."""
import json, glob, os, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = 2000
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def bm(x, rng):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; b = np.median(rng.choice(x, (B, len(x))), axis=1)
    return float(np.median(x)), tuple(np.percentile(b, [2.5, 97.5])), len(x)
blocks = re.split(r'={20,}', open(os.path.join(HERE, 'results_prereg.txt')).read())
for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); meta = rec['meta']
    blk = [b for b in blocks if b.strip().startswith(chip)][0]
    a_hi = max(float(m) for m in re.findall(r'slope a = [+-][0-9.]+ \[[+-][0-9.]+, ([+-][0-9.]+)\]', blk))
    L = rec['L_dt'] * 0.004; dt_max = max(a_hi, 0.0) * L
    props = json.load(open(os.path.join(HERE, 'backend_properties_%s.json' % chip)))
    T2 = {q: next((p['value'] for p in plist if p['name'] == 'T2'), np.nan) for q, plist in enumerate(props['qubits'])}  # fix after first run crashed: a qubit with no T2 is excluded (bm drops nan)
    S = rec['S']; nS = len(S); rng = np.random.default_rng(20261024)
    idx = lambda **kw: [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    def coh(mode, N):
        ex = {st: np.mean(np.concatenate([np.stack([1 - 2 * bits(z['fin_%d' % i], k) for k in range(nS)], 1) for i in idx(arm='S', mode=mode, N=N, setting=st)]), 0) for st in ('XX', 'YY')}
        return np.hypot(ex['XX'], ex['YY'])
    t2 = np.array([T2[s] for s in S])
    print('%-14s a_hi %+.4f  L %.3f us  dt_max %.4f us per readout  median T2 %.0f us  (qubits without a T2: %d)' % (chip, a_hi, L, dt_max, np.nanmedian(t2), np.isnan(t2).sum()))
    for N in (2, 4, 8):
        cD = coh('D', N); loss = cD - coh('M', N); bound = cD * (1 - np.exp(-N * dt_max / t2))
        m, ci, n = bm(loss - bound, rng)
        st = 'EXCLUDED' if ci[0] > 0 else 'NOT EXCLUDED'
        print('   N=%d  median loss %+.4f  median extra-time bound %.4f  loss-bound %+.4f [%+.4f,%+.4f] n=%d -> %s' % (N, np.median(loss), np.median(bound), m, *ci, n, st))
