"""Check E exactly as CHECK_E_WIRED_VS_UNWIRED.md."""
import json, glob, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = 2000
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip)); meta = rec['meta']
    props = json.load(open(os.path.join(HERE, 'backend_properties_%s.json' % chip)))
    T2 = {q: next((p['value'] for p in pl if p['name'] == 'T2'), np.nan) for q, pl in enumerate(props['qubits'])}
    adj = {}
    for g in props['gates']:
        if g['gate'] == 'cz': a, b = g['qubits']; adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)
    E = set(rec['E']); S = rec['S']; nS = len(S); rng = np.random.default_rng(20261025)
    wired = np.array([len(adj.get(s, set()) & E) > 0 for s in S]); t2 = np.array([T2[s] for s in S])
    idx = lambda **kw: [i for i, m in enumerate(meta) if all(m.get(k) == v for k, v in kw.items())]
    def coh(mode, N):
        ex = {st: np.mean(np.concatenate([np.stack([1 - 2 * bits(z['fin_%d' % i], k) for k in range(nS)], 1) for i in idx(arm='S', mode=mode, N=N, setting=st)]), 0) for st in ('XX', 'YY')}
        return np.hypot(ex['XX'], ex['YY'])
    print(chip)
    for N in (2, 4, 8):
        cD, cM = coh('D', N), coh('M', N)
        ok = (cD > 0.05) & (cM > 0.05) & np.isfinite(t2)
        tau = np.where(ok, -t2 / N * np.log(np.where(ok, cM / cD, 1)), np.nan)
        w = tau[wired & ok]; u = tau[~wired & ok]
        d = np.median(w) - np.median(u)
        bs = [np.median(rng.choice(w, len(w))) - np.median(rng.choice(u, len(u))) for _ in range(B)]; lo, hi = np.percentile(bs, [2.5, 97.5])
        st = 'ABOVE -> extra time / chip-wide cause EXCLUDED as the explanation of the wired excess' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT -> cannot tell'
        print('   N=%d  implied extra time per readout: wired %.3f us (n=%d), unwired %.3f us (n=%d); difference %+.3f [%+.3f, %+.3f] us -> %s (dropped %d)'
              % (N, np.median(w), len(w), np.median(u), len(u), d, lo, hi, st, (~ok).sum()))
