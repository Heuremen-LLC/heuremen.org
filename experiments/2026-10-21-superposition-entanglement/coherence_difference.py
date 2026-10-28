"""Plain coherence difference (analysis after the data, at an outside referee's request): for each unpaired untouched qubit,
coherence |<X>+i<Y>| after N readouts minus after N equal idles; median over qubits, 95 % bootstrap. Reads job_*.json and raw_*.npz."""
import json, os, glob, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(5)
def bits(arr, k):
    nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def med(x): b = [np.median(rng.choice(x, len(x))) for _ in range(2000)]; return np.median(x), np.percentile(b, [2.5, 97.5])
for jf in sorted(glob.glob(os.path.join(HERE, 'job_*.json'))):
    rec = json.load(open(jf)); chip = rec['chip']; z = np.load(os.path.join(HERE, 'raw_%s.npz' % chip))
    nP = len(rec['pairs']); SG = rec['singles']; meta = rec['meta']; sc = {}
    for i, m in enumerate(meta):
        if m['setting'][0] not in 'XY': continue
        arr = z['fin_%d' % i]
        for j in range(len(SG)): sc.setdefault((m['mode'], m['N'], j, m['setting'][0]), []).append((1 - 2 * bits(arr, 2 * nP + j)).mean())
    c = lambda mo, N: np.array([abs(complex(np.mean(sc[(mo, N, j, 'X')]), np.mean(sc[(mo, N, j, 'Y')]))) for j in range(len(SG))])
    out = []
    for N in (8, 16, 32):
        m, ci = med(c('M0', N) - c('D', N)); out.append('N=%2d %+.4f [%+.4f, %+.4f]' % (N, m, *ci))
    print('%-14s readout minus idle coherence (median, 95%%): ' % chip + ' | '.join(out))
    print('%-14s idle-arm coherence at N=32: median %.3f' % (chip, np.median(c('D', 32))))
