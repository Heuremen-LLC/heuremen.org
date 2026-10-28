"""Check K exactly as CHECK_K_SEGMENTS.md."""
import glob, os, numpy as np
import analyze as A
def segs(J):
    S = J.info['S']; nS = len(S); Ns = sorted({m['N'] for m in J.meta if m['N'] > 0})
    def ph(c, mode, N): return np.arctan2(J.mean_pm(c, 'fin', nS, mode=mode, N=N, setting='YY'), J.mean_pm(c, 'fin', nS, mode=mode, N=N, setting='XX'))
    def total(c):
        d = {N: A.wrap(ph(c, 'M', N) - ph(c, 'D', N)) for N in Ns}; t = {Ns[0]: d[Ns[0]]}
        for a, b in zip(Ns[:-1], Ns[1:]): t[b] = t[a] + A.wrap(d[b] - d[a])
        return {N: np.degrees(v) for N, v in t.items()}
    A.reading('RK.1 kick per readout, readouts 1 -> 9 (offset-free)', lambda c: (total(c)[9] - total(c)[1]) / 8, 'deg')
    A.reading('RK.2 kick per readout, readouts 27 -> 81 (offset-free)', lambda c: (total(c)[81] - total(c)[27]) / 54, 'deg')
for jf in sorted(glob.glob(os.path.join(A.HERE, 'job_count_*.json'))):
    chip = os.path.basename(jf)[len('job_count_'):-5]; J = A.Job('count', chip); J.props_n = J.rec.get('n_qubits', 156)
    print('=' * 100); print('count', chip); A.CUR[0] = chip; segs(J)
A.summary()
