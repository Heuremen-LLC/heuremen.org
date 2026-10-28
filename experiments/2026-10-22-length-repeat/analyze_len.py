"""analyze_len.py -- committed analysis for the two-length run (PLAN_LENGTH.md). Usage: analyze_len.py CHIP [CHIP ...] [--dir D]
Kick per readout k_L = wrap(phase_M - phase_D)/N in degrees per untouched qubit (phase = atan2(<Y>, <X>)), median over qubits, 95% bootstrap
over qubits (2000, seed 20261023). Copy gap = copy 0 minus copy 1; an interval excluding 0 makes the reading CANNOT_TELL."""
import sys, json, os, numpy as np
B = 2000; rng = np.random.default_rng(20261023)
D = sys.argv[sys.argv.index('--dir') + 1] if '--dir' in sys.argv else os.path.dirname(os.path.abspath(__file__))
def bits(arr, k): nby = arr.shape[1]; return ((arr[:, nby - 1 - k // 8] >> (k % 8)) & 1).astype(int)
def wrap(x): return (x + np.pi) % (2 * np.pi) - np.pi
def lab(lo, hi, ref=0.0): return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'
def load(chip):
    rec = json.load(open(f'{D}/job_{chip}.json')); z = np.load(f'{D}/raw_{chip}.npz'); return rec, z
def ev(rec, z, reg, nbits, copies, **kw):
    ix = [i for i, m in enumerate(rec['meta']) if m['copy'] in copies and all(m[k] == v for k, v in kw.items())]
    a = np.concatenate([z[f'{reg}_{i}'] for i in ix]); return np.array([(1 - 2 * bits(a, k)).mean() for k in range(nbits)])
def kicks(rec, z, N, length, copies):
    nS = len(rec['S']); ph = {}
    for mode in ('D', 'M'):
        x = ev(rec, z, 'fin', nS, copies, arm='S', N=N, length=length, mode=mode, basis='X'); y = ev(rec, z, 'fin', nS, copies, arm='S', N=N, length=length, mode=mode, basis='Y')
        ph[mode] = (np.arctan2(y, x), np.hypot(x, y))
    k = np.degrees(wrap(ph['M'][0] - ph['D'][0])) / N; alive = (ph['D'][1] > 5 / np.sqrt(rec['shots'])) & (ph['M'][1] > 5 / np.sqrt(rec['shots']))
    return np.where(alive, k, np.nan)
def boot(f, n):
    idx = rng.integers(0, n, (B, n)); return np.array([f(i) for i in idx])
def analyse(chip, out):
    rec, z = load(chip); dt = rec['dt']; L1 = rec['L_dt']['L1'] * dt * 1e6; L2 = rec['L_dt']['L2'] * dt * 1e6; rho = L2 / L1
    print(f'\n=== {chip}  job {rec["job_id"]}  L1 {L1:.3f} us  L2 {L2:.3f} us  L2/L1 {rho:.3f}  shots {rec["shots"]}', file=out); R = {'chip': chip, 'rho': rho, 'L1': L1, 'L2': L2}
    for N in (4, 8):
        k = {(L, c): kicks(rec, z, N, L, c) for L in ('L1', 'L2') for c in ((0,), (1,), (0, 1))}
        ok = np.isfinite(k[('L1', (0, 1))]) & np.isfinite(k[('L2', (0, 1))]); n = int(ok.sum())
        a1, a2 = k[('L1', (0, 1))][ok], k[('L2', (0, 1))][ok]
        for L, a in (('L1', a1), ('L2', a2)):
            g = (k[(L, (0,))] - k[(L, (1,))])[ok]; gb = np.percentile(boot(lambda i: np.nanmedian(g[i]), n), [2.5, 97.5])
            mb = boot(lambda i: np.median(a[i]), n); lo, hi = np.percentile(mb, [2.5, 97.5]); ct = gb[0] > 0 or gb[1] < 0
            print(f'  R1 kick N={N} {L}: {np.median(a):+.3f} [{lo:+.3f}, {hi:+.3f}] deg/readout  n={n}  copy gap [{gb[0]:+.3f}, {gb[1]:+.3f}] -> {"CANNOT_TELL (copies disagree)" if ct else lab(lo, hi)}', file=out)
            R[(N, L)] = (np.median(a), lo, hi, ct)
        rb = boot(lambda i: np.median(a2[i]) / np.median(a1[i]), n); r = np.median(a2) / np.median(a1); lo, hi = np.percentile(rb, [2.5, 97.5])
        c1 = lo <= 1 <= hi; cr = lo <= rho <= hi; ct = R[(N, 'L1')][3] or R[(N, 'L2')][3]
        st = ('CANNOT_TELL (copies disagree)' if ct else 'CANNOT_TELL (both 1 and L2/L1 inside)' if (c1 and cr) else 'PER-READOUT' if c1 else 'PROPORTIONAL-TO-LENGTH' if cr else
              'NEITHER: between' if min(1, rho) < r < max(1, rho) else 'NEITHER: outside both')
        print(f'  R2 ratio k(L2)/k(L1) N={N}: {r:.3f} [{lo:.3f}, {hi:.3f}]  (per readout 1, proportional {rho:.3f}) -> {st}', file=out)
        bb = boot(lambda i: (np.median(a1[i]) - np.median(a2[i])) / (L1 - L2), n); bpt = (np.median(a1) - np.median(a2)) / (L1 - L2)
        ab = np.array([np.median(a1[i]) for i in rng.integers(0, n, (B, n))]) - bb * L1; apt = np.median(a1) - bpt * L1
        print(f'  R3 split N={N}: per-readout part a {apt:+.3f} [{np.percentile(ab,2.5):+.3f}, {np.percentile(ab,97.5):+.3f}] deg, per-us part b {bpt:+.3f} [{np.percentile(bb,2.5):+.3f}, {np.percentile(bb,97.5):+.3f}] deg/us', file=out)
        R[(N, 'state')] = st
    E = rec['E']
    for N in (4, 8):
        for L in ('L1', 'L2'):
            pd = {c: (1 - ev(rec, z, 'efin', len(E), c, arm='B', N=N, length=L, mode='D', basis='X')) / 2 for c in ((0,), (1,), (0, 1))}
            pm = {c: (1 - ev(rec, z, 'efin', len(E), c, arm='B', N=N, length=L, mode='M', basis='X')) / 2 for c in ((0,), (1,), (0, 1))}
            d = pd[(0, 1)] - pm[(0, 1)]; g = (pd[(0,)] - pm[(0,)]) - (pd[(1,)] - pm[(1,)]); n = len(d)
            mb = boot(lambda i: np.median(d[i]), n); lo, hi = np.percentile(mb, [2.5, 97.5]); gb = np.percentile(boot(lambda i: np.median(g[i]), n), [2.5, 97.5])
            ct = gb[0] > 0 or gb[1] < 0
            print(f'  R4 measured qubit, P(1) idle minus readout, N={N} {L}: {np.median(d):+.4f} [{lo:+.4f}, {hi:+.4f}] n={n} (idle {np.median(pd[(0,1)]):.3f}, readout {np.median(pm[(0,1)]):.3f}) -> {"CANNOT_TELL (copies disagree)" if ct else lab(lo, hi)}', file=out)
    return R
if __name__ == '__main__':
    chips = [a for a in sys.argv[1:] if a.startswith('ibm_') or a.startswith('synth')]
    for c in chips: analyse(c, sys.stdout)
