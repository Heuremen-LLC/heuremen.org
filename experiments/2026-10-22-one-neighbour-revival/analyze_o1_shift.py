"""analyze_o1_shift.py -- DEFINED AFTER THE JOB, BEFORE THIS QUANTITY WAS READ (22 Oct 2026). The one-neighbour run (PLAN_ZZ.md) holds an arm
O1 with neighbour B1 in |1> beside arm Z (B1 in |0>), same A, same taus, same job; its committed analysis read only the coherence ratio
C(O1)/C(Z) (R6). This reads the PHASE: per triple, the frequency shift of A when B1 goes from 0 to 1 is the slope of
unwrap(phase_O1(tau) - phase_Z(tau)) against tau over taus where both arms' coherence exceeds 5/sqrt(shots); df = slope / (2 pi), in Hz.
That shift is the static zz of the A-B1 edge in the convention 'frequency change of A when B flips', with no beat fit and no superposed
neighbour. Readings per chip, unit = triple, 95% bootstrap over triples (2000), four outcomes:
  S1 median |df| (kHz) and median stated |zz_A-B1| (kHz)
  S2 median |df| / stated            (1 = IBM's zz in this convention; 0.5 or 2 = a factor-2 convention)
  S3 log-log slope of |df| on stated (1 = tracks edge by edge, whatever a constant convention factor)
  S4 rank correlation |df| vs stated
  copy gap: S2 and S3 recomputed per copy; their difference's interval excluding 0 -> CANNOT_TELL.
KNOWN-GOOD (fixed now): on ibm_fez the slope S3 must include 1 and S2 must lie in [0.4, 2.5]; otherwise this phase reading is not decided on
either chip."""
import sys, os, json, numpy as np
sys.argv = [sys.argv[0]] + sys.argv[1:]
import analyze_zz as Z
RNG = np.random.default_rng(20261041); B = 2000
def shifts(J, copies):
    T = np.array(J.taus) * 1e-6; ph = {}; co = {}
    for arm in ('Z', 'O1'):
        X = np.array([Z.expect(J, arm, t, 'X', copies)[0] for t in J.taus]); Y = np.array([Z.expect(J, arm, t, 'Y', copies)[0] for t in J.taus])
        ph[arm] = np.arctan2(Y, X); co[arm] = np.hypot(X, Y)
    shots = J.rec.get('shots', 1000) * len(copies); floor = 5 / np.sqrt(shots); out = np.full(J.m, np.nan)
    for j in range(J.m):
        live = (co['Z'][:, j] > floor) & (co['O1'][:, j] > floor)
        if live.sum() < 4: continue
        d = np.unwrap(ph['O1'][live, j] - ph['Z'][live, j]); t = T[live]
        out[j] = np.polyfit(t, d, 1)[0] / (2 * np.pi)
    return out
def lab(lo, hi, ref): return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'
def stats(df, st):
    ok = np.isfinite(df) & np.isfinite(st) & (st > 0) & (np.abs(df) > 0); a = np.abs(df[ok]); s = st[ok]; n = len(a)
    def f(i): x, y = a[i], s[i]; return np.median(x / y), np.polyfit(np.log(y), np.log(x), 1)[0], np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1]
    pt = f(np.arange(n)); bs = np.array([f(RNG.integers(0, n, n)) for _ in range(B)]); return n, a, s, pt, bs
if __name__ == '__main__':
    for chip in ('ibm_fez', 'ibm_phoenix'):
        J = Z.Job(chip); zz = Z.zz_table(J.props, J.edges)
        st = np.array([zz.get((min(a, b1), max(a, b1)), np.nan) for b1, a, b2 in J.trip])
        df = shifts(J, (0, 1)); d0 = shifts(J, (0,)); d1 = shifts(J, (1,))
        n, a, s, pt, bs = stats(df, st); lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
        _, _, _, p0, b0 = stats(d0, st); _, _, _, p1, b1 = stats(d1, st); g = np.percentile(b0 - b1, [2.5, 97.5], axis=0)
        ct = [(g[0, k] > 0 or g[1, k] < 0) for k in range(3)]
        print(f'=== {chip}  job {J.rec["job_id"]}  triples {J.m}, fitted {n}')
        print(f'  S1 median |df| {np.median(a)/1e3:.2f} kHz; median stated {np.median(s)/1e3:.2f} kHz')
        print(f'  S2 median |df|/stated {pt[0]:.3f} [{lo[0]:.3f}, {hi[0]:.3f}] vs 1 -> {"CANNOT_TELL (copies disagree)" if ct[0] else lab(lo[0], hi[0], 1)}')
        print(f'  S3 log-log slope {pt[1]:.3f} [{lo[1]:.3f}, {hi[1]:.3f}] vs 1 -> {"CANNOT_TELL (copies disagree)" if ct[1] else lab(lo[1], hi[1], 1)}; vs 0 -> {lab(lo[1], hi[1], 0)}')
        print(f'  S4 rank correlation {pt[2]:.3f} [{lo[2]:.3f}, {hi[2]:.3f}] vs 0 -> {lab(lo[2], hi[2], 0)}')
        print('  per triple (A, B1, stated kHz, |df| kHz):', ' '.join(f'({a_},{b_},{x/1e3:.1f},{abs(y)/1e3:.1f})' for (b_, a_, _), x, y in zip(J.trip, st, df)))
