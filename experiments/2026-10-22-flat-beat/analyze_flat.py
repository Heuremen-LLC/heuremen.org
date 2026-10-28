"""Analysis exactly as PLAN_FLAT.md. python analyze_flat.py [folder]. Reads job_flat_ibm_phoenix.json, raw_flat_ibm_phoenix.npz,
props_flat_ibm_phoenix.json. Unit = one pair (A with B1, A with B2): each A gives two pairs. Per pair, the beat frequency f is fitted
freely: C_P(tau) = C_Z(tau) |cos(pi f tau)| over f = 1.00-100.00 kHz in steps of 0.05 kHz, on taus where C_Z is above 5x its noise floor.
Readings carry 95 % bootstrap intervals (2000 resamples over pairs) and four outcomes; copies compared as in the earlier runs."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'zz_revival'))
import analyze_zz as Z   # its bits/expect/coh helpers and the fixed zz_table (pair = its reverse; no leading zeros)
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)); B = 2000; RNG = np.random.default_rng(20261034)
F = np.arange(1.0, 100.0001, 0.05) * 1e3
def boot(x, f):
    x = np.asarray(x, float); ok = np.all(np.isfinite(x), axis=-1) if x.ndim > 1 else np.isfinite(x); x = x[ok]
    if len(x) < 5: return np.nan, np.nan, np.nan, len(x)
    idx = RNG.integers(0, len(x), (B, len(x))); b = np.array([f(x[i]) for i in idx]); lo, hi = np.percentile(b, [2.5, 97.5]); return float(f(x)), lo, hi, len(x)
def oc(lo, hi, n): return 'CANNOT TELL' if n < 5 or not np.isfinite(lo) else 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
class J2(Z.Job):
    def __init__(s, chip):
        s.rec = json.load(open(os.path.join(HERE, 'job_flat_%s.json' % chip))); s.z = np.load(os.path.join(HERE, 'raw_flat_%s.npz' % chip))
        s.props = json.load(open(os.path.join(HERE, 'props_flat_%s.json' % chip))); s.trip = [tuple(t) for t in s.rec['info']['triples']]; s.m = len(s.trip)
        s.meta = s.rec['meta']; s.taus = sorted({m['tau_us'] for m in s.meta}); s.edges = {(min(a, c), max(a, c)) for (a, c) in s.rec['edges']}
def pairs(J, copies):
    T = np.array(J.taus) * 1e-6; floor = 5 * np.sqrt(2.0 / (np.pi * 2 * J.rec['shots']))
    C = {arm: np.array([Z.coh(J, arm, t, copies)[0] for t in J.taus]) for arm in ('Z', 'P1', 'P2')}
    zz = Z.zz_table(J.props, J.edges); out = []
    for j, (b1, a, b2) in enumerate(J.trip):
        for b, arm in ((b1, 'P1'), (b2, 'P2')):
            cal = zz.get((min(a, b), max(a, b)), np.nan); live = C['Z'][:, j] > floor
            if live.sum() < 8 or not np.isfinite(cal) or cal <= 0: out.append((np.nan, cal)); continue
            model = C['Z'][live, j][None, :] * np.abs(np.cos(np.pi * F[:, None] * T[live][None, :]))
            out.append((F[np.argmin(((model - C[arm][live, j][None, :]) ** 2).sum(1))], cal))
    return np.array(out)   # (pairs, 2): fitted f, calibrated zz (Hz)
med = lambda x: np.median(x[:, 0]) if x.ndim > 1 else np.median(x)
def slope(x):
    lf, lz = np.log(x[:, 0]), np.log(x[:, 1]); A = np.vstack([lz - lz.mean(), np.ones_like(lz)]).T; return np.linalg.lstsq(A, lf, rcond=None)[0][0]
def spear(x):
    r = lambda v: np.argsort(np.argsort(v)); a, b = r(x[:, 0]), r(x[:, 1]); return np.corrcoef(a, b)[0, 1]
if __name__ == '__main__':
    J = J2('ibm_phoenix'); print('=' * 100); print('flat ibm_phoenix job', J.rec['job_id'], '| triples', J.m, '| pairs', 2 * J.m)
    P = pairs(J, (0, 1)); P0 = pairs(J, (0,)); P1 = pairs(J, (1,))
    ok = np.isfinite(P[:, 0]); print('  pairs fitted: %d of %d; calibrated zz (kHz): median %.1f, range %.1f-%.1f' % (ok.sum(), len(P), np.nanmedian(P[:, 1]) / 1e3, np.nanmin(P[:, 1]) / 1e3, np.nanmax(P[:, 1]) / 1e3))
    rows = [('Q1 median fitted beat (kHz)', lambda x: med(x) / 1e3, None),
            ('Q1 median fitted beat minus 12.5 kHz (kHz)', lambda x: med(x) / 1e3 - 12.5, 'd'),
            ('Q2 log-log slope of fitted beat on calibrated zz', slope, None),
            ('Q2 slope minus 1 (beat tracks the calibration)', lambda x: slope(x) - 1, 'd'),
            ('Q3 rank correlation, fitted beat vs calibrated zz', spear, None),
            ('Q4 median of fitted beat / calibrated zz', lambda x: np.median(x[:, 0] / x[:, 1]), None)]
    for name, f, kind in rows:
        v, lo, hi, n = boot(P, f); g = boot(P0[:, 0] - P1[:, 0], lambda d: np.median(d)) if kind is None or True else None
        d0 = f(P0[np.all(np.isfinite(P0), 1)]); d1 = f(P1[np.all(np.isfinite(P1), 1)])
        print('  %-52s %+.4g [%+.4g, %+.4g] n=%d -> %s   (copy 0 %+.4g, copy 1 %+.4g; per-pair fitted-beat copy gap %+.3g kHz [%+.3g, %+.3g])'
              % (name, v, lo, hi, n, oc(lo, hi, n) if kind or name.startswith('Q2 log') or name.startswith('Q3') else 'level', d0, d1, g[0] / 1e3, g[1] / 1e3, g[2] / 1e3))
    print('  per pair (A, B, calibrated zz kHz, fitted beat kHz):')
    k = 0
    for (b1, a, b2) in J.trip:
        for b in (b1, b2):
            print('    %3d %3d  %6.2f  %6.2f' % (a, b, P[k, 1] / 1e3, P[k, 0] / 1e3)); k += 1
