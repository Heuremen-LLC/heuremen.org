"""Analysis exactly as PLAN_ZZ.md. python analyze_zz.py [folder]. Reads job_zz_<chip>.json, raw_zz_<chip>.npz, props_zz_<chip>.json.
Unit = one triple (B1, A, B2). Every reading: per triple, median over triples, 95 % bootstrap (2000), four outcomes; per-triple copy0 - copy1
interval contains 0 -> pooled value is the reading, excludes 0 -> CANNOT TELL with both copies printed. Signs carried."""
import json, glob, os, sys
import numpy as np
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)); B = 2000; RNG = np.random.default_rng(20261033)
S_GRID = np.arange(0.10, 2.5001, 0.0025)
def outcome(lo, hi, n):
    if n < 5 or not np.isfinite(lo): return 'CANNOT TELL'
    return 'ABOVE' if lo > 0 else 'BELOW' if hi < 0 else 'CONSISTENT'
def boot(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 5: return np.nan, np.nan, np.nan, len(x), 'CANNOT TELL'
    b = np.median(RNG.choice(x, (B, len(x))), 1); lo, hi = np.percentile(b, [2.5, 97.5]); return float(np.median(x)), lo, hi, len(x), outcome(lo, hi, len(x))
def bits(arr, n):
    return np.unpackbits(np.asarray(arr, np.uint8), axis=1, bitorder='big')[:, ::-1][:, :n]
def zz_table(props, edges):
    """zz_<ab> in GHz -> Hz, keyed by the coupling-map edge whose two numbers concatenate to <ab> (ambiguity -> assert)."""
    out = {}
    for g in props['general']:
        if not g['name'].startswith('zz_'): continue
        d = g['name'][3:]; hits = sorted({(min(int(d[:k]), int(d[k:])), max(int(d[:k]), int(d[k:]))) for k in range(1, len(d)) if all(x == '0' or x[0] != '0' for x in (d[:k], d[k:])) and (min(int(d[:k]), int(d[k:])), max(int(d[:k]), int(d[k:]))) in edges})   # FIX 22 Oct after data: a pair and its reverse are one edge; no qubit number has a leading zero
        assert len(hits) == 1, (g['name'], hits); a, c = hits[0]; out[(min(a, c), max(a, c))] = abs(g['value']) * 1e9
    return out
def qp(props, q, name):
    for p in props['qubits'][q]:
        if p['name'] == name: return p['value']
    return np.nan
class Job:
    def __init__(s, chip):
        s.rec = json.load(open(os.path.join(HERE, 'job_zz_%s.json' % chip))); s.z = np.load(os.path.join(HERE, 'raw_zz_%s.npz' % chip))
        s.props = json.load(open(os.path.join(HERE, 'props_zz_%s.json' % chip))); s.trip = [tuple(t) for t in s.rec['info']['triples']]; s.m = len(s.trip)
        s.meta = s.rec['meta']; s.taus = sorted({m['tau_us'] for m in s.meta})
        s.edges = {(min(a, c), max(a, c)) for (a, c) in s.rec.get('edges', [])}
    def idx(s, **kw): return [i for i, m in enumerate(s.meta) if all(m.get(k) == v for k, v in kw.items())]
    def reg(s, name, i): return bits(s.z['%s_%d' % (name, i)], s.m)
def expect(J, arm, tau, st, copies, cond=None, split=None):
    """<1-2a> per A qubit, optionally restricted to shots where register cond[0] bit == cond[1] (per qubit), or a seeded random half (split)."""
    num = np.zeros(J.m); den = np.zeros(J.m)
    for i in [i for c in copies for i in J.idx(arm=arm, tau_us=tau, setting=st, copy=c)]:
        a = J.reg('a', i); v = 1 - 2 * a.astype(float)
        if cond is not None: w = (J.reg(cond[0], i) == cond[1]).astype(float)
        elif split is not None: w = (np.random.default_rng(split + i).random(a.shape) < 0.5).astype(float)
        else: w = np.ones_like(v)
        num += (v * w).sum(0); den += w.sum(0)
    return num / np.maximum(den, 1), den
def coh(J, arm, tau, copies):
    x, n = expect(J, arm, tau, 'X', copies); y, _ = expect(J, arm, tau, 'Y', copies); return np.hypot(x, y), n
def coh_cond(J, arm, tau, copies, reg):
    """Coherence of A in each branch of the final bit of register reg, weighted by branch size."""
    tot = np.zeros(J.m); wsum = np.zeros(J.m)
    for b in (0, 1):
        x, n = expect(J, arm, tau, 'X', copies, cond=(reg, b)); y, _ = expect(J, arm, tau, 'Y', copies, cond=(reg, b))
        tot += np.hypot(x, y) * n; wsum += n
    return tot / np.maximum(wsum, 1)
def coh_half(J, arm, tau, copies):
    """Bias-matched none for coh_cond: the same statistic on two seeded random halves of the shots."""
    tot = np.zeros(J.m); wsum = np.zeros(J.m)
    for h in (0, 1):
        x0, n = expect(J, arm, tau, 'X', copies, split=11 + h * 1000); y0, _ = expect(J, arm, tau, 'Y', copies, split=11 + h * 1000)
        tot += np.hypot(x0, y0) * n; wsum += n
    return tot / np.maximum(wsum, 1)
def readings(J, copies):
    T = np.array(J.taus); taus_s = T * 1e-6
    C = {arm: np.array([coh(J, arm, t, copies)[0] for t in J.taus]) for arm in ('Z', 'P1', 'P2', 'P12', 'O1', 'EZ', 'EP1')}   # (ntau, m)
    zz = zz_table(J.props, J.edges); f1 = np.array([zz.get((min(a, b1), max(a, b1)), np.nan) for b1, a, b2 in J.trip])
    f2 = np.array([zz.get((min(a, b2), max(a, b2)), np.nan) for b1, a, b2 in J.trip])
    out = {}
    s_fit = np.full(J.m, np.nan); s_fit2 = np.full(J.m, np.nan)
    for j in range(J.m):
        for f, arm, store in ((f1[j], 'P1', s_fit), (f2[j], 'P2', s_fit2)):
            if not np.isfinite(f) or f <= 0: continue
            model = C['Z'][:, j][None, :] * np.abs(np.cos(np.pi * S_GRID[:, None] * f * taus_s[None, :]))
            store[j] = S_GRID[np.argmin(((model - C[arm][:, j][None, :]) ** 2).sum(1))]
    out['R1 fitted scale s (P1, B1 alone in +)'] = s_fit
    out['R1 fitted scale s (P2, B2 alone in +)'] = s_fit2
    out['R1 fitted scale s minus 1.0 (P1)'] = s_fit - 1.0
    out['R1 fitted scale s minus 0.5 (P1)'] = s_fit - 0.5
    T1b = np.array([qp(J.props, b1, 'T1') for b1, a, b2 in J.trip]); eps = np.array([qp(J.props, b1, 'readout_error') for b1, a, b2 in J.trip])
    h = np.full(J.m, np.nan); hp = np.full(J.m, np.nan); du = np.full(J.m, np.nan); dc = np.full(J.m, np.nan); ec = np.full(J.m, np.nan)
    pc = np.full(J.m, np.nan); o1 = np.full(J.m, np.nan)
    floor = 5 * np.sqrt(2.0 / (np.pi * 2 * J.rec['shots']))
    for j in range(J.m):
        if not np.isfinite(s_fit[j]): continue
        fr = s_fit[j] * f1[j]; ts = 1.0 / fr * 1e6; td = 0.5 / fr * 1e6
        if ts <= T[-2]:
            k = np.argsort(np.abs(T - ts))[:3]; ratio = C['P1'][k, j] / C['Z'][k, j]
            h[j] = np.mean(ratio); hp[j] = (1 + np.exp(-ts / T1b[j])) / 2
        if td <= T[-2]:
            k = int(np.argmin(np.abs(T - td))); t = J.taus[k]
            du[j] = C['P1'][k, j] / C['Z'][k, j]
            dc[j] = coh_cond(J, 'P1', t, copies, 'b1')[j] / coh_half(J, 'Z', t, copies)[j]
            ec[j] = C['EP1'][k, j] / C['EZ'][k, j] - 1
        alive = C['Z'][:, j] > floor
        if alive.sum() >= 5:
            x1 = C['P1'][alive, j] / C['Z'][alive, j]; x2 = C['P2'][alive, j] / C['Z'][alive, j]; x12 = C['P12'][alive, j] / C['Z'][alive, j]
            pc[j] = np.sum((x12 - x1 * x2) ** 2) - np.sum(((1 - x12) - ((1 - x1) + (1 - x2))) ** 2)
            o1[j] = np.mean(C['O1'][alive, j] / C['Z'][alive, j])
    out['R2 return height at revival / Z, minus ZZ+T1 prediction (1+exp(-t*/T1_B1))/2'] = h - hp
    out['R2 return height at revival / Z (as measured)'] = h
    out['R3a unconditional coherence at the dip / Z'] = du
    out['R3b coherence split by B1 final bit at the dip / same statistic on Z'] = dc
    out['R3c R3b - R3a (coherence held in the correlation with B1)'] = dc - du
    out['R4 echo: EP1/EZ - 1 at the dip (static ZZ predicts 0)'] = ec
    out['R5 SSE(product) - SSE(sum) for both in + (BELOW -> product)'] = pc
    out['R6 B1 in 1: mean C(O1)/C(Z) over live taus (reference)'] = o1
    return out, dict(f1=f1, f2=f2, n_rev=int(np.isfinite(h).sum()), n_dip=int(np.isfinite(du).sum()))
def run_chip(chip):
    J = Job(chip)
    print('=' * 100); print('zz', chip, 'job', J.rec['job_id'], '| triples', J.m)
    pooled, info = readings(J, (0, 1)); c0, _ = readings(J, (0,)); c1, _ = readings(J, (1,))
    print('  props ZZ (kHz) B1 edge: median %.2f, range %.2f-%.2f | revival inside grid for %d triples, dip for %d'
          % (np.nanmedian(info['f1']) / 1e3, np.nanmin(info['f1']) / 1e3, np.nanmax(info['f1']) / 1e3, info['n_rev'], info['n_dip']))
    res = {}
    for k in pooled:
        v, lo, hi, n, oc = boot(pooled[k]); g = boot(c0[k] - c1[k])
        if np.isfinite(g[1]) and (g[1] > 0 or g[2] < 0):
            print('  %-78s copies shifted %+.4g [%+.4g, %+.4g] | copy0 %+.4g copy1 %+.4g -> CANNOT TELL' % (k, g[0], g[1], g[2], np.nanmedian(c0[k]), np.nanmedian(c1[k]))); res[k] = None
        else:
            print('  %-78s %+.4g [%+.4g, %+.4g] n=%d -> %s   (copy gap %+.4g [%+.4g, %+.4g])' % (k, v, lo, hi, n, oc, g[0], g[1], g[2])); res[k] = (v, lo, hi)
    return res
if __name__ == '__main__':
    R = {}
    for jf in sorted(glob.glob(os.path.join(HERE, 'job_zz_*.json'))):
        chip = json.load(open(jf))['chip']
        if not os.path.exists(os.path.join(HERE, 'raw_zz_%s.npz' % chip)): print(chip, 'NOTHING RETURNED'); continue
        R[chip] = run_chip(chip)
    k = 'R1 fitted scale s (P1, B1 alone in +)'
    if len(R) == 2 and all(R[c].get(k) for c in R):
        (a, b) = list(R)
        va, la, ha = R[a][k]; vb, lb, hb = R[b][k]; sa = (ha - la) / 3.92; sb = (hb - lb) / 3.92; d = va - vb; sd = np.hypot(sa, sb)
        print('=' * 100); print('  R7 scale s, %s minus %s: %+.4f [%+.4f, %+.4f] -> %s' % (a, b, d, d - 1.96 * sd, d + 1.96 * sd, outcome(d - 1.96 * sd, d + 1.96 * sd, 99)))
