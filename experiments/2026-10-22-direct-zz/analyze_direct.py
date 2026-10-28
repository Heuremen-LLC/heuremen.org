"""analyze_direct.py -- committed analysis (PLAN_DIRECT.md). Usage: analyze_direct.py CHIP [--dir D] [--props P]
Per edge (A-B1 from arm O1, A-B2 from arm O2): phase of A = atan2(<Y>, <X>) per tau; O shots kept only where that B reads 1 at the end,
Z shots only where it reads 0 (heralded); df = slope of unwrap(phase_O - phase_Z) vs tau / 2 pi (Hz) over taus where both arms'
coherence > 5/sqrt(kept shots). Unit = edge; 95% bootstrap over TRIPLES (2000, the two edges of a triple resampled together).
Readings: D1 median |df| and stated |zz|; D2 log-log slope of |df| on stated (edges with stated > 1 kHz); D3 median |df|/stated;
D4 rank correlation; D5 edges stated 0.00: their |df|. Copies: D2 and D3 per copy; gap interval excluding 0 -> CANNOT_TELL."""
import sys, os, json, numpy as np
D = sys.argv[sys.argv.index('--dir') + 1] if '--dir' in sys.argv else os.path.dirname(os.path.abspath(__file__))
B = 2000; RNG = np.random.default_rng(20261026)
def bits(arr, n): return np.unpackbits(np.asarray(arr, np.uint8), axis=1, bitorder='big')[:, ::-1][:, :n]
def zz_table(props, edges):
    out = {}
    for g in props['general']:
        if not g['name'].startswith('zz_'): continue
        d = g['name'][3:]; hits = sorted({(min(int(d[:k]), int(d[k:])), max(int(d[:k]), int(d[k:]))) for k in range(1, len(d)) if all(x == '0' or x[0] != '0' for x in (d[:k], d[k:])) and (min(int(d[:k]), int(d[k:])), max(int(d[:k]), int(d[k:]))) in edges})
        if len(hits) == 1: out[hits[0]] = abs(g['value']) * 1e9
    return out
def load(chip):
    rec = json.load(open(f'{D}/job_direct_{chip}.json')); z = np.load(f'{D}/raw_direct_{chip}.npz'); return rec, z
def phase(rec, z, arm, tau, copies, breg, bval):
    m = len(rec['triples']); num = {s: np.zeros(m) for s in 'XY'}; den = {s: np.zeros(m) for s in 'XY'}
    for i, mt in enumerate(rec['meta']):
        if mt['arm'] != arm or mt['tau_us'] != tau or mt['copy'] not in copies: continue
        a = bits(z[f"a_{i}"], m).astype(float); bb = bits(z[f'{breg}_{i}'], m); w = (bb == bval)
        num[mt['setting']] += ((1 - 2 * a) * w).sum(0); den[mt['setting']] += w.sum(0)
    x = num['X'] / np.maximum(den['X'], 1); y = num['Y'] / np.maximum(den['Y'], 1)
    return np.arctan2(y, x), np.hypot(x, y), np.minimum(den['X'], den['Y'])
def shifts(rec, z, copies):
    taus = sorted({m['tau_us'] for m in rec['meta']}); T = np.array(taus) * 1e-6; out = {}
    for arm, breg in (('O1', 'b1'), ('O2', 'b2')):
        po = [phase(rec, z, arm, t, copies, breg, 1) for t in taus]; pz = [phase(rec, z, 'Z', t, copies, breg, 0) for t in taus]
        for j in range(len(rec['triples'])):
            co = np.array([p[1][j] for p in po]); cz = np.array([p[1][j] for p in pz]); no = np.array([p[2][j] for p in po]); nz = np.array([p[2][j] for p in pz])
            live = (co > 5 / np.sqrt(np.maximum(no, 1))) & (cz > 5 / np.sqrt(np.maximum(nz, 1)))
            if live.sum() < 5: out[(j, arm)] = np.nan; continue
            d = np.unwrap(np.array([p[0][j] for p in po])[live] - np.array([p[0][j] for p in pz])[live])
            out[(j, arm)] = abs(np.polyfit(T[live], d, 1)[0] / (2 * np.pi))
    return out
def table(rec, sh, zz):
    rows = []
    for j, (b1, a, b2) in enumerate(rec['triples']):
        for arm, b in (('O1', b1), ('O2', b2)): rows.append((j, zz.get((min(a, b), max(a, b)), np.nan), sh[(j, arm)], a, b))
    return rows
def stats(rows, idx):
    r = [rows[k] for k in idx]; st = np.array([x[1] for x in r]); df = np.array([x[2] for x in r]); ok = np.isfinite(st) & np.isfinite(df) & (st > 1e3) & (df > 0)
    st, df = st[ok], df[ok]
    if len(st) < 5: return np.nan, np.nan, np.nan
    rk = np.corrcoef(np.argsort(np.argsort(st)), np.argsort(np.argsort(df)))[0, 1]
    return np.polyfit(np.log(st), np.log(df), 1)[0], np.median(df / st), rk
def boot(rows, ntrip):
    by = {}
    for k, x in enumerate(rows): by.setdefault(x[0], []).append(k)
    trips = list(by); out = []
    for _ in range(B):
        pick = RNG.choice(trips, len(trips)); idx = [k for t in pick for k in by[t]]; out.append(stats(rows, idx))
    return np.array(out)
def lab(lo, hi, ref): return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'
def analyse(chip, props_path=None):
    rec, z = load(chip); props = json.load(open(props_path or f'{D}/props_direct_{chip}.json')); edges = {tuple(e) for e in rec['edges']}; zz = zz_table(props, edges)
    rows = table(rec, shifts(rec, z, (0, 1)), zz); r0 = table(rec, shifts(rec, z, (0,)), zz); r1 = table(rec, shifts(rec, z, (1,)), zz)
    pt = stats(rows, range(len(rows))); bs = boot(rows, len(rec['triples'])); lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
    g0 = boot(r0, 0); g1 = boot(r1, 0); gg = np.nanpercentile(g0 - g1, [2.5, 97.5], axis=0); ct = [(gg[0, k] > 0 or gg[1, k] < 0) for k in range(2)]
    st = np.array([x[1] for x in rows]); df = np.array([x[2] for x in rows])
    print(f'=== {chip}  job {rec["job_id"]}  edges {len(rows)}, fitted {int(np.isfinite(df).sum())}')
    print(f'  D1 median |df| {np.nanmedian(df)/1e3:.2f} kHz; median stated {np.nanmedian(st)/1e3:.2f} kHz')
    print(f'  D2 log-log slope {pt[0]:.3f} [{lo[0]:.3f}, {hi[0]:.3f}] vs 1 -> {"CANNOT_TELL (copies disagree)" if ct[0] else lab(lo[0], hi[0], 1)}; vs 0 -> {lab(lo[0], hi[0], 0)}')
    print(f'  D3 median |df|/stated {pt[1]:.3f} [{lo[1]:.3f}, {hi[1]:.3f}] vs 1 -> {"CANNOT_TELL (copies disagree)" if ct[1] else lab(lo[1], hi[1], 1)}')
    print(f'  D4 rank correlation {pt[2]:.3f} [{lo[2]:.3f}, {hi[2]:.3f}] vs 0 -> {lab(lo[2], hi[2], 0)}')
    zero = [(a, b, d) for (_, s, d, a, b) in rows if np.isfinite(s) and s <= 1e3]
    print('  D5 edges stated <= 1 kHz (A, B, |df| kHz):', ' '.join(f'({a},{b},{d/1e3:.1f})' for a, b, d in zero) or 'none')
    print('  per edge (A, B, stated kHz, |df| kHz):', ' '.join(f'({a},{b},{s/1e3:.1f},{d/1e3:.1f})' for (_, s, d, a, b) in rows))
    return pt, lo, hi
if __name__ == '__main__':
    for c in [x for x in sys.argv[1:] if x.startswith('ibm_') or x.startswith('synth')]: analyse(c)
