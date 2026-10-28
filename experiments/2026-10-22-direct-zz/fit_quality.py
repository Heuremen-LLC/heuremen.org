"""fit_quality.py -- defined after the data (outside referee, 23 Oct): residual rms of each edge's linear phase fit, with the
committed analysis's own phase/unwrap code; edges whose fit residual exceeds 0.5 rad (unwrap slip) are dropped and D2-D4 recomputed
with the committed stats/bootstrap. The same rule on both chips. Not in the plan."""
import sys, json, numpy as np
sys.argv = ['x', 'ibm_phoenix']
import analyze_direct as A
def rms_table(rec, z):
    taus = sorted({m['tau_us'] for m in rec['meta']}); T = np.array(taus) * 1e-6; out = {}
    for arm, breg in (('O1', 'b1'), ('O2', 'b2')):
        po = [A.phase(rec, z, arm, t, (0, 1), breg, 1) for t in taus]; pz = [A.phase(rec, z, 'Z', t, (0, 1), breg, 0) for t in taus]
        for j in range(len(rec['triples'])):
            co = np.array([p[1][j] for p in po]); cz = np.array([p[1][j] for p in pz]); no = np.array([p[2][j] for p in po]); nz = np.array([p[2][j] for p in pz])
            live = (co > 5 / np.sqrt(np.maximum(no, 1))) & (cz > 5 / np.sqrt(np.maximum(nz, 1)))
            if live.sum() < 5: out[(j, arm)] = np.nan; continue
            d = np.unwrap(np.array([p[0][j] for p in po])[live] - np.array([p[0][j] for p in pz])[live]); c = np.polyfit(T[live], d, 1)
            out[(j, arm)] = float(np.sqrt(np.mean((d - np.polyval(c, T[live])) ** 2)))
    return out
for chip, d in (('ibm_fez', '/Users/heuremen/bones-forge/heuremen.org/experiments/2026-10-22-direct-zz'), ('ibm_phoenix', A.D)):
    A.D = d; rec, z = A.load(chip); props = json.load(open(f'{d}/props_direct_{chip}.json')); zz = A.zz_table(props, {tuple(e) for e in rec['edges']})
    rows = A.table(rec, A.shifts(rec, z, (0, 1)), zz); rm = rms_table(rec, z)
    keys = [(j, arm) for j in range(len(rec['triples'])) for arm in ('O1', 'O2')]
    bad = [k for k, (key) in enumerate(keys) if not (rm[key] <= 0.5)]
    print(f'=== {chip}: median fit rms {np.nanmedian(list(rm.values())):.3f} rad; edges with rms > 0.5 rad or no fit: ' + ' '.join(f'({rows[k][3]},{rows[k][4]},rms {rm[keys[k]]:.2f})' for k in bad))
    keep = [r for k, r in enumerate(rows) if k not in bad]; A.RNG = np.random.default_rng(7)
    pt = A.stats(keep, range(len(keep))); bs = A.boot(keep, 0); lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
    n = sum(1 for r in keep if np.isfinite(r[1]) and r[1] > 1e3 and np.isfinite(r[2]) and r[2] > 0)
    print(f'  without them, {n} edges: slope {pt[0]:.3f} [{lo[0]:.3f}, {hi[0]:.3f}]; median |df|/stated {pt[1]:.3f} [{lo[1]:.3f}, {hi[1]:.3f}]; rank {pt[2]:.3f} [{lo[2]:.3f}, {hi[2]:.3f}]')
