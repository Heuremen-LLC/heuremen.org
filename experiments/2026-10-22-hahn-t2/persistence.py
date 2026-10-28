"""persistence.py -- are the same qubits below their stated T2 every run, or different ones each time? (23 Oct 2026, business analysis,
written after the runs were read.) Uses only files on disk; zero QPU seconds. Fits are hahn_pool.fit (FIXED and FREE), unchanged.

Per chip with repeat runs, per qubit:
  WITHIN  : copy 1 vs copy 2 of the SAME job (delivered T2). Its agreement is the ceiling set by measurement noise alone.
  BETWEEN : run a vs run b (delivered T2, each run's copies pooled). If BETWEEN ~ WITHIN, each qubit keeps its T2 between jobs;
            if BETWEEN << WITHIN, qubits move between jobs by more than the measurement can explain.
  RATIO   : delivered/stated per run (each run's own record); agreement of ratios between runs = persistence of "below stated".
  TRAP    : Spearman of ratio against the stated value; strongly negative means high stated values read low partly by luck
            (a stated T2 that is itself a noisy single measurement), which a per-qubit bill line must account for.
  BELOW-ALL: qubits below stated in every run, against the count expected if runs were independent.
Agreement = Spearman rank correlation, 95 % interval from 2000 bootstrap draws over qubits."""
import json, sys, itertools, numpy as np
sys.path.insert(0, '.'); import hahn as R
from hahn_pool import fit
RUNS = {'ibm_kingston': ['kingston_run1', 'kingston_run2', 'kingston_fresh'],
        'ibm_marrakesh': ['marrakesh_run1', 'marrakesh_run2'],
        'ibm_fez': ['fez_run1', 'fez_run2', 'fez_clock1', 'fez_clock2']}
rng = np.random.default_rng(7)

def rank(a): return np.argsort(np.argsort(a)).astype(float)
def spear(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 8: return (np.nan, np.nan, np.nan)
    r = np.corrcoef(rank(a), rank(b))[0, 1]; bs = []
    for _ in range(2000):
        i = rng.integers(0, len(a), len(a)); bs.append(np.corrcoef(rank(a[i]), rank(b[i]))[0, 1])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5]); return (r, lo, hi)
def fmt(t): return 'n/a' if np.isnan(t[0]) else f'{t[0]:+.2f} [{t[1]:+.2f}, {t[2]:+.2f}]'

def per_qubit(tag, free):
    """delivered T2 per qubit: pooled, copy 1, copy 2; and the stated T2 from this run's own record."""
    rec = json.load(open(f'job_hahn_ibm_{tag}.json')); z = np.load(f'raw_hahn_ibm_{tag}.npz')
    props = json.load(open(f'props_hahn_ibm_{tag}.json'))
    stated = {q: [x['value'] for x in ps if x['name'] == 'T2'][0] for q, ps in enumerate(props['qubits']) if any(x['name'] == 'T2' for x in ps)}
    p = {}
    for i, m in enumerate(rec['meta']):
        S = rec['sets'][m['set']]; b = z[f'c_{i}'][:, :len(S)]
        for j, q in enumerate(S): p.setdefault(q, {}).setdefault(m['copy'], {})[m['tau_us']] = float(b[:, j].mean())
    out = {}
    for q, cp in p.items():
        if q not in stated or len(cp) < 2: continue
        vals = {}
        for c in cp:
            T2, A, _ = fit(R.TAUS, [cp[c][t] for t in R.TAUS], free); vals[c] = T2 if A > 0.2 else np.nan
        y = [np.mean([cp[c][t] for c in cp]) for t in R.TAUS]; T2, A, _ = fit(R.TAUS, y, free)
        c0, c1 = sorted(cp)[:2]
        out[q] = dict(pooled=T2 if A > 0.2 else np.nan, c0=vals[c0], c1=vals[c1], stated=stated[q])
    return out

for free in (False, True):
    print(f'\n######## {"FREE asymptote" if free else "FIXED 1/2"} fit ########')
    for chip, tags in RUNS.items():
        d = {t: per_qubit(t, free) for t in tags}
        print(f'\n== {chip}: runs {", ".join(tags)}')
        for t in tags:
            qs = [q for q in d[t] if np.isfinite(d[t][q]['c0']) and np.isfinite(d[t][q]['c1'])]
            w = spear([d[t][q]['c0'] for q in qs], [d[t][q]['c1'] for q in qs])
            qr = [q for q in d[t] if np.isfinite(d[t][q]['pooled'])]
            trap = spear([d[t][q]['pooled'] / d[t][q]['stated'] for q in qr], [d[t][q]['stated'] for q in qr])
            print(f'  WITHIN {t:16} copy vs copy (delivered T2): {fmt(w)}  n={len(qs)} | TRAP ratio vs stated: {fmt(trap)}')
        for a, b in itertools.combinations(tags, 2):
            qs = [q for q in d[a] if q in d[b] and np.isfinite(d[a][q]['pooled']) and np.isfinite(d[b][q]['pooled'])]
            bt = spear([d[a][q]['pooled'] for q in qs], [d[b][q]['pooled'] for q in qs])
            br = spear([d[a][q]['pooled'] / d[a][q]['stated'] for q in qs], [d[b][q]['pooled'] / d[b][q]['stated'] for q in qs])
            print(f'  BETWEEN {a} vs {b}: delivered T2 {fmt(bt)} | ratio {fmt(br)}  n={len(qs)}')
        qs = [q for q in d[tags[0]] if all(q in d[t] and np.isfinite(d[t][q]['pooled']) for t in tags)]
        below = [[d[t][q]['pooled'] < d[t][q]['stated'] for q in qs] for t in tags]
        allb = int(np.sum(np.all(below, axis=0))); exp = len(qs) * np.prod([np.mean(b) for b in below])
        print(f'  BELOW-ALL: {allb} of {len(qs)} qubits below stated in every run; {exp:.1f} expected if runs were independent')
