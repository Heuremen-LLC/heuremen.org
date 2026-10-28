"""combine_t1.py -- PLAN_LADDER_R2.md rule, written before run 2 was read. Marrakesh delivered/stated T1: each run alone (its own record),
the two runs pooled by qubit (P(1) per wait averaged over runs, one fit, stated = run 1's record), and chi2/dof between the run values."""
import json, sys, numpy as np
sys.path.insert(0, '.'); import ladder as L
def run(sfx):
    rec = json.load(open(f'job_ladder_ibm_marrakesh{sfx}.json')); z = np.load(f'raw_ladder_ibm_marrakesh{sfx}.npz')
    st = L.stated_of(json.load(open(f'props_ladder_ibm_marrakesh{sfx}.json')))
    cv = L.curves(rec['meta'], lambda i: z[f'c_{i}'][:, :len(rec['sets'][rec['meta'][i]['set']])], rec['sets'])
    return {q: cv[q]['T1']['pool'] for q in cv}, st
rng = np.random.default_rng(1); R = [run(''), run('_r2')]; vals = []
for k, (c, st) in enumerate(R):
    r = []
    for q, sig in c.items():
        T, A, _ = L.fit(L.TAUS['T1'], [sig[t] for t in L.TAUS['T1']], True)
        if A > 0.4 and np.isfinite(T) and T < 1990 and q in st['T1']: r.append(T / st['T1'][q])
    m, lo, hi, n = L.med_ci(r, rng); vals.append((m, (hi - lo) / 3.92, L.four(lo, hi))); print(f'run {k + 1}: T1 delivered/stated median {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {L.four(lo, hi)}')
st = R[0][1]; r = []
for q in R[0][0]:
    if q in R[1][0] and q in st['T1']:
        y = [(R[0][0][q][t] + R[1][0][q][t]) / 2 for t in L.TAUS['T1']]; T, A, _ = L.fit(L.TAUS['T1'], y, True)
        if A > 0.4 and np.isfinite(T) and T < 1990: r.append(T / st['T1'][q])
m, lo, hi, n = L.med_ci(r, rng); w = np.array([1 / v[1] ** 2 for v in vals]); mu = np.sum(w * [v[0] for v in vals]) / w.sum()
chi = float(np.sum(w * (np.array([v[0] for v in vals]) - mu) ** 2))
print(f'pooled by qubit (stated from run 1): {m:.3f} [{lo:.3f}, {hi:.3f}] n {n} -> {L.four(lo, hi)}; chi2/dof between runs {chi:.2f}')
r2_below = vals[1][2] == 'BELOW'
print('BILL LINE (run 2 alone BELOW, pooled BELOW, chi2/dof < 2):', 'YES' if (r2_below and hi < 1 and chi < 2) else 'NO')
