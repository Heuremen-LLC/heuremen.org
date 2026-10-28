"""hahn_pool.py CHIP RAW1 JOB1 [RAW2 JOB2 ...] -- ADDENDUM after the referee (22 Oct night), defined after runs 1-2 were read, before any
fresh-record run is read. Runs of one chip are pooled at the QUBIT level (per qubit, P(1) per tau averaged over runs, equal shots), never
combined by inverse variance (same qubits, same stated values: not independent). Two fits per qubit:
  FIXED  P1 = 1/2 + A exp(-t/T2)            (the committed hahn.py fit)
  FREE   P1 = B + A exp(-t/T2), B free     (the model of IBM's open-source T2HahnAnalysis: base free, bounded [0, 1])
Reading: median delivered/stated over qubits, 4000-draw bootstrap over qubits, four outcomes; stated T2 from the FIRST run's record."""
import sys, json, numpy as np
sys.path.insert(0, '.'); import hahn as R
GRID = np.exp(np.linspace(np.log(1), np.log(2000), 1200))
def fit(t, y, free):
    t = np.asarray(t, float); y = np.asarray(y, float); best = (np.nan, np.nan, np.inf)
    for T2 in GRID:
        e = np.exp(-t / T2)
        if free:
            X = np.column_stack([np.ones_like(e), e]); c = np.linalg.lstsq(X, y, rcond=None)[0]; r = y - X @ c; A = c[1]
        else:
            yy = y - 0.5; A = float(e @ yy / (e @ e)); r = yy - A * e
        s = float(np.sqrt(np.mean(r ** 2)))
        if s < best[2]: best = (float(T2), float(A), s)
    return best
def p1_of(raw, job):
    rec = json.load(open(job)); z = np.load(raw); out = {}
    for i, m in enumerate(rec['meta']):
        S = rec['sets'][m['set']]; b = z[f'c_{i}'][:, :len(S)]
        for j, q in enumerate(S): out.setdefault(q, {}).setdefault(m['tau_us'], []).append(float(b[:, j].mean()))
    return out
if __name__ == '__main__':
    chip = sys.argv[1]; pairs = list(zip(sys.argv[2::2], sys.argv[3::2]))
    props = json.load(open(pairs[0][1].replace('job_hahn', 'props_hahn')))
    stated = {q: [x['value'] for x in ps if x['name'] == 'T2'][0] for q, ps in enumerate(props['qubits']) if any(x['name'] == 'T2' for x in ps)}
    runs = [p1_of(r, j) for r, j in pairs]; qs = sorted(set.intersection(*[set(r) for r in runs]))
    rng = np.random.default_rng(1)
    for free in (False, True):
        x = []
        for q in qs:
            y = [np.mean([np.mean(r[q][t]) for r in runs]) for t in R.TAUS]; T2, A, _ = fit(R.TAUS, y, free)
            if A > 0.2 and q in stated and np.isfinite(T2): x.append(T2 / stated[q])
        x = np.array(x); med = float(np.median(x)); lo, hi = np.percentile([np.median(rng.choice(x, len(x))) for _ in range(4000)], [2.5, 97.5])
        out = 'ABOVE' if lo > 1 else 'BELOW' if hi < 1 else 'CONSISTENT'
        print(f'{chip} {len(pairs)} run(s) pooled by qubit, {"FREE asymptote" if free else "FIXED 1/2   "}: median {med:.3f} [{lo:.3f}, {hi:.3f}] over {len(x)} qubits -> {out}  (below stated {int((x < 1).sum())})')
