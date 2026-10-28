"""trap.py -- do the qubits IBM states as best deliver least, relative to their stated T2? (23 Oct 2026, business analysis, after the data.)
Per run: qubits split into quarters by STATED T2 (that run's own record). Reading: median delivered/stated in the top quarter minus the
bottom quarter (top - bottom; negative = the qubits stated best under-deliver most), 4000-draw bootstrap over qubits, both fits.
Rebuttal tests built in:
  FIT BIAS : the same analysis on synthetic qubits whose true T2 IS the stated value (P1 = 1/2 + A exp(-t/T2), A from the qubit's
             own fit, binomial shots matching the job: 2 copies x shots), through the same fit. Its median and 95 % range is the null.
  DRIFT    : the kingston fresh-record run (stated values 0.6 h old) is read on its own.
Outcome: BELOW if the 95 % interval of (top - bottom) is under 0 AND under the synthetic null's 2.5 % point, under BOTH fits."""
import json, sys, numpy as np
sys.path.insert(0, '.'); import hahn as R
from hahn_pool import fit
RUNS = ['kingston_run1', 'kingston_run2', 'kingston_fresh', 'marrakesh_run1', 'marrakesh_run2', 'fez_run1', 'fez_run2',
        'fez_clock1', 'fez_clock2', 'pittsburgh', 'phoenix', 'boston']
TAUS = np.array(R.TAUS, float)

def load(tag):
    rec = json.load(open(f'job_hahn_ibm_{tag}.json')); z = np.load(f'raw_hahn_ibm_{tag}.npz'); props = json.load(open(f'props_hahn_ibm_{tag}.json'))
    stated = {q: [x['value'] for x in ps if x['name'] == 'T2'][0] for q, ps in enumerate(props['qubits']) if any(x['name'] == 'T2' for x in ps)}
    p, n = {}, {}
    for i, m in enumerate(rec['meta']):
        S = rec['sets'][m['set']]; b = z[f'c_{i}'][:, :len(S)]
        for j, q in enumerate(S):
            p.setdefault(q, {}).setdefault(m['tau_us'], []).append(float(b[:, j].mean())); n.setdefault(q, {}).setdefault(m['tau_us'], []).append(b.shape[0])
    return p, n, stated

def quart(ratio, stated, rng, B=4000):
    q = np.array(sorted(ratio)); r = np.array([ratio[k] for k in q]); s = np.array([stated[k] for k in q])
    lo_cut, hi_cut = np.percentile(s, [25, 75]); top, bot = r[s >= hi_cut], r[s <= lo_cut]
    d = np.median(top) - np.median(bot)
    bs = [np.median(rng.choice(top, len(top))) - np.median(rng.choice(bot, len(bot))) for _ in range(B)]
    return d, np.percentile(bs, 2.5), np.percentile(bs, 97.5), np.median(top), np.median(bot), len(top), len(bot)

def run(tag, free, rng, nsim=40):
    p, n, stated = load(tag); ratio, amp = {}, {}
    for q in p:
        if q not in stated: continue
        y = [np.mean(p[q][t]) for t in R.TAUS]; T2, A, _ = fit(R.TAUS, y, free)
        if A > 0.2 and np.isfinite(T2): ratio[q] = T2 / stated[q]; amp[q] = A
    real = quart(ratio, stated, rng)
    nulls = []
    for _ in range(nsim):                       # synthetic: true T2 = stated, same shots, same fit
        sr = {}
        for q in ratio:
            shots = np.array([sum(n[q][t]) for t in R.TAUS]); pr = np.clip(0.5 + amp[q] * np.exp(-TAUS / stated[q]), 0, 1)
            y = rng.binomial(shots, pr) / shots; T2, A, _ = fit(R.TAUS, y, free)
            if A > 0.2 and np.isfinite(T2): sr[q] = T2 / stated[q]
        nulls.append(quart(sr, stated, rng, B=1)[0])
    return real, np.percentile(nulls, [2.5, 50, 97.5]), len(ratio)

if __name__ == '__main__':
    rng = np.random.default_rng(11)
    for tag in RUNS:
        line = []
        for free in (False, True):
            (d, lo, hi, mt, mb, nt, nb), (n025, n50, n975), nq = run(tag, free, rng)
            out = 'BELOW' if (hi < 0 and hi < n025) else 'ABOVE' if (lo > 0 and lo > n975) else 'CONSISTENT'
            line.append(f'{"FREE " if free else "FIXED"} top-bottom {d:+.3f} [{lo:+.3f}, {hi:+.3f}] (top q {mt:.3f}, bottom q {mb:.3f}; n {nt}/{nb}) '
                        f'null {n50:+.3f} [{n025:+.3f}, {n975:+.3f}] -> {out}')
        print(f'{tag}  ({nq} qubits)\n   ' + '\n   '.join(line), flush=True)
