#!/usr/bin/env python3
"""sim_gate_nonexp.py -- how large a negative slope per gate can NON-EXPONENTIAL decay produce with PERFECT (zero-cost) X gates?

Reads gate_cost/data/gate_cost.json (read-only). Needs numpy: /Users/heuremen/bones-forge/venv_quantum/bin/python sim_gate_nonexp.py
Estimator as in verify_gate_cost.py:
  circuit  : start |0>, n x [ wait 30 us ; X ]  (n = 2,6,10), read right after the last X  (so n-1 excited waits)
  reference: start |1> (X), wait 30 / 90 / 300 us, read
  q        : own = clip(R30); r90 = clip(R90)^(1/3); r300 = clip(R300)^(1/10), clip to [0.02,1]
  expected : recursion on<-off+on*(1-q), off<-on*q  from on=0
  residual = N_n - expected(q, n), mean over readable qubits; slope = least-squares slope over n=2,6,10.
Synthetic qubits reuse each real run's readable-qubit list, its R300 as per-qubit anchor, and its per-level shot errors.
"""
import json, sys, os
import numpy as np

DATA = '/Users/heuremen/bones-forge/heuremen.org/experiments/2026-10-21-measurements/gate_cost/data/gate_cost.json'
NS = (2, 6, 10)
RUNS = json.load(open(DATA))['runs']
SEL = [371, 392, 721, 722, 723, 441, 711, 712, 713]          # fez, kingston x4, boston x4 (the nine runs analysed)
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(60); GH_W = GH_W / GH_W.sum()   # standard-normal quadrature

def book(q, n):
    on, off = np.zeros_like(q), np.ones_like(q)
    for _ in range(n): on, off = off + on * (1 - q), on * q
    return on

def book_var(qs, n):           # per-wait survival list qs[0..n-1] (wait 1 acts on ground state, as in the real circuit)
    on, off = 0.0, 1.0
    for i in range(n): on, off = off + on * (1 - qs[i]), on * qs[i]
    return on

def clipq(x): return np.clip(x, 0.02, 1.0)

def slope_est(N, R30, R90, R300):
    """N: dict n->(...,Q) levels; R*: (...,Q). returns dict ruler -> slope per gate (mean over qubits first)"""
    out = {}
    for nm, q in (('own', clipq(R30)), ('r90', clipq(R90) ** (1 / 3)), ('r300', clipq(R300) ** (1 / 10))):
        q = clipq(q)
        d = np.stack([np.nanmean(N[n] - book(q, n), axis=-1) for n in NS], -1)         # (..., 3)
        x = np.array(NS, float); x = x - x.mean()
        out[nm] = (d * x).sum(-1) / (x ** 2).sum()
    return out

# ---------------- decay models: each returns true levels for ON refs and the three train levels ----------------
def solve(fun, a, lo=1e-9, hi=1.0):
    """per-qubit scale parameter k so that fun(k) = a (fun decreasing in k)  -- bisection in log k"""
    lo = np.full_like(a, np.log(1e-7)); hi = np.full_like(a, np.log(1.0))
    for _ in range(70):
        m = (lo + hi) / 2; f = fun(np.exp(m)); big = f > a
        lo = np.where(big, m, lo); hi = np.where(big, hi, m)
    return np.exp((lo + hi) / 2)

def model_exp(a, gate=0.0, jit=None, rng=None):
    k = -np.log(a) / 300.0
    return mix_core(k, 0.0, a, gate, jit, rng)

def mix_core(k0, sigma, a, gate, jit, rng):
    """per-shot rate k = k0*exp(sigma*z - sigma^2/2); rate persists through the whole circuit. gate = X duration (us) added to each excited segment."""
    z = GH_X if sigma > 0 else np.zeros(1); w = GH_W if sigma > 0 else np.ones(1)
    k = k0[None, :] * np.exp(sigma * z[:, None] - 0.5 * sigma ** 2)           # (G,Q)
    S = lambda t: (w[:, None] * np.exp(-k * t)).sum(0)
    q = np.exp(-k * (30.0 + gate))
    tr = {n: (w[:, None] * book(q, n)).sum(0) for n in NS}
    return dict(R30=S(30.0), R90=S(90.0), R300=S(300.0), **{'N%d' % n: tr[n] for n in NS})

def model_mix(a, sigma, gate=0.0):
    f = lambda k0: (GH_W[:, None] * np.exp(-k0[None, :] * np.exp(sigma * GH_X[:, None] - 0.5 * sigma ** 2) * 300.0)).sum(0)
    k0 = solve(f, a)
    return mix_core(k0, sigma, a, gate, None, None)

def model_stretch(a, beta, clock):
    """Kohlrausch S(t)=exp(-(t/tau)^beta), tau per qubit so S(300)=a.  clock:
       'abs0'   hazard keyed to time since circuit start (reference X at t=0; train X's at 30,60,..)
       'abs30'  hazard keyed to time since the FIRST X (train clock offset by 30 us relative to 'abs0')
       'renew'  decay clock restarts at every X (each excited wait has survival S(30): pure Markov)"""
    beta = np.broadcast_to(np.asarray(beta, float), a.shape)
    tau = 300.0 / (-np.log(a)) ** (1.0 / beta)
    S = lambda t: np.exp(-(t / tau) ** beta)
    out = dict(R30=S(30.0), R90=S(90.0), R300=S(300.0))
    for n in NS:
        if clock == 'renew': qs = [S(30.0)] * n
        else:
            c = 0.0 if clock == 'abs0' else 30.0
            # wait i (1-based) is on ground for i=1; for i>=2 excited from X at time 30*(i-1) until 30*i
            qs = [None] + [S(30.0 * i - c) / S(30.0 * (i - 1) - c) for i in range(2, n + 1)]
            qs[0] = np.ones_like(a)
        out['N%d' % n] = book_var(qs, n)
    return out

def model_mix_redraw(a, sigma):
    """CONTROL: rate is re-drawn independently at EVERY wait (fluctuation faster than 30 us): Markov with q = mean survival. Own interval must be immune."""
    f = lambda k0: (GH_W[:, None] * np.exp(-k0[None, :] * np.exp(sigma * GH_X[:, None] - 0.5 * sigma ** 2) * 300.0)).sum(0)
    k0 = solve(f, a); k = k0[None, :] * np.exp(sigma * GH_X[:, None] - 0.5 * sigma ** 2)
    S = lambda t: (GH_W[:, None] * np.exp(-k * t)).sum(0)
    q = S(30.0)
    return dict(R30=q, R90=S(90.0), R300=S(300.0), **{'N%d' % n: book(q, n) for n in NS})

def jitter(truth_fn, a, sig_j, rng, S):
    """per-circuit T1 drift: every circuit (N2,N6,N10,R30,R90,R300) gets its own rate multiplier exp(sig_j*z) -- exponential decay otherwise"""
    k = -np.log(a) / 300.0
    out = {}
    for key in ('R30', 'R90', 'R300', 'N2', 'N6', 'N10'):
        kk = k[None, :] * np.exp(sig_j * rng.standard_normal((S, len(a))) - 0.5 * sig_j ** 2)
        if key[0] == 'R': out[key] = np.exp(-kk * float(key[1:]))
        else: out[key] = book(np.exp(-kk * 30.0), int(key[1:]))
    return out

# ---------------- run machinery ----------------
def run_arrays(r):
    A = {k: np.array(r[k]) for k in ('N2', 'N6', 'N10', 'R30', 'R90', 'R300', 'se_N2', 'se_N6', 'se_N10', 'se_R30', 'se_R90', 'se_R300')}
    return A

NOISE = [1.0]

def synth_levels(truth, A, S, rng):
    """truth: dict of (Q,) or (S,Q) arrays; add Gaussian shot noise with the real per-qubit-per-level errors"""
    L = {}
    for k in ('N2', 'N6', 'N10', 'R30', 'R90', 'R300'):
        t = np.asarray(truth[k]); t = np.broadcast_to(t, (S, len(A[k])))
        L[k] = t + NOISE[0] * A['se_' + k][None, :] * rng.standard_normal((S, len(A[k])))
    return L

def evaluate(L):
    N = {2: L['N2'], 6: L['N6'], 10: L['N10']}
    sl = slope_est(N, L['R30'], L['R90'], L['R300'])
    share = (L['R30'] < clipq(L['R300']) ** 0.1).mean(-1)
    return sl, share, L['R30'].mean(-1), L['R90'].mean(-1)

def simulate(label, build, nseeds=300, seed=1, runs=SEL, jit=False):
    """build(a, A, rng, S) -> truth dict. returns dict ruler -> (mean over runs of slope, sd over seeds of run-mean), shares."""
    rng = np.random.default_rng(seed)
    per = {k: [] for k in ('own', 'r90', 'r300')}; sh = []; m30 = []; m90 = []
    for r in RUNS:
        if r['run'] not in runs: continue
        A = run_arrays(r); a = np.clip(A['R300'], 0.02, 0.98)
        truth = build(a, A, rng, nseeds)
        L = synth_levels(truth, A, nseeds, rng)
        sl, share, a30, a90 = evaluate(L)
        for k in per: per[k].append(sl[k])
        sh.append(share); m30.append(a30); m90.append(a90)
    res = {}
    for k in per:
        M = np.mean(per[k], axis=0)      # (seeds,) average over the runs
        res[k] = (M.mean(), M.std(), np.mean([p.std() for p in per[k]]))
    return res, np.mean(sh), np.mean(m30), np.mean(m90)

def real_numbers(runs=SEL):
    per = {k: [] for k in ('own', 'r90', 'r300')}; sh = []; m30 = []; m90 = []
    for r in RUNS:
        if r['run'] not in runs: continue
        A = run_arrays(r)
        sl, share, a30, a90 = evaluate({k: A[k][None, :] for k in ('N2', 'N6', 'N10', 'R30', 'R90', 'R300')})
        for k in per: per[k].append(float(sl[k][0]))
        sh.append(float(share[0])); m30.append(float(a30[0])); m90.append(float(a90[0]))
    return {k: np.mean(v) for k, v in per.items()}, np.mean(sh), np.mean(m30), np.mean(m90), per, runs

def fmt(res, sh, a30, a90):
    return ('own %+.5f (seed sd %.5f)  r90 %+.5f  r300 %+.5f | below-extrap share %.2f  mean R30 %.3f R90 %.3f' %
            (res['own'][0], res['own'][1], res['r90'][0], res['r300'][0], sh, a30, a90))

def match_runs(P, name, make, grid, nseeds=100):
    """for each real run find the model parameter whose synthetic below-extrapolation share equals the run's real share, then report the model slope there"""
    P('%s -- parameter matched PER RUN so the synthetic share below the 300-us extrapolation equals that run\'s real share:' % name)
    P('   run  chip        param | real share | real slope own/r90/r300            | model slope own/r90/r300         | mean N2,N6,N10 real vs model')
    rows = []
    for r in RUNS:
        if r['run'] not in SEL: continue
        A = run_arrays(r); a = np.clip(A['R300'], 0.02, 0.98)
        real = evaluate({k: A[k][None, :] for k in ('N2', 'N6', 'N10', 'R30', 'R90', 'R300')})
        rs = float(real[1][0]); rsl = {k: float(v[0]) for k, v in real[0].items()}
        rng = np.random.default_rng(100 + r['run']); sh = []; sl = []
        for g in grid:
            tr = make(a, g); L = synth_levels(tr, A, nseeds, rng); e = evaluate(L)
            sh.append(float(np.mean(e[1]))); sl.append({k: float(np.mean(v)) for k, v in e[0].items()})
        sh = np.array(sh); i = int(np.argmin(np.abs(sh - rs)))
        j = i
        for c in (i - 1, i + 1):
            if 0 <= c < len(sh) and min(sh[i], sh[c]) <= rs <= max(sh[i], sh[c]) and sh[c] != sh[i]: j = c
        w = 0 if j == i else (rs - sh[i]) / (sh[j] - sh[i])
        gp = grid[i] + w * (grid[j] - grid[i]); ms = {k: sl[i][k] + w * (sl[j][k] - sl[i][k]) for k in sl[i]}
        tr = make(a, gp); mN = [float(np.mean(tr['N%d' % n])) for n in NS]; rN = [float(np.mean(A['N%d' % n])) for n in NS]
        P('   %3d  %-11s %.3f |   %.2f     | %+.5f %+.5f %+.5f | %+.5f %+.5f %+.5f | real %.4f %.4f %.4f  model %.4f %.4f %.4f' %
          (r['run'], r['chip'], gp, rs, rsl['own'], rsl['r90'], rsl['r300'], ms['own'], ms['r90'], ms['r300'], *rN, *mN))
        rows.append((rsl, ms))
    f = lambda i, k: np.mean([x[i][k] for x in rows])
    P('   MEAN of 9 runs: real own %+.5f r90 %+.5f r300 %+.5f   model own %+.5f r90 %+.5f r300 %+.5f   (model/real own = %.2f)' % (f(0,'own'), f(0,'r90'), f(0,'r300'), f(1,'own'), f(1,'r90'), f(1,'r300'), f(1,'own')/f(0,'own')))
    P()

def per_qubit_binned(P):
    """Test: a true gate cost should not depend on how non-exponential a qubit's decay is; the mixture effect should.
    Pool the 9 runs, per-qubit own-interval slope (3-point), binned by d = R300^(1/10) - R30 (positive = below its extrapolation)."""
    x = np.array(NS, float); x = x - x.mean()
    def pooled(get):
        D, Sl = [], []
        for r in RUNS:
            if r['run'] not in SEL: continue
            L = get(r)
            q = clipq(L['R30']); res = np.stack([L['N%d' % n] - book(q, n) for n in NS], -1)
            Sl.append((res * x).sum(-1) / (x ** 2).sum()); D.append(clipq(L['R300']) ** 0.1 - L['R30'])
        return np.concatenate(D), np.concatenate(Sl)
    edges = [-9, -0.03, 0.0, 0.03, 0.06, 0.10, 9]
    def table(D, Sl):
        return [(int(((D >= lo) & (D < hi)).sum()), float(Sl[(D >= lo) & (D < hi)].mean())) for lo, hi in zip(edges[:-1], edges[1:])]
    rng = np.random.default_rng(7)
    real = pooled(lambda r: {k: np.array(r[k]) for k in ('N2', 'N6', 'N10', 'R30', 'R300')})
    def synth(builder):
        def g(r):
            A = run_arrays(r); a = np.clip(A['R300'], 0.02, 0.98); out = {k: 0 for k in ('N2', 'N6', 'N10', 'R30', 'R300')}
            acc = None
            for _ in range(40):
                L = synth_levels(builder(a), A, 1, rng); L = {k: v[0] for k, v in L.items()}
                acc = L if acc is None else {k: np.concatenate([acc[k], L[k]]) for k in L}
            return acc
        return g
    # replicate per-qubit: 40 noise draws per qubit
    def pooled_s(builder):
        D, Sl = [], []
        for r in RUNS:
            if r['run'] not in SEL: continue
            A = run_arrays(r); a = np.clip(A['R300'], 0.02, 0.98); L = synth_levels(builder(a), A, 40, rng)
            q = clipq(L['R30']); res = np.stack([L['N%d' % n] - book(q, n) for n in NS], -1)
            Sl.append(((res * x).sum(-1) / (x ** 2).sum()).ravel()); D.append((clipq(L['R300']) ** 0.1 - L['R30']).ravel())
        return np.concatenate(D), np.concatenate(Sl)
    P('PER-QUBIT own-interval slope, binned by non-exponentiality d = R300^(1/10) - R30 (9 runs pooled). cells: n qubits, mean slope')
    P('   bins of d: ' + '  '.join('[%s,%s)' % (('-inf' if lo < -1 else '%.2f' % lo), ('inf' if hi > 1 else '%.2f' % hi)) for lo, hi in zip(edges[:-1], edges[1:])))
    for nm, (D, Sl) in (('MEASURED', real), ('synthetic exp (perfect X)', pooled_s(lambda a: model_exp(a))), ('synthetic mixture sigma=0.43', pooled_s(lambda a: model_mix(a, 0.43)))):
        P('   %-30s' % nm + '  '.join('%d:%+.4f' % c for c in table(D, Sl)))
    P()

def main2(P):
    NOISE[0] = 0.0
    P('NOISE-FREE (infinite shots, exact expectations; slope of the estimator on the true curves; exp baseline must be 0):')
    for nm, b in (('exp', lambda a, A, rng, S: model_exp(a)),
                  ('mixture sigma=0.4', lambda a, A, rng, S: model_mix(a, 0.4)),
                  ('mixture sigma=0.5', lambda a, A, rng, S: model_mix(a, 0.5)),
                  ('mixture sigma=1.2', lambda a, A, rng, S: model_mix(a, 1.2)),
                  ('CONTROL mixture re-drawn every wait sigma=0.5', lambda a, A, rng, S: model_mix_redraw(a, 0.5)),
                  ('CONTROL mixture re-drawn every wait sigma=1.2', lambda a, A, rng, S: model_mix_redraw(a, 1.2)),
                  ('stretch beta=0.9 renew', lambda a, A, rng, S: model_stretch(a, 0.9, 'renew')),
                  ('stretch beta=0.9 abs0', lambda a, A, rng, S: model_stretch(a, 0.9, 'abs0')),
                  ('exp gate g=0.035', lambda a, A, rng, S: model_exp(a, gate=0.035)),
                  ('exp gate g=0.1', lambda a, A, rng, S: model_exp(a, gate=0.1))):
        res, sh, a30, a90 = simulate(nm, b, nseeds=1); P('%-52s own %+.5f  r90 %+.5f  r300 %+.5f' % (nm, res['own'][0], res['r90'][0], res['r300'][0]))
    NOISE[0] = 1.0
    P()
    grid = [round(x, 3) for x in np.arange(0.0, 1.01, 0.05)]
    per_qubit_binned(P)
    match_runs(P, 'MIXTURE (T1 varies shot to shot, sigma = lognormal width of the rate)', lambda a, g: model_mix(a, max(g, 1e-9)), grid)
    match_runs(P, 'STRETCHED, clock since circuit start (beta<1 side only)', lambda a, g: model_stretch(a, g, 'abs0'), [round(x, 3) for x in np.arange(1.0, 0.79, -0.02)])
    match_runs(P, 'STRETCHED, clock restarted at each X (physically: each flip re-excites)', lambda a, g: model_stretch(a, g, 'renew'), [round(x, 3) for x in np.arange(1.0, 0.79, -0.02)])
    match_runs(P, 'CONTROL: mixture re-drawn at every wait', lambda a, g: model_mix_redraw(a, max(g, 1e-9)), grid)

def main():
    out = []
    P = lambda s='': (out.append(s), print(s))
    NSEED = 300
    rl, rsh, r30, r90, per, runs = real_numbers()
    P('MEASURED DATA (9 runs: ibm_fez 371, ibm_kingston 392/721/722/723, ibm_boston 441/711/712/713; mean of runs):')
    P('   own %+.5f  r90 %+.5f  r300 %+.5f | below-extrap share %.2f  mean R30 %.3f R90 %.3f' % (rl['own'], rl['r90'], rl['r300'], rsh, r30, r90))
    idx = {r['run']: i for i, r in enumerate([x for x in RUNS if x['run'] in SEL])}
    for r in [x for x in RUNS if x['run'] in SEL]:
        i = idx[r['run']]; P('   run %d %s own %+.5f r90 %+.5f r300 %+.5f' % (r['run'], r['chip'], per['own'][i], per['r90'][i], per['r300'][i]))
    P('   (note: mean over runs of single-run slopes)')
    P()
    P('SYNTHETIC: PERFECT X gates (zero cost). %d noise seeds per row; noise = real per-qubit shot errors; anchor = each qubit\'s R300.' % NSEED)
    P('Slope is mean over the 9 runs; "seed sd" is the sd of that 9-run mean over seeds.')
    P()
    def row(name, build, **kw):
        res, sh, a30, a90 = simulate(name, build, nseeds=NSEED, **kw); P('%-52s %s' % (name, fmt(res, sh, a30, a90))); return res
    P('(a) single exponential')
    row('exp', lambda a, A, rng, S: model_exp(a))
    P('(b) T1 mixture: rate varies shot-to-shot (lognormal sigma), persists through a whole circuit')
    for s in (0.3, 0.6, 0.9, 1.2, 1.5, 2.0):
        row('mixture sigma=%.1f' % s, lambda a, A, rng, S, s=s: model_mix(a, s))
    P('(b2) T1 drifts circuit-to-circuit (every circuit has its own T1), exponential otherwise')
    for s in (0.1, 0.3, 0.6):
        row('drift per circuit sigma=%.1f' % s, lambda a, A, rng, S, s=s: jitter(None, a, s, rng, S))
    P('(c) stretched exponential, hazard keyed to time since circuit start (abs0) / since first X (abs30) / restarted at each X (renew)')
    for clock in ('abs0', 'abs30', 'renew'):
        for b in (0.8, 0.9, 1.0, 1.2, 1.4):
            row('stretch beta=%.1f clock=%s' % (b, clock), lambda a, A, rng, S, b=b, c=clock: model_stretch(a, b, c))
    P('(c2) stretched with each qubit\'s OWN measured beta (3-point fit of ln(-ln R) vs ln t, clipped 0.8-1.4)')
    def bq(A):
        t = np.log([30, 90, 300]); L = np.log(-np.log(np.clip(np.vstack([A['R30'], A['R90'], A['R300']]), .02, .98)))
        return np.clip(np.array([np.polyfit(t, L[:, i], 1)[0] for i in range(L.shape[1])]), 0.8, 1.4)
    for clock in ('abs0', 'abs30', 'renew'):
        row('stretch per-qubit beta clock=%s' % clock, lambda a, A, rng, S, c=clock: model_stretch(a, bq(A), c))
    P('(d) timing: finite X-gate duration (each excited segment is 30+g us, reference waits are 30), exponential and mixture decay')
    for g in (0.035, 0.1, 0.5, 1.0):
        row('exp, gate g=%.3f us' % g, lambda a, A, rng, S, g=g: model_exp(a, gate=g))
    for g in (0.035, 0.1):
        row('mixture sigma=1.2, gate g=%.3f us' % g, lambda a, A, rng, S, g=g: model_mix(a, 1.2, gate=g))
    P()
    main2(P)
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results.text'), 'w').write('\n'.join(out) + '\n')

if __name__ == '__main__':
    main()
