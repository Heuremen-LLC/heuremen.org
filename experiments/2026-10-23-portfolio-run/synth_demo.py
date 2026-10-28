"""synth_demo.py -- planted worlds for demo.py's readings (run before any job). Success per placement planted from a world's true error set
through demo.predict, binomial shots (2 copies x SHOTS), then R1 OURS/IBM read exactly as demo.py --analyze reads it.
  ibm_true  : the hardware's true errors are IBM's published figures   -> OURS must not read ABOVE IBM
  ours_true : the true errors are our measured figures                  -> OURS should read ABOVE IBM when the predictions differ
  equal     : every placement has the same true success                 -> CONSISTENT (the floor; 200 draws, count of false ABOVE/BELOW)"""
import sys, json, numpy as np
sys.argv = ['demo.py', 'x']; import demo as D

chip = 'ibm_kingston'; props = json.load(open(f'../sweep_run/props_sweep_{chip}.json'))
b = D.service().backend(chip); n = b.num_qubits   # the live backend object: coupling map and compiler target only, no job
edges = sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}); paths = D.lines(edges, n)
sxc, czc = D.counts_template(b, paths[0]); ibm = D.ibm_errors(props); ours = D.our_errors(chip)
(pi_, p_ibm), _ = D.best(paths, sxc, czc, *ibm); (po_, p_our), _ = D.best(paths, sxc, czc, *ours)
rng = np.random.default_rng(3); N = 2 * D.SHOTS
def r1(t_ibm, t_our):
    x = rng.random(N) < t_our; y = rng.random(N) < t_ibm; bs = [rng.choice(x, N).mean() / rng.choice(y, N).mean() for _ in range(2000)]
    lo, hi = np.percentile(bs, [2.5, 97.5]); return x.mean() / y.mean(), lo, hi, 'ABOVE' if lo > 1 else 'BELOW' if hi < 1 else 'CONSISTENT'
print(f'picks: IBM {p_ibm} OURS {p_our}')
for w, es in (('ibm_true', ibm), ('ours_true', ours)):
    ti, to = D.predict(p_ibm, sxc, czc, *es), D.predict(p_our, sxc, czc, *es)
    print(f'{w}: true success IBM pick {ti:.4f}, OURS pick {to:.4f} -> R1 {r1(ti, to)}')
out = [r1(0.75, 0.75)[3] for _ in range(200)]
print('equal world, 200 draws: ' + str({k: out.count(k) for k in ('ABOVE', 'BELOW', 'CONSISTENT')}))
