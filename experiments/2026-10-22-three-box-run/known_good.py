"""known_good.py -- the committed fit and kept-shot reader applied to the 11 Oct ibm_fez three-box jobs (real hardware, already read on
22 Oct before this plan: corner P = -0.530 +- 0.053, A P = +0.817 +- 0.054). PASS if this code reproduces both to within 0.01."""
import glob, json, numpy as np, os
from tb_common import xp_from_bits, fit_SP
D = os.path.expanduser('~/bones-forge/arxiv_one_readout/history_rule/oct11_threebox')
eps_of = {'circuit-62': .4, 'circuit-90': .3, 'circuit-118': .25, 'circuit-146': .2, 'circuit-174': .15, 'circuit-202': .1}
det = json.load(open(D + '/oct11_unrecorded_detail.json')); data = {}
for f in glob.glob(D + '/db5t*.npz'):
    jid = os.path.basename(f)[:-4]; e = eps_of[det[jid]['circuits'][0][0]]; z = np.load(f)
    for i, n in enumerate(['A_00', 'B_01', 'corner_11', 'unused_10']): data.setdefault(n, {}).setdefault(e, []).append(xp_from_bits(z[f'pub{i}_c'][:, 0]))
ok = True
for n, want in [('A_00', 0.817), ('corner_11', -0.530)]:
    E = sorted(data[n]); x = [np.concatenate(data[n][e]) for e in E]; m = [a.mean() for a in x]; se = [a.std(ddof=1) / np.sqrt(len(a)) for a in x]
    S, P, chi = fit_SP(E, m, se); good = abs(P - want) < 0.01; ok &= good
    print(f'{n}: P {P:+.3f} (want {want:+.3f}) S {S:+.3f} chi2/dof {chi:.2f} -> {"PASS" if good else "FAIL"}')
print('KNOWN-GOOD', 'PASS' if ok else 'FAIL')
