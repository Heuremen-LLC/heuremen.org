"""synth_addendum.py -- planted worlds for the four-state addendum, two synthetic chips each (as the real run): NULL, SIGN, CORNER, BOTH, POS00."""
import json, numpy as np
from qiskit.quantum_info import Statevector
from tb_common import VARIANTS, circuit, keys
import analyze_threebox as A, addendum_fourstate as F
rng = np.random.default_rng(11); SHOTS = 6000
P = {kk: Statevector(circuit(*kk).remove_final_measurements(inplace=False)).probabilities([0, 1, 2]) for kk in keys()}
def world(name, scale):
    out = [f'######## WORLD {name}']
    for chip in ('synthA', 'synthB'):
        raw, mod, ide = {}, {}, {}
        for kk in keys():
            v, b, e = kk; p = P[kk]; keep = p[0] + p[4]; xi = (p[0] - p[4]) / keep; key = f'{v}|{b}|{e:.2f}'
            ide[key] = float(xi); mod[key] = float(0.95 * xi); xh = 0.95 * xi * scale(v, b)
            for cp in (0, 1):
                kept = rng.binomial(SHOTS, keep); up = rng.binomial(kept, (1 + xh) / 2)
                raw[f'{key}|{cp}'] = rng.permutation(np.concatenate([np.zeros(up, np.uint8), np.full(kept - up, 4, np.uint8), np.full(SHOTS - kept, 1, np.uint8)]))
        np.savez_compressed(f'/tmp/_tba_{chip}.npz', **raw)
        out.append(F.fourstate(A.analyze(f'/tmp/_tba_{chip}.npz', {'chip': chip, 'model': mod, 'ideal': ide}), f'/tmp/_tba_{chip}.npz'))
    return '\n'.join(out)
neg = lambda v, b: VARIANTS[v]['sign'].get(b, 1) < 0
print(world('NULL (expect NEITHER)', lambda v, b: 1.0))
print(world('SIGN x0.6 (expect SIGN)', lambda v, b: 0.6 if neg(v, b) else 1.0))
print(world('CORNER x0.6 (expect CORNER)', lambda v, b: 0.6 if b == '11' else 1.0))
print(world('BOTH x0.75 each (expect BOTH)', lambda v, b: (0.75 if neg(v, b) else 1.0) * (0.75 if b == '11' else 1.0)))
print(world('|00> POSITION x0.6 (expect NEITHER on a and b, c BELOW)', lambda v, b: 0.6 if b == '00' else 1.0))
