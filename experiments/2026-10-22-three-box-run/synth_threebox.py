"""synth_threebox.py -- runs the committed analysis on planted data before any job. Three planted worlds, each on three synthetic chips:
  NULL     hardware = IBM model exactly (sampled shots): R1 and R2 must read CONSISTENT
  SIGN     negative box's probe signal scaled by 0.55, positive boxes by 0.9: R1 BELOW-in-size (hardware P less negative => ABOVE), R2 BELOW in every variant
  CORNER   whichever box sits at |11> scaled by 0.55, others 0.9: R2 BELOW for V11, ABOVE for VA (its |11> box is positive), CONSISTENT-or-ABOVE for V10
The 'IBM model' here is the noiseless ideal times 0.95 (a stand-in); shots 6000 per copy (the smaller of the two planned chips)."""
import json, numpy as np, sys
from qiskit.quantum_info import Statevector
from tb_common import EPS, BOXES, VARIANTS, circuit, keys
import analyze_threebox as A
rng = np.random.default_rng(7); SHOTS = 6000
def probs(kk):
    p = Statevector(circuit(*kk).remove_final_measurements(inplace=False)).probabilities([0, 1, 2]); return p
def world(name, scale):
    raws, models = [], []
    for chip in ('synthA', 'synthB', 'synthC'):
        raw, mod, ide = {}, {}, {}
        for kk in keys():
            v, b, e = kk; p = probs(kk); keep = p[0] + p[4]; xi = (p[0] - p[4]) / keep; key = f'{v}|{b}|{e:.2f}'
            ide[key] = float(xi); mod[key] = float(0.95 * xi); xh = 0.95 * xi * scale(v, b)
            for cp in (0, 1):
                n = SHOTS; kept = rng.binomial(n, keep); up = rng.binomial(kept, (1 + xh) / 2)
                bits = np.concatenate([np.zeros(up, np.uint8), np.full(kept - up, 4, np.uint8), np.full(n - kept, 1, np.uint8)])
                raw[f'{key}|{cp}'] = rng.permutation(bits)
        np.savez_compressed(f'/tmp/_tb_{chip}.npz', **raw); json.dump({'chip': chip, 'model': mod, 'ideal': ide}, open(f'/tmp/_tb_{chip}.json', 'w'))
        raws.append(f'/tmp/_tb_{chip}.npz'); models.append(json.load(open(f'/tmp/_tb_{chip}.json')))
    res = [A.analyze(r, m) for r, m in zip(raws, models)]
    return f'\n######## WORLD {name} ########' + A.report(res) + A.combined(res)
neg = lambda v, b: VARIANTS[v]['sign'].get(b, 1) < 0
print(world('NULL', lambda v, b: 1.0))
print(world('SIGN', lambda v, b: 0.55 if neg(v, b) else 0.9))
print(world('CORNER', lambda v, b: 0.55 if b == '11' else 0.9))
