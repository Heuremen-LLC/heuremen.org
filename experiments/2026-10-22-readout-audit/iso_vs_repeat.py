"""iso_vs_repeat.py -- ADDENDUM (/q B13), written and pushed BEFORE the isolated or repeat results were opened: on marrakesh, compare the
isolated read with the all-qubit REPEAT fired in the same minute (raw_ro_ibm_marrakesh.npz, run 2), not only with run 1 (an hour earlier).
Same reading as readout_isolated.compare: per-qubit isolated - all-qubit P(1|0), median, bootstrap, four outcomes."""
import json, numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import readout_isolated as RI
H = RI.HERE; C = 'ibm_marrakesh'
rec = json.load(open(f'{H}/job_iso_{C}.json')); z = np.load(f'{H}/raw_iso_{C}.npz'); Q = rec['qubits']
iso = {q: 0 for q in Q}; ni = 0
for i, m in enumerate(rec['meta']):
    if m['state'] == 0:
        b = z[f'c_{i}']; ni += b.shape[0]
        for j, q in enumerate(Q): iso[q] += int(b[:, j].sum())
r2 = json.load(open(f'{H}/job_ro_{C}.json')); z2 = np.load(f'{H}/raw_ro_{C}.npz'); sim = {q: 0 for q in Q}; ns = 0
for i, m in enumerate(r2['meta']):
    if m['state'] == 0:
        b = z2[f'c_{i}']; ns += b.shape[0]
        for q in Q: sim[q] += int(b[:, q].sum())
print(f'=== {C}: isolated ({rec["job_id"]}) vs all-qubit repeat ({r2["job_id"]}), fired in the same minute')
for q in Q: print(f'    q{q}: isolated {iso[q] / ni:.4f}  all-qubit repeat {sim[q] / ns:.4f}')
RI.compare(iso, sim, ni, ns, 'vs repeat ')
