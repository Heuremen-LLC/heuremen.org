"""Edges IBM states as zz = 0 on ibm_phoenix (record in force at submission): per-copy direct shift with the committed analysis's own
functions, and the two-qubit gates IBM's same record carries on those edges. Defined after the data (the plan's D5 lists these edges)."""
import sys, json
sys.argv = ['x', 'ibm_phoenix']
import analyze_direct as A
rec, z = A.load('ibm_phoenix'); props = json.load(open('props_direct_ibm_phoenix.json'))
zero = {tuple(sorted(int(x) for x in (k, ))) for k in []}
for e in ((100, 101), (70, 80), (86, 96)):
    gates = [(g['gate'], [x['value'] for x in g['parameters'] if x['name'] == 'gate_error'][0]) for g in props['gates'] if g['qubits'] == list(e)]
    vals = []
    for cp in ((0,), (1,), (0, 1)):
        sh = A.shifts(rec, z, cp)
        for j, (b1, a, b2) in enumerate(rec['triples']):
            for arm, b in (('O1', b1), ('O2', b2)):
                if sorted((a, b)) == list(e): vals.append(sh[(j, arm)] / 1e3)
    print(f'edge {e[0]}-{e[1]}: stated zz 0.00 kHz; shift copy 0 {vals[0]:.2f}, copy 1 {vals[1]:.2f}, pooled {vals[2]:.2f} kHz; '
          f'IBM record gates on this edge: ' + ', '.join(f'{g} error {v:.4f}' for g, v in gates))
