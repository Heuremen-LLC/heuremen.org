"""Synthetic test of analyze_excitation.py on FakeFez: planted G (kick in 1 minus kick in 0) = -0.6 deg/step, and a null (G = 0)."""
import os, sys, json, warnings; warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import excitation, runs, analyze as A, analyze_excitation as AE
from qiskit_ibm_runtime.fake_provider import FakeFez
b = FakeFez(); n = b.num_qubits
for label, G in (('planted', -0.6), ('null', 0.0)):
    rng = np.random.default_rng(9); OUT = os.path.join(HERE, 'synth_exc_' + label); os.makedirs(OUT, exist_ok=True)
    pubs, meta, shots, info = excitation.build_excitation(b, 'ibm_fez'); meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; z = {}
    nS = len(info['S'])
    for i, m in enumerate(meta):
        e = m['mode'][1]; per = (2.0 if m['mode'][0] == 'M' else 0.0) + (3.0 if e == '1' else 0.0) + (G if (m['mode'] == 'M1') else 0.0)
        phi = np.radians(per * m['N']) * np.ones(nS); ex = 0.8 * (np.cos(phi) if m['setting'] == 'XX' else np.sin(phi))
        z['fin_%d' % i] = (rng.random((shots, nS)) < (1 - ex) / 2).astype(np.int8)
    props = b.properties().to_dict(); json.dump(props, open(os.path.join(OUT, 'props_excitation_ibm_fez.json'), 'w'), default=str)
    json.dump(dict(run='excitation', chip='ibm_fez', n_qubits=n, job_id='SYNTH-' + label, meta=meta, info=info, shots=shots), open(os.path.join(OUT, 'job_excitation_ibm_fez.json'), 'w'), default=str)
    np.savez_compressed(os.path.join(OUT, 'raw_excitation_ibm_fez.npz'), **z)
    print('#', label, 'G planted', G); J = A.Job('excitation', 'ibm_fez', OUT); J.props_n = n; AE.excitation(J)
