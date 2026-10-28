"""Synthetic test of analyze_freshclock.py on FakeFez. Three worlds: NO FADE (kick 2.0 deg per readout whatever came before), FADE (2.0 after
none, 1.4 after 27, 0.8 after 72 readouts, 2.0 after 72 idles), and NULL (no kick at all). Fresh clocks, coherence 0.85."""
import os, sys, json, warnings; warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import freshclock, analyze as A, analyze_freshclock as AF
from qiskit_ibm_runtime.fake_provider import FakeFez
b = FakeFez(); n = b.num_qubits; pubs, meta0, shots, info = freshclock.build_freshclock(b, 'ibm_fez')
for k, (label, size) in enumerate((('NO FADE', {'M0': 2.0, 'M27': 2.0, 'M72': 2.0, 'D72': 2.0}), ('FADE', {'M0': 2.0, 'M27': 1.4, 'M72': 0.8, 'D72': 2.0}), ('NULL', {'M0': 0., 'M27': 0., 'M72': 0., 'D72': 0.}))):
    rng = np.random.default_rng(100 + k); OUT = os.path.join(HERE, 'synth_fresh_%d' % k); os.makedirs(OUT, exist_ok=True)
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta0]; z = {}; nS = len(info['S'])
    for i, m in enumerate(meta):
        phi = np.radians(-size[m['pre']] * 9 if m['win'] == 'M' else 0.0) * np.ones(nS); ex = 0.85 * (np.cos(phi) if m['setting'] == 'XX' else np.sin(phi))
        z['fin_%d' % i] = (rng.random((shots, nS)) < (1 - ex) / 2).astype(np.int8)
    json.dump(b.properties().to_dict(), open(os.path.join(OUT, 'props_freshclock_ibm_fez.json'), 'w'), default=str)
    json.dump(dict(run='freshclock', chip='ibm_fez', n_qubits=n, job_id='SYNTH-' + label, meta=meta, info=info, shots=shots), open(os.path.join(OUT, 'job_freshclock_ibm_fez.json'), 'w'), default=str)
    np.savez_compressed(os.path.join(OUT, 'raw_freshclock_ibm_fez.npz'), **z)
    print('#', label, size); A.LOG.clear(); J = A.Job('freshclock', 'ibm_fez', OUT); J.props_n = n; AF.freshclock(J)
