"""Raw-readout (IQ) run: closes the last ordinary cause of the superposition (T2) cost, a measured neighbour leaked to a level that
reads as 0. PREREG_IQ.md, sha256 in PREREG_IQ.sha256. Dry-run by default; --fire submits. Darkwave Disco instance.
Same roles and circuits as the control run's arm S (run.py), plus two calibration circuits; every readout returned as a raw IQ point."""
import os, sys, json, time, random, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_ibm_runtime import SamplerV2
import run as base
HERE = os.path.dirname(os.path.abspath(__file__)); CHIP = base.CHIP; SHOTS = 2000; SEED = 20261022
def build(b):
    r = json.load(open(base.ROLES % CHIP)); n = b.num_qubits; E = sorted(r['emitters']); S = [q for q in range(n) if q not in set(E)]
    base.L_DT = int(round(b.target['measure'][(E[0],)].duration / b.target.dt))
    pubs, meta = [], []
    for state in ('0', '1'):            # calibration: every qubit prepared in 0, or in 1, and read once
        qc = QuantumCircuit(n, n)
        if state == '1': qc.x(range(n))
        qc.barrier(); qc.measure(range(n), range(n)); pubs.append(qc); meta.append(dict(arm='CAL', state=state))
    for st in ('XX', 'YY'): pubs.append(base.spectator_circuit(n, E, [], S, E, 'D', 0, st, '0', False)); meta.append(dict(arm='S', mode='D', N=0, setting=st))
    for N in (2, 4, 8):
        for mode in ('D', 'M'):
            for st in ('XX', 'YY'):
                pubs.append(base.spectator_circuit(n, E, [], S, E, mode, N, st, '0', False)); meta.append(dict(arm='S', mode=mode, N=N, setting=st))
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED).shuffle(order)
    return [pubs[i] for i in order], [meta[i] for i in order], order, E, S
if __name__ == '__main__':
    fire = '--fire' in sys.argv; s = base.service(); rem = s.usage().get('usage_remaining_seconds') or 0; print('remaining seconds:', rem)
    b = s.backend(CHIP); pubs, meta, order, E, S = build(b)
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6))
    print('%s: %d pubs, E %d, S %d, L=%d dt, circuits %.1f-%.1f us, total shots %d, ESTIMATED %.0f s, queue %s'
          % (CHIP, len(tq), len(E), len(S), base.L_DT, dur.min() * 1e6, dur.max() * 1e6, SHOTS * len(tq), est, b.status().pending_jobs))
    if not fire: sys.exit()
    if rem < est + 3: print('STOP: estimated %.0f s, %.0f s left' % (est, rem)); sys.exit()
    props = b.properties(); json.dump(props.to_dict() if props is not None else None, open(HERE + '/iq_backend_properties_%s.json' % CHIP, 'w'), default=str)
    smp = SamplerV2(mode=b); smp.options.execution.meas_type = 'kerneled'
    job = smp.run(tq, shots=SHOTS)
    rec = {'chip': CHIP, 'job_id': job.job_id(), 'submitted': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'meta': meta, 'order': order, 'E': E, 'S': S,
           'L_dt': base.L_DT, 'shots': SHOTS, 'meas_type': 'kerneled', 'prereg_sha256': open(HERE + '/PREREG_IQ.sha256').read().split()[0]}
    json.dump(rec, open(HERE + '/iqjob_%s.json' % CHIP, 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res):
        for name in pr.data: out['%s_%d' % (name, i)] = np.asarray(getattr(pr.data, name)).astype(np.complex64)
    np.savez_compressed(HERE + '/iqraw_%s.npz' % CHIP, **out)
    print('DONE %s: remaining now %s' % (CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
