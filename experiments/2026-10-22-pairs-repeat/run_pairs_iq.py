"""Pair run with raw readout: closes the two open causes of the entanglement cost (a leaked neighbour that reads as 0; extra elapsed time).
PREREG_PAIRS.md, sha256 in PREREG_PAIRS.sha256. Dry-run by default; --fire submits. Darkwave Disco instance. Uses run.py's circuit builder
unchanged. One job per copy per processor (two jobs per processor) to keep each raw result a manageable size."""
import os, sys, json, time, random, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_ibm_runtime import SamplerV2
import run as base
HERE = os.path.dirname(os.path.abspath(__file__)); CHIP = base.CHIP; SEED = 20261023; SHOTS = {'CAL': 2000, 'P': 4000}
def build(b):
    r = json.load(open(base.ROLES % CHIP)); n = b.num_qubits; E = sorted(r['emitters']); S = [q for q in range(n) if q not in set(E)]
    P = base.pairs_for(b, S); paired = {q for pr in P for q in pr}; SG = [q for q in S if q not in paired]
    base.L_DT = int(round(b.target['measure'][(E[0],)].duration / b.target.dt))
    K1 = [sorted(random.Random(20261021).sample(E, 10))[0]]          # the same single qubit as the control run's dose arm
    pubs, meta = [], []
    for state in ('0', '1'):
        qc = QuantumCircuit(n, n)
        if state == '1': qc.x(range(n))
        qc.barrier(); qc.measure(range(n), range(n)); pubs.append(qc); meta.append(dict(arm='CAL', state=state))
    for st in base.SETTINGS: pubs.append(base.spectator_circuit(n, E, P, SG, E, 'D', 0, st, '0', False)); meta.append(dict(arm='P', mode='D', N=0, K=39, setting=st))
    for N in (2, 4):
        for mode in ('D', 'M'):
            for st in base.SETTINGS: pubs.append(base.spectator_circuit(n, E, P, SG, E, mode, N, st, '0', False)); meta.append(dict(arm='P', mode=mode, N=N, K=39, setting=st))
    for st in base.SETTINGS: pubs.append(base.spectator_circuit(n, E, P, SG, K1, 'M', 4, st, '0', False)); meta.append(dict(arm='P', mode='M', N=4, K=1, setting=st))
    return pubs, meta, E, S, P, SG, K1
if __name__ == '__main__':
    fire = '--fire' in sys.argv; s = base.service(); rem = s.usage().get('usage_remaining_seconds') or 0; print('remaining seconds:', rem)
    b = s.backend(CHIP); pubs, meta, E, S, P, SG, K1 = build(b)
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); shots = np.array([SHOTS[m['arm']] for m in meta])
    est1 = float(np.sum(shots * (dur + 290e-6)))
    print('%s: %d pubs per copy, pairs %d, singles %d, K1 %s, L=%d dt, total shots per copy %d, ESTIMATED %.0f s per copy, %.0f s per processor (two copies)'
          % (CHIP, len(tq), len(P), len(SG), K1, base.L_DT, shots.sum(), est1, 2 * est1))
    if not fire: sys.exit()
    props = b.properties(); json.dump(props.to_dict() if props is not None else None, open(HERE + '/pairs_backend_properties_%s.json' % CHIP, 'w'), default=str)
    for copy in (0, 1):
        rem = s.usage().get('usage_remaining_seconds') or 0
        if rem < est1 + 3: print('STOP before copy %d: estimated %.0f s, %.0f s left' % (copy, est1, rem)); sys.exit()
        order = list(range(len(tq))); random.Random(SEED + copy).shuffle(order)
        smp = SamplerV2(mode=b); smp.options.execution.meas_type = 'kerneled'
        job = smp.run([(tq[i], None, int(shots[i])) for i in order])
        rec = {'chip': CHIP, 'copy': copy, 'job_id': job.job_id(), 'submitted': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'meta': [dict(meta[i], copy=copy) for i in order],
               'E': E, 'S': S, 'pairs': P, 'singles': SG, 'K1': K1, 'L_dt': base.L_DT, 'meas_type': 'kerneled', 'prereg_sha256': open(HERE + '/PREREG_PAIRS.sha256').read().split()[0]}
        json.dump(rec, open(HERE + '/pairsjob_%s_copy%d.json' % (CHIP, copy), 'w')); print('SUBMITTED', CHIP, 'copy', copy, job.job_id(), flush=True)
        res = job.result(); out = {}
        for i, pr in enumerate(res):
            for name in pr.data: out['%s_%d' % (name, i)] = np.asarray(getattr(pr.data, name)).astype(np.complex64)
        np.savez_compressed(HERE + '/pairsraw_%s_copy%d.npz' % (CHIP, copy), **out)
        print('DONE %s copy %d: remaining now %s' % (CHIP, copy, s.usage().get('usage_remaining_seconds')), flush=True)
