"""X-FLIP CHECK (PLAN_FLIP.md, sha256 in PLAN_FLIP.sha256). Heuremen LLC, Wayfinder and Bones, 22 October 2026.
Same roles and circuits as control_run/run_iq.py (control_run/run.py's spectator_circuit, unchanged, via a thin wrapper), N = 0, 4, 8,
every readout raw IQ (meas_type kerneled). NEW: after the last step, X on every measured qubit (barrier-fenced), then the measured
qubits' final readout (register 'efin') alongside the spectators'. Calibrations: all 0, all 1, and '0X' (read in 0, X, read again).
Dry-run by default (builds, transpiles, verifies the X survives, estimates seconds; nothing submitted). --fire submits. Darkwave Disco."""
import os, sys, json, time, random, warnings; warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); CTRL = os.path.join(os.path.dirname(HERE), 'control_run'); sys.path.insert(0, CTRL)
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import XGate, Barrier
from qiskit_ibm_runtime import SamplerV2
import run as base
CHIP = base.CHIP; SHOTS = int(base.arg('--shots', 3000)); SEED = 20261023; NS = (4, 8); OVERHEAD = 290e-6

def flip_wrap(qc, n, E):
    """Take spectator_circuit's output; right after its LAST barrier (end of the N steps, or end of preparation when N = 0) insert
    X on every measured qubit fenced by a barrier on all qubits; then append the measured qubits' final readout into 'efin'."""
    last = max(i for i, ci in enumerate(qc.data) if ci.operation.name == 'barrier')
    out = qc.copy_empty_like(); efin = ClassicalRegister(len(E), 'efin'); out.add_register(efin)
    for ci in qc.data[:last + 1]: out.append(ci.operation, ci.qubits, ci.clbits)
    out.x(E); out.barrier()
    for ci in qc.data[last + 1:]: out.append(ci.operation, ci.qubits, ci.clbits)
    out.measure(E, efin); return out

def build(b):
    r = json.load(open(base.ROLES % CHIP)); n = b.num_qubits; E = sorted(r['emitters']); S = [q for q in range(n) if q not in set(E)]
    base.L_DT = int(round(b.target['measure'][(E[0],)].duration / b.target.dt))
    pubs, meta = [], []
    for state in ('0', '1'):            # calibration as run_iq.py: every qubit in 0, or in 1, read once
        qc = QuantumCircuit(n, n)
        if state == '1': qc.x(range(n))
        qc.barrier(); qc.measure(range(n), range(n)); pubs.append(qc); meta.append(dict(arm='CAL', state=state))
    # NEW calibration '0X': every qubit in 0, read once (register c0), X on every qubit (fenced), read again (register c1):
    # the post-flip cluster after one readout, taken only from shots whose first read was verified 0
    qc = QuantumCircuit(n); r0 = ClassicalRegister(n, 'c0'); r1 = ClassicalRegister(n, 'c1'); qc.add_register(r0); qc.add_register(r1)
    qc.barrier(); qc.measure(range(n), r0); qc.barrier(); qc.x(range(n)); qc.barrier(); qc.measure(range(n), r1)
    pubs.append(qc); meta.append(dict(arm='CAL', state='0X'))
    for st in ('XX', 'YY'):
        pubs.append(flip_wrap(base.spectator_circuit(n, E, [], S, E, 'D', 0, st, '0', False), n, E)); meta.append(dict(arm='S', mode='D', N=0, setting=st))
    for N in NS:
        for mode in ('D', 'M'):
            for st in ('XX', 'YY'):
                pubs.append(flip_wrap(base.spectator_circuit(n, E, [], S, E, mode, N, st, '0', False), n, E)); meta.append(dict(arm='S', mode=mode, N=N, setting=st))
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED).shuffle(order)
    return [pubs[i] for i in order], [meta[i] for i in order], order, E, S

def verify(tq, meta, E):
    """The flip must survive transpilation: in every S circuit, each measured qubit carries exactly one x, after its last mid readout
    and before its final readout; in CAL 0X every qubit carries one x between its two readouts."""
    bad = []
    for c, m in zip(tq, meta):
        if m['arm'] == 'CAL' and m['state'] != '0X': continue
        Q = E if m['arm'] == 'S' else list(range(c.num_qubits)); seq = {q: [] for q in Q}
        for ci in c.data:
            nm = ci.operation.name
            if nm in ('x', 'measure', 'sx', 'rz', 'cz', 'h', 'sdg'):
                for qb in ci.qubits:
                    q = c.find_bit(qb).index
                    if q in seq: seq[q].append(nm)
        for q, s in seq.items():
            ops = [o for o in s if o != 'rz']
            if ops.count('x') != 1 or ops[-1] != 'measure' or ops[-2] != 'x' or ('sx' in ops):
                bad.append((m, q, ops[-4:])); break
    return bad

if __name__ == '__main__':
    fire = '--fire' in sys.argv; s = base.service(); rem = s.usage().get('usage_remaining_seconds') or 0; print('remaining seconds:', rem)
    b = s.backend(CHIP); pubs, meta, order, E, S = build(b)
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    bad = verify(tq, meta, E); print('flip check in transpiled circuits:', 'OK (one x per measured qubit, right before its final readout)' if not bad else 'FAILED %s' % bad[:3])
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + OVERHEAD))
    nval = sum(SHOTS * sum(cr.size for cr in c.cregs) for c in tq)
    print('%s: %d pubs, E %d, S %d, L=%d dt, circuits %.1f-%.1f us, shots/circuit %d, total shots %d, ESTIMATED %.1f s, raw IQ ~%.0f MB (5.4 B/value as 21 Oct), queue %s'
          % (CHIP, len(tq), len(E), len(S), base.L_DT, dur.min() * 1e6, dur.max() * 1e6, SHOTS, SHOTS * len(tq), est, nval * 5.43 / 1e6, b.status().pending_jobs))
    if not fire: sys.exit()
    if bad: print('STOP: flip did not survive transpilation'); sys.exit()
    if rem < est + 3: print('STOP: estimated %.0f s, %.0f s left' % (est, rem)); sys.exit()
    props = b.properties(); json.dump(props.to_dict() if props is not None else None, open(HERE + '/flip_backend_properties_%s.json' % CHIP, 'w'), default=str)
    smp = SamplerV2(mode=b); smp.options.execution.meas_type = 'kerneled'
    job = smp.run(tq, shots=SHOTS)
    rec = {'chip': CHIP, 'job_id': job.job_id(), 'submitted': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'meta': meta, 'order': order, 'E': E, 'S': S,
           'L_dt': base.L_DT, 'shots': SHOTS, 'meas_type': 'kerneled', 'plan_sha256': open(HERE + '/PLAN_FLIP.sha256').read().split()[0]}
    json.dump(rec, open(HERE + '/flipjob_%s.json' % CHIP, 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res):
        for name in pr.data: out['%s_%d' % (name, i)] = np.asarray(getattr(pr.data, name)).astype(np.complex64)
    np.savez_compressed(HERE + '/flipraw_%s.npz' % CHIP, **out)
    print('DONE %s: remaining now %s' % (CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
