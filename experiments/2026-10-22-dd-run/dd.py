"""DD run, 22 Oct 2026 (Heuremen LLC, Wayfinder and Bones). PLAN_DD.md, sha256 in PLAN_DD.sha256.
Question: does a symmetric X-X echo on the untouched qubits during each step remove (a) the readout phase kick and (b) the readout coherence
cost? Like for like: the echo is applied identically in the idle arm (D) and the readout arm (M).
Dry-run by default: builds, transpiles against the live chip, VERIFIES (X count on every spectator, every step exactly L, DD on/off scheduled
durations equal), estimates seconds, submits nothing. --fire submits. Darkwave Disco instance (control_run/run.py service()).
Usage: python dd.py CHIP [--fire]     Every circuit twice per job (copy 0, copy 1), shuffled together."""
import os, sys, json, time, math, random, warnings; warnings.filterwarnings('ignore')
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import SamplerV2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.append(os.path.join(HERE, '..', 'control_run'))
import run as base
SEED = 20261022; SHOTS = 4000; NS = (0, 2, 4, 8, 16); RUN = 'dd'
def timing(b, E, S):
    """L = published measure duration (dt), identical on every measured qubit (asserted); g = X duration, identical on every spectator.
    Echo per step: delay d1, X, delay d2, X, delay d3 with d1 = d3 = (L/4 - g/2) rounded DOWN to the alignment, d2 = L - 2g - 2 d1."""
    t = b.target; dt = t.dt
    Ls = {int(round(t['measure'][(q,)].duration / dt)) for q in E}; assert len(Ls) == 1, ('measure durations differ on E', Ls)
    gs = {int(round(t['x'][(q,)].duration / dt)) for q in S}; assert len(gs) == 1, ('X durations differ on S', gs)
    L, g = Ls.pop(), gs.pop()
    al = math.lcm(int(getattr(t, 'granularity', 1) or 1), int(getattr(t, 'pulse_alignment', 1) or 1))
    d1 = int((L / 4 - g / 2) // al) * al; d2 = L - 2 * g - 2 * d1; d3 = d1
    assert d1 > 0 and d2 > 0 and L % al == 0 and g % al == 0 and d2 % al == 0, (L, g, al, d1, d2)
    assert d1 + g + d2 + g + d3 == L
    return dict(L_dt=L, g_dt=g, align=al, granularity=int(getattr(t, 'granularity', 1) or 1), pulse_alignment=int(getattr(t, 'pulse_alignment', 1) or 1),
                acquire_alignment=int(getattr(t, 'acquire_alignment', 1) or 1), d1=d1, d2=d2, d3=d3, dt=dt)
def dd_circuit(n, E, S, mode, N, st, dd, T):
    """Spectators S on the equator (SX). N steps, each fenced by barriers on ALL qubits: E measured (M, held in 0) or idle exactly L (D);
    S either wait (DD off: no instruction, the scheduler pads) or run the echo filling exactly L (DD on). Then S rotated to X or Y,
    every qubit measured into 'fin' (qubit q -> bit q)."""
    L = T['L_dt']; qc = QuantumCircuit(n)
    regs = [ClassicalRegister(len(E), 'm%d' % k) for k in range(N if mode == 'M' else 0)]
    for r in regs: qc.add_register(r)
    fin = ClassicalRegister(n, 'fin'); qc.add_register(fin)
    qc.sx(S); qc.barrier()
    for k in range(N):
        if mode == 'M': qc.measure(E, regs[k])
        else: qc.delay(L, E, unit='dt')
        if dd == 'on':
            qc.delay(T['d1'], S, unit='dt'); qc.x(S); qc.delay(T['d2'], S, unit='dt'); qc.x(S); qc.delay(T['d3'], S, unit='dt')
        qc.barrier()
    for q in S: base.rot(qc, q, st[0])
    qc.measure(range(n), fin); return qc
def build(b, chip):
    n = b.num_qubits; E = sorted(json.load(open(base.ROLES % chip))['emitters']); S = [q for q in range(n) if q not in set(E)]
    T = timing(b, E, S); base.L_DT = T['L_dt']; pubs, meta = [], []
    for dd in ('off', 'on'):
        for st in ('XX', 'YY'): pubs.append(dd_circuit(n, E, S, 'D', 0, st, dd, T)); meta.append(dict(dd=dd, mode='D', N=0, setting=st))
        for N in NS[1:]:
            for mode in ('D', 'M'):
                for st in ('XX', 'YY'): pubs.append(dd_circuit(n, E, S, mode, N, st, dd, T)); meta.append(dict(dd=dd, mode=mode, N=N, setting=st))
    return pubs, meta, SHOTS, dict(E=E, S=S, **T)
def verify(tq, meta, info, n, target):
    """On the TRANSPILED, SCHEDULED circuits: (1) X count per spectator = 2N with DD on, 0 with DD off; no X on measured qubits;
    (2) every step (full barrier to full barrier) lasts exactly L; (3) scheduled duration identical between DD on and off for the same
    (mode, N, setting). Raises on any failure; returns a printable summary."""
    S = set(info['S']); L = info['L_dt']; durs = {}
    for c, m in zip(tq, meta):
        xs = {}; bars = []
        for inst, t0 in zip(c.data, c.op_start_times):
            qs = [c.find_bit(q).index for q in inst.qubits]
            if inst.operation.name == 'x':
                for q in qs: xs[q] = xs.get(q, 0) + 1
            if inst.operation.name == 'barrier' and len(qs) == n: bars.append(t0)
        want = 2 * m['N'] if m['dd'] == 'on' else 0
        assert all(xs.get(q, 0) == want for q in S), ('X count', m, sorted(set(xs.values())))
        assert not any(q not in S for q in xs), ('X on a measured qubit', m)
        steps = np.diff(bars); assert len(steps) == m['N'] and np.all(steps == L), ('step lengths', m, steps)
        durs[(m['dd'], m['mode'], m['N'], m['setting'], m['copy'])] = (int(bars[-1] - bars[0]), int(round(c.estimate_duration(target, 'dt'))))
    for k, v in durs.items():
        if k[0] == 'on': assert durs[('off',) + k[1:]] == v, ('DD on/off scheduled durations differ (step block, total)', k, v, durs[('off',) + k[1:]])
    return 'VERIFIED: X per spectator = 2N (DD on) / 0 (DD off), none on measured qubits; every step exactly L = %d dt; step block (N x L) and total scheduled duration identical DD on vs off for every (mode, N, setting, copy)' % L
if __name__ == '__main__':
    CHIP = sys.argv[1]; fire = '--fire' in sys.argv
    s = base.service(); rem = s.usage().get('usage_remaining_seconds') or 0; b = s.backend(CHIP); n = b.num_qubits
    pubs, meta, shots, info = build(b, CHIP)
    meta = [dict(m, copy=c) for c in (0, 1) for m in meta]; pubs = pubs + pubs
    order = list(range(len(pubs))); random.Random(SEED + 11).shuffle(order); pubs = [pubs[i] for i in order]; meta = [meta[i] for i in order]
    tq = transpile(pubs, b, optimization_level=1, initial_layout=list(range(n)), scheduling_method='alap', seed_transpiler=7)
    print('%s timing: dt %.3g s, L %d dt, X %d dt, granularity %d, pulse_alignment %d, acquire_alignment %d -> echo d1 %d, X, d2 %d, X, d3 %d dt'
          % (CHIP, info['dt'], info['L_dt'], info['g_dt'], info['granularity'], info['pulse_alignment'], info['acquire_alignment'], info['d1'], info['d2'], info['d3']))
    print(verify(tq, meta, info, n, b.target))
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(shots * np.sum(dur + 290e-6))
    for dd in ('off', 'on'):
        for N in NS:
            ds = sorted({round(d * 1e6, 3) for d, m in zip(dur, meta) if m['dd'] == dd and m['N'] == N})
            print('  DD %-3s N=%2d circuit durations (us): %s' % (dd, N, ds))
    print('%s %s: %d pubs, %d shots each, E %d, S %d, circuits %.1f-%.1f us, total shots %d, ESTIMATED %.0f s, remaining %s, queue %s'
          % (RUN, CHIP, len(tq), shots, len(info['E']), len(info['S']), dur.min() * 1e6, dur.max() * 1e6, shots * len(tq), est, rem, b.status().pending_jobs), flush=True)
    if not fire: sys.exit()
    if rem < est + 3: print('STOP: estimated %.0f s, %.0f s left' % (est, rem)); sys.exit()
    props = b.properties(); json.dump(props.to_dict() if props is not None else None, open(HERE + '/props_%s_%s.json' % (RUN, CHIP), 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=shots)
    rec = dict(run=RUN, chip=CHIP, n_qubits=n, job_id=job.job_id(), submitted=time.strftime('%Y-%m-%dT%H:%M:%S%z'), meta=meta, order=order, shots=shots,
               info=info, plan_sha256=open(HERE + '/PLAN_DD.sha256').read().split()[0])
    json.dump(rec, open(HERE + '/job_%s_%s.json' % (RUN, CHIP), 'w'), default=str); print('SUBMITTED', RUN, CHIP, job.job_id(), flush=True)
    res = job.result(); out = {}
    for i, pr in enumerate(res):
        for name in pr.data: out['%s_%d' % (name, i)] = getattr(pr.data, name).array
    np.savez_compressed(HERE + '/raw_%s_%s.npz' % (RUN, CHIP), **out)
    print('DONE %s %s: remaining now %s' % (RUN, CHIP, s.usage().get('usage_remaining_seconds')), flush=True)
