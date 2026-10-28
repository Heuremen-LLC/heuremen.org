"""readout_audit.py -- delivered readout error against IBM's stated readout fields, same hour (PLAN_READOUT.md).
  readout_audit.py CHIP [--fire]     dry run / submit (open instance of .ibm_quantum_key_3)
  readout_audit.py --synth WORLD     analysis on simulated data (match, worse, zeros_real)
  readout_audit.py --analyze CHIP    read raw_ro_<chip>.npz + props_ro_<chip>.json + job_ro_<chip>.json
Circuits: P0 = every qubit idle in |0>, measure all; P1 = X on every qubit, measure all (IBM: "preparing the qubit in the |0> (|1>) state").
Each twice (copy 0 / copy 1), shuffled, 10000 shots -> 20000 shots per state per qubit."""
import os, sys, json, time, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge')
SHOTS = 10000; SEED = 20261025; KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'open-instance'
ZERO_UB = 3.0 / 4096   # 95 % upper bound on a rate when 0 of 4096 shots err (IBM's stated values are multiples of 1/4096)

def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h

def boot_med(x, seed=1):
    r = np.random.default_rng(seed); x = np.asarray(x); m = [np.median(r.choice(x, len(x))) for _ in range(4000)]
    return float(np.median(x)), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

def outcome(lo, hi, ref):
    return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'

def readings(k10, k01, n, stated, copies_gap=None):
    """k10[q]: counts of 1 after |0>; k01[q]: counts of 0 after |1>; n shots per state; stated[q] = (readout_error, p10, p01)."""
    qs = sorted(q for q in k10 if q in stated)
    d10 = {q: k10[q] / n for q in qs}; d01 = {q: k01[q] / n for q in qs}; dre = {q: (d10[q] + d01[q]) / 2 for q in qs}
    P0 = np.array([dre[q] for q in qs]); kg = float(np.median([1 - d10[q] for q in qs]))
    print(f'  KNOWN-GOOD median P(0 | |0>) = {kg:.4f} (must be >= 0.90): {"MET" if kg >= 0.90 else "NOT MET"}')
    rat = [dre[q] / stated[q][0] for q in qs if stated[q][0] > 0]; m, lo, hi = boot_med(rat)
    print(f'  R1 delivered / stated readout_error, median over {len(rat)} qubits: {m:.3f} [{lo:.3f}, {hi:.3f}] vs 1 -> {outcome(lo, hi, 1)}  (above stated: {sum(r > 1 for r in rat)})')
    for name, dd, j in (('P(1|0)', d10, 1), ('P(0|1)', d01, 2)):
        rr = [dd[q] / stated[q][j] for q in qs if stated[q][j] > 0]; m, lo, hi = boot_med(rr, 2)
        print(f'     {name}: delivered / stated median {m:.3f} [{lo:.3f}, {hi:.3f}] -> {outcome(lo, hi, 1)}')
    bad = [q for q in qs if abs(stated[q][0] - (stated[q][1] + stated[q][2]) / 2) > 0.001]
    if bad:
        near_re = sum(abs(dre[q] - stated[q][0]) < abs(dre[q] - (stated[q][1] + stated[q][2]) / 2) for q in bad)
        print(f'  R2 qubits whose stated readout_error differs from the mean of its two components by > 0.001: {len(bad)}; delivered closer to readout_error on {near_re}, to the mean on {len(bad) - near_re}')
    else:
        print('  R2 no qubit in this record has readout_error differing from the mean of its components by > 0.001')
    z = [q for q in qs if stated[q][1] == 0]
    above = [q for q in z if wilson(k10[q], n)[0] > ZERO_UB]
    print(f'  R3 qubits stated prob_meas1_prep0 = 0: {len(z)}; delivered P(1|0) with 95 % lower bound above {ZERO_UB:.5f} (0 of 4096): {len(above)}'
          + (f'; their delivered P(1|0) median {np.median([k10[q] / n for q in above]):.4f}' if above else ''))
    if copies_gap is not None:
        g, glo, ghi = copies_gap; print(f'  copy gap in readout error (median per qubit, copy 0 - copy 1): {g:+.4f} [{glo:+.4f}, {ghi:+.4f}] -> {"CANNOT_TELL (copies shifted)" if glo > 0 or ghi < 0 else "copies agree"}')

def synth(world, seed=4):
    rng = np.random.default_rng(seed); n = 2 * SHOTS; stated = {}; k10 = {}; k01 = {}
    for q in range(150):
        p10 = 0.0 if q % 10 == 0 else rng.uniform(0.002, 0.02); p01 = rng.uniform(0.005, 0.05)
        re = (p10 + p01) / 2 if q % 15 else re_bad(rng, p10, p01)
        stated[q] = (re, p10, p01)
        f = 1.5 if world == 'worse' else 1.0
        t10 = (0.003 if (p10 == 0 and world == 'zeros_real') else p10 * f); t01 = p01 * f
        k10[q] = rng.binomial(n, t10); k01[q] = rng.binomial(n, t01)
    return k10, k01, n, stated

def re_bad(rng, a, b):
    return (a + b) / 2 + 0.01   # planted: stated readout_error off its own definition on one qubit in 15

def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)

if __name__ == '__main__':
    if '--synth' in sys.argv:
        w = sys.argv[sys.argv.index('--synth') + 1]; print(f'=== synthetic world {w} (match: delivered = stated; worse: delivered 1.5 x stated; zeros_real: stated zeros are really 0.3 %)')
        readings(*synth(w)); sys.exit(0)
    if '--analyze' in sys.argv:
        CHIP = sys.argv[sys.argv.index('--analyze') + 1]
        rec = json.load(open(f'{HERE}/job_ro_{CHIP}.json')); z = np.load(f'{HERE}/raw_ro_{CHIP}.npz'); props = json.load(open(f'{HERE}/props_ro_{CHIP}.json'))
        stated = {}
        for q, ps in enumerate(props['qubits']):
            v = {x['name']: x.get('value') for x in ps}
            if all(v.get(k) is not None for k in ('readout_error', 'prob_meas1_prep0', 'prob_meas0_prep1')): stated[q] = (v['readout_error'], v['prob_meas1_prep0'], v['prob_meas0_prep1'])
        k10 = {}; k01 = {}; per = {0: {}, 1: {}}; n = 0
        for i, m in enumerate(rec['meta']):
            b = z[f'c_{i}']; cnt = b.sum(axis=0) if m['state'] == 0 else (~b).sum(axis=0)
            tgt = k10 if m['state'] == 0 else k01
            for q in range(b.shape[1]):
                tgt[q] = tgt.get(q, 0) + int(cnt[q]); per[m['copy']].setdefault(q, [0, 0])[m['state']] += int(cnt[q])
            n_state = b.shape[0]
        n = 2 * n_state
        g = [(per[0][q][0] + per[0][q][1]) / (2 * n_state) - (per[1][q][0] + per[1][q][1]) / (2 * n_state) for q in per[0] if q in per[1]]
        r = np.random.default_rng(3); gb = [np.median(r.choice(g, len(g))) for _ in range(4000)]
        print(f'=== {CHIP} job {rec["job_id"]}; IBM record {props["last_update_date"]}; {n} shots per state per qubit')
        readings(k10, k01, n, stated, (float(np.median(g)), *np.percentile(gb, [2.5, 97.5]))); sys.exit(0)
    CHIP = sys.argv[1]; FIRE = '--fire' in sys.argv
    from qiskit import QuantumCircuit, transpile
    s = service(); b = s.backend(CHIP); nq = b.num_qubits
    c0 = QuantumCircuit(nq, nq); c0.measure(range(nq), range(nq))
    c1 = QuantumCircuit(nq, nq); c1.x(range(nq)); c1.measure(range(nq), range(nq))
    items = [(c, dict(state=st, copy=cp)) for cp in (0, 1) for st, c in ((0, c0), (1, c1))]; random.Random(SEED).shuffle(items)
    tq = transpile([c for c, _ in items], b, optimization_level=0, initial_layout=list(range(nq)), seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6)); rem = s.usage().get('usage_remaining_seconds')
    print(f'{CHIP}: {nq} qubits, {len(tq)} circuits, shots {SHOTS}, ESTIMATED {est:.0f} s, remaining {rem}')
    if not FIRE: sys.exit(0)
    if rem is not None and rem < est + 3: print('STOP: not enough seconds'); sys.exit(1)
    from qiskit_ibm_runtime import SamplerV2
    json.dump(b.properties().to_dict(), open(f'{HERE}/props_ro_{CHIP}.json', 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = dict(chip=CHIP, job_id=job.job_id(), submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), meta=[m for _, m in items], shots=SHOTS, keyfile=KEYFILE, instance=INSTANCE)
    json.dump(rec, open(f'{HERE}/job_ro_{CHIP}.json', 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
