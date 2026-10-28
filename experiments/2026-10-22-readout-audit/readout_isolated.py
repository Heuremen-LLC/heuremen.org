"""readout_isolated.py -- the qubits IBM states P(1|0) = 0, read ALONE (every other qubit untouched and unmeasured), against the same qubits
in the all-qubit run 1 (PLAN_READOUT_REPEAT.md).
  readout_isolated.py CHIP [--fire] | --synth WORLD (same, crosstalk) | --analyze CHIP
Circuits: P0 = the listed qubits measured from |0>; P1 = X on them, measured. Each twice (copies), shuffled, 10000 shots (20000 per state)."""
import os, sys, json, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import readout_audit as RA
HERE = RA.HERE; SHOTS = 10000; SEED = 20261026

def zeros_of(props):
    return [qi for qi, q in enumerate(props['qubits']) if any(x['name'] == 'prob_meas1_prep0' and x.get('value') == 0 for x in q)]

def compare(iso, sim, n_iso, n_sim, label=''):
    """iso[q], sim[q]: counts of 1 after |0>. Per-qubit difference isolated - simultaneous, median over qubits, bootstrap; plus R3 on isolated."""
    qs = sorted(iso); d = np.array([iso[q] / n_iso - sim[q] / n_sim for q in qs])
    r = np.random.default_rng(2); b = [np.median(r.choice(d, len(d))) for _ in range(4000)]; lo, hi = np.percentile(b, [2.5, 97.5])
    above = [q for q in qs if RA.wilson(iso[q], n_iso)[0] > RA.ZERO_UB]
    print(f'  {label}R1 isolated - simultaneous P(1|0), median over {len(qs)} qubits: {np.median(d):+.4f} [{lo:+.4f}, {hi:+.4f}] -> {RA.outcome(lo, hi, 0)}')
    print(f'  {label}R3 isolated: stated-zero qubits with delivered P(1|0) 95 % lower bound above {RA.ZERO_UB:.5f}: {len(above)} of {len(qs)}; delivered median {np.median([iso[q] / n_iso for q in qs]):.4f}')

if __name__ == '__main__':
    if '--synth' in sys.argv:
        w = sys.argv[sys.argv.index('--synth') + 1]; rng = np.random.default_rng(6); n = 2 * SHOTS
        true = {q: rng.uniform(0.002, 0.02) for q in range(14)}
        sim = {q: rng.binomial(n, true[q]) for q in true}; iso = {q: rng.binomial(n, true[q] * (0.5 if w == 'crosstalk' else 1.0)) for q in true}
        print(f'=== synthetic world {w} (same: isolated = simultaneous; crosstalk: simultaneous readout doubles P(1|0))'); compare(iso, sim, n, n); sys.exit(0)
    if '--analyze' in sys.argv:
        CHIP = sys.argv[sys.argv.index('--analyze') + 1]
        rec = json.load(open(f'{HERE}/job_iso_{CHIP}.json')); z = np.load(f'{HERE}/raw_iso_{CHIP}.npz'); Q = rec['qubits']
        iso = {q: 0 for q in Q}; n_iso = 0
        for i, m in enumerate(rec['meta']):
            if m['state'] == 0:
                b = z[f'c_{i}']; n_iso += b.shape[0]
                for j, q in enumerate(Q): iso[q] += int(b[:, j].sum())
        r1 = json.load(open(f'{HERE}/job_ro_{CHIP}_run1.json')); z1 = np.load(f'{HERE}/raw_ro_{CHIP}_run1.npz'); sim = {q: 0 for q in Q}; n_sim = 0
        for i, m in enumerate(r1['meta']):
            if m['state'] == 0:
                b = z1[f'c_{i}']; n_sim += b.shape[0]
                for q in Q: sim[q] += int(b[:, q].sum())
        print(f'=== {CHIP} job {rec["job_id"]}; qubits stated P(1|0)=0 at firing: {Q}; {n_iso} isolated shots vs {n_sim} all-qubit shots (run 1)')
        for q in Q: print(f'    q{q}: isolated {iso[q] / n_iso:.4f}  all-qubit {sim[q] / n_sim:.4f}')
        k1 = sum(int((~z[f'c_{i}']).sum()) for i, m in enumerate(rec['meta']) if m['state'] == 1); n1 = sum(z[f'c_{i}'].shape[0] for i, m in enumerate(rec['meta']) if m['state'] == 1) * len(Q)
        kg = 1 - k1 / n1; print(f'  KNOWN-GOOD mean P(1 | |1>) over listed qubits {kg:.4f} (must be >= 0.80): {"MET" if kg >= 0.80 else "NOT MET"}')
        compare(iso, sim, n_iso, n_sim); sys.exit(0)
    CHIP = sys.argv[1]; FIRE = '--fire' in sys.argv
    from qiskit import QuantumCircuit, transpile
    s = RA.service(); b = s.backend(CHIP); props = b.properties().to_dict(); Q = zeros_of(json.loads(json.dumps(props, default=str)))
    k = len(Q); c0 = QuantumCircuit(b.num_qubits, k); c0.measure(Q, range(k)); c1 = QuantumCircuit(b.num_qubits, k); c1.x(Q); c1.measure(Q, range(k))
    items = [(c, dict(state=st, copy=cp)) for cp in (0, 1) for st, c in ((0, c0), (1, c1))]; random.Random(SEED).shuffle(items)
    tq = transpile([c for c, _ in items], b, optimization_level=0, initial_layout=list(range(b.num_qubits)), seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6)); rem = s.usage().get('usage_remaining_seconds')
    print(f'{CHIP}: stated-zero qubits {Q}, {len(tq)} circuits, shots {SHOTS}, ESTIMATED {est:.0f} s, remaining {rem}')
    if not FIRE or not Q: sys.exit(0)
    if rem is not None and rem < est + 3: print('STOP: not enough seconds'); sys.exit(1)
    from qiskit_ibm_runtime import SamplerV2
    json.dump(props, open(f'{HERE}/props_iso_{CHIP}.json', 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    json.dump(dict(chip=CHIP, job_id=job.job_id(), submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), qubits=Q, meta=[m for _, m in items], shots=SHOTS),
              open(f'{HERE}/job_iso_{CHIP}.json', 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
