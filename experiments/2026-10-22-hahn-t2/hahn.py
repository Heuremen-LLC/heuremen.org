"""hahn.py -- delivered Hahn-echo T2 against IBM's stated T2 on ibm_fez, same hour (PLAN_HAHN.md).
  hahn.py ibm_fez                dry run: build, transpile, estimate seconds
  hahn.py ibm_fez --fire         submit one job (open instance of .ibm_quantum_key_3)
  hahn.py --synth WORLD          analysis on simulated data (WORLD = match, short, mixed)
  hahn.py --analyze ibm_fez      read raw_hahn_ibm_fez.npz + props_hahn_ibm_fez.json and print the readings
Two sparse sets of qubits (no two qubits of a set coupled), every other qubit idle in |0>. Per measured qubit: SX, wait tau/2, X, wait tau/2,
SX-dagger, measure; P(1) decays from about 1 toward 1/2 with the echo T2. tau = 0, 10, 20, 40, 60, 90, 120, 160, 220, 300, 400 us. Each circuit twice
(copy 0 / copy 1), shuffled, 2000 shots."""
import os, sys, json, time, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge')
TAUS = [0, 10, 20, 40, 60, 90, 120, 160, 220, 300, 400]; SHOTS = 2000; SEED = 20261024
KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'open-instance'


def sparse_sets(n, edges):
    """Greedy two sets with no edge inside a set: set A by index order, set B from the rest."""
    adj = {q: set() for q in range(n)}
    for a, b in edges: adj[a].add(b); adj[b].add(a)
    sets = []
    rest = list(range(n))
    for _ in range(2):
        s = []
        for q in rest:
            if all(x not in adj[q] for x in s): s.append(q)
        sets.append(s); rest = [q for q in rest if q not in s]
    return sets


def fit_t2(t_us, p1):
    """P(1) = 1/2 + A exp(-t/T2), A and T2 by least squares on a T2 grid (1 to 2000 us). Returns T2, A, rms."""
    t = np.asarray(t_us, float); y = np.asarray(p1, float) - 0.5; best = (np.nan, np.nan, np.inf)
    for T2 in np.exp(np.linspace(np.log(1), np.log(2000), 1200)):
        e = np.exp(-t / T2); A = float(e @ y / (e @ e)); r = float(np.sqrt(np.mean((y - A * e) ** 2)))
        if r < best[2]: best = (float(T2), A, r)
    return best


def readings(p1, stated, qubits):
    """p1: {copy: {q: [P1 per tau]}}; stated: {q: IBM T2 us}. Ratio delivered/stated per qubit, median over qubits, bootstrap 4000."""
    rng = np.random.default_rng(1)
    def ratios(copies):
        out = {}
        for q in qubits:
            y = np.mean([p1[c][q] for c in copies], axis=0); T2, A, r = fit_t2(TAUS, y)
            if A > 0.2 and np.isfinite(T2) and stated.get(q): out[q] = T2 / stated[q]
        return out
    pooled = ratios((0, 1)); c0 = ratios((0,)); c1 = ratios((1,))
    qs = sorted(pooled); x = np.array([pooled[q] for q in qs]); med = float(np.median(x))
    bs = [np.median(rng.choice(x, len(x))) for _ in range(4000)]; lo, hi = np.percentile(bs, [2.5, 97.5])
    common = [q for q in qs if q in c0 and q in c1]; g = np.array([c0[q] - c1[q] for q in common])
    gb = [np.median(rng.choice(g, len(g))) for _ in range(4000)]; glo, ghi = np.percentile(gb, [2.5, 97.5])
    shifted = glo > 0 or ghi < 0
    out = 'CANNOT_TELL' if shifted else ('ABOVE' if lo > 1 else 'BELOW' if hi < 1 else 'CONSISTENT')
    print(f'  delivered / stated T2, median over {len(qs)} qubits: {med:.3f} [{lo:.3f}, {hi:.3f}] vs 1 -> {out}'
          f'   (copy gap {np.median(g):+.3f} [{glo:+.3f}, {ghi:+.3f}]; qubits below stated {int((x < 1).sum())} of {len(x)})')
    return med, lo, hi, out


def synth(world, seed=3):
    rng = np.random.default_rng(seed); qubits = list(range(60)); stated = {q: float(rng.uniform(80, 300)) for q in qubits}
    true = {q: stated[q] * (1.0 if world == 'match' else 0.7 if world == 'short' else (0.7 if q % 2 else 1.0)) for q in qubits}
    p1 = {c: {q: [float(np.clip(0.5 + 0.48 * np.exp(-t / true[q]) + rng.normal(0, np.sqrt(0.25 / SHOTS)), 0, 1)) for t in TAUS] for q in qubits} for c in (0, 1)}
    return p1, stated, qubits


def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)


def build(b):
    from qiskit import QuantumCircuit
    n = b.num_qubits; edges = sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}); sets = sparse_sets(n, edges)
    circs, meta = [], []
    for si, S in enumerate(sets):
        for tau in TAUS:
            qc = QuantumCircuit(n, len(S)); qc.sx(S); qc.barrier()
            if tau: qc.delay(tau / 2, S, unit='us')
            qc.x(S)
            if tau: qc.delay(tau / 2, S, unit='us')
            qc.barrier(); qc.sxdg(S); qc.measure(S, list(range(len(S))))
            circs.append(qc); meta.append(dict(set=si, tau_us=tau))
    return circs, meta, sets, edges


if __name__ == '__main__':
    if '--synth' in sys.argv:
        w = sys.argv[sys.argv.index('--synth') + 1]; print(f'=== synthetic world {w} (match: true = stated; short: true = 0.7 x stated; mixed: half 0.7)')
        p1, st, qs = synth(w); readings(p1, st, qs); sys.exit(0)
    if '--analyze' in sys.argv:
        CHIP = sys.argv[sys.argv.index('--analyze') + 1]
        rec = json.load(open(f'{HERE}/job_hahn_{CHIP}.json')); z = np.load(f'{HERE}/raw_hahn_{CHIP}.npz'); props = json.load(open(f'{HERE}/props_hahn_{CHIP}.json'))
        stated = {}
        for q, ps in enumerate(props['qubits']):
            v = [x['value'] for x in ps if x['name'] == 'T2']
            if v: stated[q] = v[0]
        p1 = {0: {}, 1: {}}
        for i, m in enumerate(rec['meta']):
            S = rec['sets'][m['set']]; bits = z[f'c_{i}'][:, :len(S)]
            for j, q in enumerate(S): p1[m['copy']].setdefault(q, {})[m['tau_us']] = float(bits[:, j].mean())
        p1 = {c: {q: [p1[c][q][t] for t in TAUS] for q in p1[c]} for c in p1}
        allq = sorted(p1[0]); t0 = np.median([p1[0][q][0] for q in allq])
        print(f'=== {CHIP} job {rec["job_id"]}; IBM record {props["last_update_date"]}; KNOWN-GOOD median P(1) at tau 0 = {t0:.3f} (must be >= 0.90): {"MET" if t0 >= 0.90 else "NOT MET"}')
        readings(p1, stated, allq); sys.exit(0)
    CHIP = sys.argv[1]; FIRE = '--fire' in sys.argv
    from qiskit import transpile
    s = service(); b = s.backend(CHIP); circs, meta, sets, edges = build(b)
    items = [(c, dict(m, copy=cp)) for cp in (0, 1) for c, m in zip(circs, meta)]; random.Random(SEED).shuffle(items)
    tq = transpile([c for c, _ in items], b, optimization_level=1, initial_layout=list(range(b.num_qubits)), scheduling_method='alap', seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6)); rem = s.usage().get('usage_remaining_seconds')
    print(f'{CHIP}: sets of {len(sets[0])} and {len(sets[1])} qubits, {len(tq)} circuits, shots {SHOTS}, ESTIMATED {est:.0f} s, remaining {rem}')
    if not FIRE: sys.exit(0)
    if rem is not None and rem < est + 3: print('STOP: not enough seconds'); sys.exit(1)
    from qiskit_ibm_runtime import SamplerV2
    json.dump(b.properties().to_dict(), open(f'{HERE}/props_hahn_{CHIP}.json', 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = dict(chip=CHIP, job_id=job.job_id(), submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), meta=[m for _, m in items], sets=sets, shots=SHOTS, keyfile=KEYFILE, instance=INSTANCE)
    json.dump(rec, open(f'{HERE}/job_hahn_{CHIP}.json', 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
