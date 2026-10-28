"""probe_split.py -- where does the three-box probe offset come from? (PLAN_PROBESPLIT.md)
  probe_split.py CHIP                 dry run: build, transpile, estimate seconds
  probe_split.py CHIP --fire          submit one job (open instance of .ibm_quantum_key_3)
  probe_split.py --synth WORLD        run the analysis on a simulated job (WORLD = null, readout, system, both)
  probe_split.py --analyze CHIP       read raw_split_<chip>.npz and print the readings
Circuits (qubits s0, s1, probe in the 22 Oct / 11 Oct layout), each twice (copy 0 / copy 1), shuffled, 10000 shots:
  K   probe |+>, read in X                      known-good: <X> must be within 0.03 of 1 - 2 * IBM's stated prob_meas1_prep0 of the probe
  P0  probe |0>, read in X, system idle          the probe's own X-readout bias
  P2_<var>  system pre-selected psi(var), post-rotation phi(var) -> |00>, NO controlled gate, probe |0> read in X; kept s0 = s1 = 0
The three-box offset S (fit intercept of the unused box, tb_common.fit_SP) is read from the three-box run of the same chip."""
import os, sys, json, time, random
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tb_common as T
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge')
LAYOUT = {'ibm_marrakesh': [6, 7, 8], 'ibm_fez': [142, 143, 144]}
SHOTS = 10000; SEED = 20261023; KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'open-instance'
UNUSED = {'V11': '10', 'V10': '11', 'VA': '10'}


def circuits():
    out = []
    k = QuantumCircuit(3, 3, name='K'); k.h(2); k.h(2); k.measure([0, 1, 2], [0, 1, 2]); out.append(('K', k))      # |+>, then X read
    p0 = QuantumCircuit(3, 3, name='P0'); p0.h(2); p0.measure([0, 1, 2], [0, 1, 2]); out.append(('P0', p0))
    for var, V in T.VARIANTS.items():
        qc = QuantumCircuit(3, 3, name='P2_' + var)
        qc.append(StatePreparation(T.vec(V['support'])), [0, 1])
        qc.append(StatePreparation(T.vec(V['support'], V['sign'])).inverse(), [0, 1])
        qc.h(2); qc.measure([0, 1, 2], [0, 1, 2]); out.append(('P2_' + var, qc))
    return out


def xp(bits, kept):
    b = np.asarray(bits)
    if kept:
        return T.xp_from_bits(b)
    return 1.0 - 2.0 * ((b >> 2) & 1)


def boot(x, B=4000, seed=1):
    r = np.random.default_rng(seed); x = np.asarray(x); m = [r.choice(x, len(x)).mean() for _ in range(B)]
    return float(x.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def outcome(lo, hi, ref=0.0):
    return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'


def analyze(data, readout_error, offsets):
    """data: {name: {copy: bits array}}; readout_error: IBM's stated probe readout_error; offsets: {var: (S, lo, hi)} from the three-box fit."""
    res = {}
    for name in data:
        kept = name.startswith('P2')
        c = {cp: xp(data[name][cp], kept) for cp in data[name]}
        pooled = np.concatenate([c[0], c[1]]); m, lo, hi = boot(pooled)
        g = c[0].mean() - c[1].mean(); gw = 1.96 * np.sqrt(c[0].var() / len(c[0]) + c[1].var() / len(c[1]))
        res[name] = dict(mean=m, lo=lo, hi=hi, n=len(pooled), gap=g, gap_lo=g - gw, gap_hi=g + gw)
    k = res['K']; target = 1 - 2 * readout_error; kg = abs(k['mean'] - target) < 0.03
    print(f"KNOWN-GOOD  K <X> {k['mean']:+.4f} [{k['lo']:+.4f}, {k['hi']:+.4f}] against 1 - 2 x stated prob_meas1_prep0 {target:+.4f}: {'MET' if kg else 'NOT MET'}")
    for name, r in res.items():
        gap = 'copies shifted -> CANNOT_TELL' if (r['gap_lo'] > 0 or r['gap_hi'] < 0) else 'copies agree'
        print(f"  {name:<7} <X> {r['mean']:+.4f} [{r['lo']:+.4f}, {r['hi']:+.4f}] n={r['n']}  vs 0: {outcome(r['lo'], r['hi'])}  ({gap}, gap {r['gap']:+.4f} [{r['gap_lo']:+.4f}, {r['gap_hi']:+.4f}])")
    p0 = res['P0']
    for var, (S, slo, shi) in offsets.items():
        p2 = res['P2_' + var]
        d1 = p2['mean'] - p0['mean']; w1 = np.hypot((p2['hi'] - p2['lo']) / 2, (p0['hi'] - p0['lo']) / 2)
        d2 = S - p2['mean']; w2 = np.hypot((shi - slo) / 2, (p2['hi'] - p2['lo']) / 2)
        print(f"  {var} unused box {UNUSED[var]}: offset S {S:+.4f} = probe readout P0 {p0['mean']:+.4f} + system ops (P2 - P0) {d1:+.4f} [{d1 - w1:+.4f}, {d1 + w1:+.4f}] {outcome(d1 - w1, d1 + w1)}"
              f" + coupling gate (S - P2) {d2:+.4f} [{d2 - w2:+.4f}, {d2 + w2:+.4f}] {outcome(d2 - w2, d2 + w2)}")
    return res, kg


def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)


def synth(world, seed=7):
    """Simulated job: ideal circuit statistics plus planted effects. null: stated-like readout 1 % symmetric; readout: probe X-readout bias
    +0.06; system: system ops shift the probe by +0.06 on kept shots; both: both. Coupling-gate part planted at +0.03 in every world."""
    rng = np.random.default_rng(seed); data = {}
    for name, _ in circuits():
        data[name] = {}
        for cp in (0, 1):
            if name == 'K': p1 = 0.01
            else:
                x = 0.0
                if world in ('readout', 'both'): x += 0.06
                if world in ('system', 'both') and name.startswith('P2'): x += 0.06
                p1 = (1 - x) / 2
            probe = (rng.random(SHOTS) < p1).astype(int)
            sys_bits = (rng.random(SHOTS) < (0.6 if name.startswith('P2') else 0.0)).astype(int) * rng.integers(1, 4, SHOTS)
            data[name][cp] = (probe << 2) | sys_bits
    off = {}
    for var in T.VARIANTS:
        base = (0.06 if world in ('readout', 'both') else 0) + (0.06 if world in ('system', 'both') else 0) + 0.03
        off[var] = (base, base - 0.02, base + 0.02)
    return data, off


if __name__ == '__main__':
    if '--synth' in sys.argv:
        w = sys.argv[sys.argv.index('--synth') + 1]; print(f'=== synthetic world {w}'); d, off = synth(w); analyze(d, 0.005, off); sys.exit(0)
    if '--analyze' in sys.argv:
        CHIP = sys.argv[sys.argv.index('--analyze') + 1]
        rec = json.load(open(f'{HERE}/job_split_{CHIP}.json')); z = np.load(f'{HERE}/raw_split_{CHIP}.npz')
        data = {}
        for i, m in enumerate(rec['meta']):
            data.setdefault(m['name'], {})[m['copy']] = z[f'c_{i}']
        props = json.load(open(f'{HERE}/props_split_{CHIP}.json')); pq = LAYOUT[CHIP][2]
        ro = [x['value'] for x in props['qubits'][pq] if x['name'] == 'prob_meas1_prep0'][0]
        off = json.load(open(f'{HERE}/offsets_{CHIP}.json'))
        print(f'=== {CHIP} job {rec["job_id"]}  probe q{pq}, stated prob_meas1_prep0 {ro:.4f}'); analyze(data, ro, {k: tuple(v) for k, v in off.items()}); sys.exit(0)
    CHIP = sys.argv[1]; FIRE = '--fire' in sys.argv
    s = service(); b = s.backend(CHIP); cs = circuits()
    items = [(n, c, cp) for cp in (0, 1) for n, c in cs]; random.Random(SEED).shuffle(items)
    tq = transpile([c for _, c, _ in items], b, optimization_level=1, initial_layout=LAYOUT[CHIP], seed_transpiler=7)
    dur = np.array([c.estimate_duration(b.target, 's') for c in tq]); est = float(SHOTS * np.sum(dur + 290e-6)); rem = s.usage().get('usage_remaining_seconds')
    print(f'{CHIP}: layout {LAYOUT[CHIP]}, {len(tq)} circuits, shots {SHOTS}, ESTIMATED {est:.0f} s, remaining {rem}')
    if not FIRE: sys.exit(0)
    if rem is not None and rem < est + 3: print('STOP: not enough seconds'); sys.exit(1)
    from qiskit_ibm_runtime import SamplerV2
    json.dump(b.properties().to_dict(), open(f'{HERE}/props_split_{CHIP}.json', 'w'), default=str)
    job = SamplerV2(mode=b).run(tq, shots=SHOTS)
    rec = dict(chip=CHIP, job_id=job.job_id(), submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), layout=LAYOUT[CHIP],
               meta=[dict(name=n, copy=cp) for n, _, cp in items], shots=SHOTS, keyfile=KEYFILE, instance=INSTANCE)
    json.dump(rec, open(f'{HERE}/job_split_{CHIP}.json', 'w')); print('SUBMITTED', CHIP, job.job_id(), flush=True)
