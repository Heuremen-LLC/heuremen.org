"""tutorial_qaoa.py -- IBM's QAOA tutorial (quantum.cloud.ibm.com/docs/en/tutorials/quantum-approximate-optimization-algorithm), small
example, run as written, then checked against the forecast Qiskit gives from IBM's published figures (PLAN_TUTORIAL_QAOA.md).
As written: 5-node graph, edges (0,1),(0,2),(0,4),(1,2),(2,3),(3,4) weight 1; MaxCut cost Hamiltonian (to_ising of the QUBO: 0.5 Z_iZ_j
per edge); QAOAAnsatz(reps=2) + measure_all; generate_preset_pass_manager(optimization_level=3); init [pi/2, pi/2, pi, pi];
Estimator default_shots 1000, dynamical_decoupling XY4, twirling enable_gates, num_randomizations 'auto'; COBYLA tol 1e-2;
Sampler default_shots 10000 with the same suppression; answer = most likely bitstring (reversed to node order), cut value.
Disclosed difference: job mode instead of a Session (pay-as-you-go bills a Session's open time); backend named instead of least_busy.
  tutorial_qaoa.py CHIP --fire | --synth (the same code on IBM's local testing mode with Aer built from the chip's record)"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge')
KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'Darkwave Disco'
EDGES = [(0, 1), (0, 2), (0, 4), (1, 2), (2, 3), (3, 4)]; N = 5

def problem():
    from qiskit.quantum_info import SparsePauliOp
    from qiskit.circuit.library import QAOAAnsatz
    terms = []
    for a, b in EDGES:
        z = ['I'] * N; z[a] = 'Z'; z[b] = 'Z'; terms.append((''.join(z)[::-1], 0.5))
    H_ = SparsePauliOp.from_list(terms); qc = QAOAAnsatz(cost_operator=H_, reps=2); qc.measure_all(); return H_, qc
def cut(bits): return sum(1 for a, b in EDGES if bits[a] != bits[b])
OPT = max(cut([(k >> i) & 1 for i in range(N)]) for k in range(2 ** N))

def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)

def run_workflow(mode, backend, log):
    """the tutorial's steps 2-4, verbatim in substance; returns params, final counts, job ids"""
    from qiskit.transpiler import generate_preset_pass_manager
    from qiskit_ibm_runtime import EstimatorV2 as Estimator, SamplerV2 as Sampler
    from scipy.optimize import minimize
    H_, circuit = problem()
    pm = generate_preset_pass_manager(optimization_level=3, backend=backend); candidate = pm.run(circuit)
    init = [np.pi / 2, np.pi / 2, np.pi, np.pi]; vals, jobs = [], []
    est = Estimator(mode=mode); est.options.default_shots = 1000
    est.options.dynamical_decoupling.enable = True; est.options.dynamical_decoupling.sequence_type = 'XY4'
    est.options.twirling.enable_gates = True; est.options.twirling.num_randomizations = 'auto'
    def cost(params):
        isa_h = H_.apply_layout(candidate.layout); job = est.run([(candidate, isa_h, params)]); r = job.result()[0]
        c = float(r.data.evs); vals.append(c); jobs.append(job.job_id() if hasattr(job, 'job_id') else None); log(f'   eval {len(vals)}: cost {c:+.4f}'); return c
    res = minimize(cost, init, args=(), method='COBYLA', tol=1e-2)
    smp = Sampler(mode=mode); smp.options.default_shots = 10000
    smp.options.dynamical_decoupling.enable = True; smp.options.dynamical_decoupling.sequence_type = 'XY4'
    smp.options.twirling.enable_gates = True; smp.options.twirling.num_randomizations = 'auto'
    opt_circ = candidate.assign_parameters(res.x); job = smp.run([(opt_circ,)], shots=int(1e4)); counts = job.result()[0].data.meas.get_counts()
    jobs.append(job.job_id() if hasattr(job, 'job_id') else None)
    return dict(params=list(map(float, res.x)), nfev=int(res.nfev), costs=vals, counts=counts, jobs=jobs, layout=list(candidate.layout.final_index_layout()), ops={k: int(v) for k, v in candidate.count_ops().items()}), candidate

def to_nodes(key): return [int(c) for c in key[::-1]]   # tutorial: to_bitstring then reverse -> node order
def summarise(counts):
    tot = sum(counts.values()); ml = max(counts, key=counts.get); p_opt = sum(v for k, v in counts.items() if cut(to_nodes(k)) == OPT) / tot
    ev = sum(v * cut(to_nodes(k)) for k, v in counts.items()) / tot; return ml, cut(to_nodes(ml)), p_opt, ev

def forecast_final(backend, candidate, params, shots=200000):
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel
    sim = AerSimulator(noise_model=NoiseModel.from_backend(backend, thermal_relaxation=True))
    return sim.run(candidate.assign_parameters(params), shots=shots, seed_simulator=11).result().get_counts()

if __name__ == '__main__':
    a = sys.argv; out = []
    log = lambda m: (print(m, flush=True), out.append(m))
    if '--synth' in a:
        from qiskit_aer import AerSimulator
        s = service(); b = s.backend('ibm_kingston'); fake = AerSimulator.from_backend(b)
        r, cand = run_workflow(fake, fake, log)
        fc = forecast_final(b, cand, r['params'])
        for lab, c in (('SYNTH hardware (local testing mode)', r['counts']), ('forecast', fc)):
            ml, mc, po, ev = summarise(c); log(f'{lab}: most likely {ml} -> cut {mc} (optimum {OPT}) | P(optimal cut) {po:.3f} | mean cut {ev:.3f}')
        json.dump(dict(run=r), open(f'{HERE}/synth_tutorial_qaoa.json', 'w'), default=str); sys.exit(0)
    CHIP = a[-1]; s = service(); b = s.backend(CHIP)
    json.dump(b.properties().to_dict(), open(f'{HERE}/props_tqaoa_{CHIP}.json', 'w'), default=str)
    log(f'=== {CHIP}: IBM QAOA tutorial as written, started {time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, remaining {s.usage().get("usage_remaining_seconds")} s')
    r, cand = run_workflow(b, b, log)
    json.dump(r, open(f'{HERE}/run_tqaoa_{CHIP}.json', 'w'), default=str)
    fc = forecast_final(b, cand, r['params']); json.dump(fc, open(f'{HERE}/forecast_tqaoa_{CHIP}.json', 'w'))
    for lab, c in (('HARDWARE', r['counts']), ('FORECAST (Aer from IBM record, same params)', fc)):
        ml, mc, po, ev = summarise(c); log(f'{lab}: most likely {ml} -> cut {mc} (optimum {OPT}) | P(optimal cut) {po:.3f} | mean cut {ev:.3f}')
    log(f'nfev {r["nfev"]}, layout {r["layout"]}, ops {r["ops"]}, remaining {s.usage().get("usage_remaining_seconds")} s')
    open(f'{HERE}/results_tqaoa_{CHIP}.text', 'w').write('\n'.join(out) + '\n')
