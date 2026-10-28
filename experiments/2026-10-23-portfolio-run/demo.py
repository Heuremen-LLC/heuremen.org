"""demo.py -- does picking qubits by IBM's published error figures cost a customer-style workload? (PLAN_DEMO.md)
Workload: a 10-qubit QAOA-shaped circuit (MaxCut on a line, two layers of ZZ interactions and X mixers, fixed angles) followed by its exact
inverse behind a barrier, so the correct output is all zeros and success = P(0000000000) is known exactly (a mirror circuit).
Placements, same job, same compile settings (only the physical qubits differ):
  QISKIT  the layout Qiskit's preset pass manager (optimization level 2, the default a company gets) picks from IBM's published figures
  IBM     the best 10-qubit line under IBM's published sx, cz, readout figures, scored by one fixed formula
  OURS    the best 10-qubit line under OUR measured sx (per pulse), cz (layer), readout figures (sweep run 1 + repeat, same formula)
  KG_*    known-good: the same three placements with nothing but measurement (success = readout only)
Score (same formula for IBM and OURS): predicted success = prod over gates (1 - error) x prod over the 10 readouts (1 - readout error),
gate counts taken from the compiled circuit on a template line. Readings in the plan."""
import os, sys, json, time, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge')
sys.path.insert(0, os.path.join(HERE, '..', 'sweep_run'))
KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'Darkwave Disco'
NQ = 10; SHOTS = 4000; GAMMA = [0.71, 1.13]; BETA = [0.42, 0.27]

def workload(n=NQ):
    from qiskit import QuantumCircuit
    u = QuantumCircuit(n)
    u.h(range(n))
    for g, b in zip(GAMMA, BETA):
        for i in range(0, n - 1, 2): u.rzz(2 * g, i, i + 1)
        for i in range(1, n - 1, 2): u.rzz(2 * g, i, i + 1)
        u.rx(2 * b, range(n))
    qc = QuantumCircuit(n, n); qc.compose(u, inplace=True); qc.barrier(); qc.compose(u.inverse(), inplace=True); qc.barrier()
    qc.measure(range(n), range(n)); return qc

def known_good(n=NQ):
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(n, n); qc.measure(range(n), range(n)); return qc

def lines(edges, n_phys, L=NQ):
    adj = {q: set() for q in range(n_phys)}
    for a, b in edges: adj[a].add(b); adj[b].add(a)
    out = []
    def dfs(path):
        if len(path) == L:
            if path[0] < path[-1]: out.append(tuple(path))
            return
        for v in adj[path[-1]]:
            if v not in path: dfs(path + [v])
    for s in range(n_phys): dfs([s])
    return out

def counts_template(b, path):
    from qiskit import transpile
    t = transpile(workload(), b, optimization_level=1, initial_layout=list(path), seed_transpiler=7)
    sx = {i: 0 for i in range(NQ)}; cz = {i: 0 for i in range(NQ - 1)}
    pos = {q: i for i, q in enumerate(path)}
    for ins in t.data:
        qs = [t.find_bit(x).index for x in ins.qubits]
        if ins.operation.name in ('sx', 'x'): sx[pos[qs[0]]] += 1
        elif ins.operation.name == 'cz': i, j = sorted(pos[q] for q in qs); assert j == i + 1; cz[i] += 1
    return sx, cz

def predict(path, sxc, czc, e_sx, e_cz, e_ro):
    lp = 0.0
    for i, q in enumerate(path):
        lp += sxc[i] * np.log1p(-e_sx[q]) + np.log1p(-e_ro[q])
    for i in range(NQ - 1):
        a, b = path[i], path[i + 1]; e = e_cz.get((a, b), e_cz.get((b, a)))
        lp += czc[i] * np.log1p(-e)
    return float(np.exp(lp))

def best(paths, sxc, czc, e_sx, e_cz, e_ro):
    ok = [p for p in paths if all(q in e_sx and q in e_ro for q in p) and all(((p[i], p[i + 1]) in e_cz or (p[i + 1], p[i]) in e_cz) for i in range(NQ - 1))]
    sc = [(predict(p, sxc, czc, e_sx, e_cz, e_ro), p) for p in ok]; sc.sort(reverse=True); return sc[0], len(ok)

def ibm_errors(props):
    e_sx, e_cz, e_ro = {}, {}, {}
    for x in props['gates']:
        v = {y['name']: y['value'] for y in x['parameters']}.get('gate_error')
        if v is None or v >= 1: continue
        if x['gate'] == 'sx': e_sx[x['qubits'][0]] = v
        if x['gate'] == 'cz': e_cz[tuple(x['qubits'])] = v
    for q, ps in enumerate(props['qubits']):
        for y in ps:
            if y['name'] == 'readout_error' and y['value'] < 1: e_ro[q] = y['value']
    return e_sx, e_cz, e_ro

def our_errors(chip):
    """mean of sweep run 1 and repeat r2: sx per pulse, cz layer error, readout = midpoint of the delivered bracket"""
    import io, contextlib, sweep as S
    runs = []
    for sfx in ('', '_r2'):
        d = os.path.join(HERE, '..', 'sweep_run'); rec = json.load(open(f'{d}/job_sweep_{chip}{sfx}.json')); z = np.load(f'{d}/raw_sweep_{chip}{sfx}.npz')
        st_q, st_g = S.stated(json.load(open(f'{d}/props_sweep_{chip}{sfx}.json')))
        data = [{r: z[f'{r}_{i}'] for r in ('a', 'b', 'c') if f'{r}_{i}' in z} for i in range(len(rec['meta']))]
        with contextlib.redirect_stdout(io.StringIO()): res = S.analyze(rec['meta'], data, rec['n'], rec['pulses'], st_q, st_g, chip)
        acc = {}
        for i, m in enumerate(rec['meta']):
            if m['kind'] in ('I1', 'I0', 'IQ'):
                for r in ('a', 'b'): acc.setdefault((m['kind'], r), []).append(data[i][r].mean(axis=0))
        P = {k: np.mean(v, axis=0) for k, v in acc.items()}
        p10 = (P[('I0', 'a')] + P[('I1', 'a')]) / 2; p01 = ((1 - P[('IQ', 'b')]) + (1 - P[('I1', 'b')])) / 2
        runs.append((res['epg'], res['epl'], {q: (p10[q] + p01[q]) / 2 for q in range(rec['n'])}))
    e_sx = {q: np.mean([r[0][q] for r in runs if q in r[0]]) for q in set(runs[0][0]) | set(runs[1][0])}
    e_cz = {e: np.mean([r[1][e] for r in runs if e in r[1]]) for e in set(runs[0][1]) | set(runs[1][1])}
    e_ro = {q: np.mean([r[2][q] for r in runs]) for q in runs[0][2]}
    return e_sx, e_cz, e_ro

def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)

def choose(b, props, chip):
    from qiskit.transpiler import generate_preset_pass_manager
    n = b.num_qubits; edges = sorted({tuple(sorted(e)) for e in b.coupling_map.get_edges()}); paths = lines(edges, n)
    sxc, czc = counts_template(b, paths[0])
    ibm = ibm_errors(props); ours = our_errors(chip)
    (pi, p_ibm), n_ok_i = best(paths, sxc, czc, *ibm); (po, p_our), n_ok_o = best(paths, sxc, czc, *ours)
    pm = generate_preset_pass_manager(optimization_level=2, backend=b, seed_transpiler=7); t = pm.run(workload())
    p_qk = tuple(t.layout.initial_index_layout(filter_ancillas=True)[:NQ]) if t.layout else None
    return dict(paths=len(paths), sxc=sxc, czc=czc, IBM=list(p_ibm), OURS=list(p_our), QISKIT=list(p_qk),
                pred_IBM_by_ibm=pi, pred_OURS_by_ours=po,
                pred_IBM_by_ours=predict(p_ibm, sxc, czc, *ours), pred_OURS_by_ibm=predict(p_our, sxc, czc, *ibm),
                pred_QISKIT_by_ibm=predict(p_qk, sxc, czc, *ibm) if all(q in ibm[0] for q in p_qk) else None,
                pred_QISKIT_by_ours=predict(p_qk, sxc, czc, *ours) if all(q in ours[0] for q in p_qk) else None)

if __name__ == '__main__':
    a = sys.argv; CHIP = a[-1]; FIRE = '--fire' in a
    if '--fetch' in a:
        rec = json.load(open(f'{HERE}/job_demo_{CHIP}.json')); res = service().job(rec['job_id']).result()
        np.savez_compressed(f'{HERE}/raw_demo_{CHIP}.npz', **{f'c_{i}': pr.data.c.to_bool_array(order='little') for i, pr in enumerate(res)}); print('saved', CHIP); sys.exit(0)
    if '--analyze' in a:
        rec = json.load(open(f'{HERE}/job_demo_{CHIP}.json')); z = np.load(f'{HERE}/raw_demo_{CHIP}.npz'); ch = rec['choice']
        succ = {}
        for i, m in enumerate(rec['meta']): succ.setdefault(m['arm'], []).append((~z[f'c_{i}'].any(axis=1)).astype(float))
        print(f'=== {CHIP} job {rec["job_id"]} | placements IBM {ch["IBM"]} OURS {ch["OURS"]} QISKIT {ch["QISKIT"]}')
        S = {k: np.concatenate(v) for k, v in succ.items()}
        for k in sorted(S):
            cp = [v.mean() for v in succ[k]]; print(f'   {k:10} success {S[k].mean():.4f} (copies {cp[0]:.4f}, {cp[1]:.4f}; n {len(S[k])})')
        rng = np.random.default_rng(9)
        def ratio(x, y):
            r = x.mean() / y.mean(); bs = [rng.choice(x, len(x)).mean() / rng.choice(y, len(y)).mean() for _ in range(4000)]; lo, hi = np.percentile(bs, [2.5, 97.5])
            return r, lo, hi, 'ABOVE' if lo > 1 else 'BELOW' if hi < 1 else 'CONSISTENT'
        kg = all(S[f'KG_{k}'].mean() >= 0.80 for k in ('IBM', 'OURS', 'QISKIT'))
        print(f'   KNOWN-GOOD (readout-only success >= 0.80 on all three placements): {"MET" if kg else "NOT MET"}')
        for lab, x, y in (('R1 OURS / IBM (same search, only the numbers differ)', 'OURS', 'IBM'), ('R2 OURS / QISKIT default', 'OURS', 'QISKIT'), ('R3 IBM / QISKIT', 'IBM', 'QISKIT')):
            r, lo, hi, o = ratio(S[x], S[y]); print(f'   {lab}: {r:.3f} [{lo:.3f}, {hi:.3f}] -> {o if kg else "CANNOT_TELL"}')
        for arm, key_i, key_o in (('IBM', 'pred_IBM_by_ibm', 'pred_IBM_by_ours'), ('OURS', 'pred_OURS_by_ibm', 'pred_OURS_by_ours'), ('QISKIT', 'pred_QISKIT_by_ibm', 'pred_QISKIT_by_ours')):
            d = S[arm].mean(); se = np.sqrt(d * (1 - d) / len(S[arm]))
            print(f'   R4 {arm:6} delivered {d:.4f} +- {1.96 * se:.4f} | predicted from IBM figures {ch[key_i]:.4f} (delivered/predicted {d / ch[key_i]:.3f}) | from ours {ch[key_o]:.4f} ({d / ch[key_o]:.3f})')
        sys.exit(0)
    s = service(); b = s.backend(CHIP); props = b.properties().to_dict(); ch = choose(b, props, CHIP)
    print(f'{CHIP}: {ch["paths"]} candidate 10-qubit lines | IBM pick {ch["IBM"]} | OURS pick {ch["OURS"]} | QISKIT pick {ch["QISKIT"]}')
    print(f'   predicted success  IBM line: {ch["pred_IBM_by_ibm"]:.4f} by IBM figures, {ch["pred_IBM_by_ours"]:.4f} by ours | OURS line: {ch["pred_OURS_by_ibm"]:.4f} by IBM, {ch["pred_OURS_by_ours"]:.4f} by ours')
    from qiskit import transpile
    items = []
    for arm in ('IBM', 'OURS', 'QISKIT'):
        for cp in (0, 1):
            items.append((transpile(workload(), b, optimization_level=1, initial_layout=ch[arm], seed_transpiler=7), dict(arm=arm, copy=cp)))
            items.append((transpile(known_good(), b, optimization_level=1, initial_layout=ch[arm], seed_transpiler=7), dict(arm='KG_' + arm, copy=cp)))
    ops = {m['arm']: dict(c.count_ops()) for c, m in items if m['copy'] == 0}
    print('   compiled ops per arm:', {k: {g: v for g, v in o.items() if g in ('sx', 'x', 'cz')} for k, o in ops.items() if not k.startswith('KG')})
    dur = np.array([c.estimate_duration(b.target, 's') for c, _ in items]); est = float(SHOTS * np.sum(dur + 290e-6)); rem = s.usage().get('usage_remaining_seconds')
    print(f'   {len(items)} circuits, ESTIMATED {est:.0f} s, remaining {rem}')
    if not FIRE: sys.exit(0)
    if rem is not None and rem < est + 10: print('STOP: not enough seconds'); sys.exit(1)
    from qiskit_ibm_runtime import SamplerV2
    json.dump(props, open(f'{HERE}/props_demo_{CHIP}.json', 'w'), default=str)
    job = SamplerV2(mode=b).run([c for c, _ in items], shots=SHOTS)
    json.dump(dict(chip=CHIP, job_id=job.job_id(), submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), meta=[m for _, m in items], choice=ch), open(f'{HERE}/job_demo_{CHIP}.json', 'w'), default=str)
    print('SUBMITTED', CHIP, job.job_id(), flush=True)
