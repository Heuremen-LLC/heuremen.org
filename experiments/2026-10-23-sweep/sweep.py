"""sweep.py -- every per-qubit and per-edge figure IBM publishes, measured the way IBM's documentation says IBM measures it, one job per
chip (PLAN_SWEEP.md). T1/T2 are in ../ladder_run and ../hahn_t2; zz in the direct-zz run; this job covers the rest:
  INIT    circuits I1 = [measure a; X; measure b] and I0 = [measure a; measure b], all qubits, shots back to back.
          'a' in I1 reads each qubit after the PREVIOUS shot ended in |1> (IBM's init_error definition, claims row 25); 'a' in I0 after
          one ended in |0>. b in I1 = prepare |1> and read (prob_meas0_prep1); a in I0 = prepare |0> and read (prob_meas1_prep0).
  RB1     simultaneous single-qubit Clifford RB on ALL qubits (claims row 60: SX error 'measured simultaneously on all qubits'),
          independent random Cliffords per qubit, lengths M1, then the inverse Clifford, measure. Error per Clifford (1-p)/2; per pulse
          = EPC / mean pulses per Clifford (SX, X counted as one pulse each; row 55: SX, ID and X errors assumed equal).
  RB2     per-edge RB in isolation (claims row 59): edges in batches whose nearest qubits are at least 3 links apart, layers of
          [random single-qubit Clifford on both qubits, then CZ], lengths M2, then the inverse two-qubit Clifford, measure both.
          Error per layer 3/4 (1-p): IBM's 2Q error 'includes the error of the layer of single-qubit Cliffords'.
  sweep.py CHIP | --fire CHIP | --fetch CHIP | --analyze CHIP | --synth WORLD"""
import os, sys, json, time, random, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); H = os.path.expanduser('~/bones-forge')
KEYFILE, INSTANCE = '.ibm_quantum_key_3', 'Darkwave Disco'
SEED = 20261024; M1 = [1, 150, 500, 1000]; S1 = 3; M2 = [1, 10, 30, 60, 100]; S2 = 3
SHOTS_RO = 4000; SHOTS_RB = 300

# ---------- the single-qubit Clifford group in rz/sx ----------
def clifford_table():
    from qiskit import QuantumCircuit, transpile
    from qiskit.quantum_info import Clifford, Operator
    els, ops, seen = [], [], []
    frontier = [QuantumCircuit(1)]
    while frontier and len(els) < 24:
        nxt = []
        for qc in frontier:
            c = Clifford(qc); key = c.tableau.tobytes()
            if key in seen: continue
            seen.append(key); els.append(c)
            t = transpile(qc, basis_gates=['rz', 'sx'], optimization_level=3)
            ops.append([(i.operation.name, snap(float(i.operation.params[0])) if i.operation.params else None) for i in t.data])
            for g in ('h', 's'):
                q2 = qc.copy(); getattr(q2, g)(0); nxt.append(q2)
        frontier = nxt
    assert len(els) == 24, len(els)
    idx = {c.tableau.tobytes(): i for i, c in enumerate(els)}
    mul = np.array([[idx[els[i].compose(els[j]).tableau.tobytes()] for j in range(24)] for i in range(24)])   # i then j
    inv = np.array([idx[els[i].adjoint().tableau.tobytes()] for i in range(24)])
    pulses = np.mean([sum(1 for n, _ in o if n in ('sx', 'x')) for o in ops])
    return els, ops, mul, inv, pulses

def snap(p):
    k = round(p / (np.pi / 2)); assert abs(p - k * np.pi / 2) < 1e-6, p; return float(k * np.pi / 2)

def put(qc, ops, q):
    for n, p in ops:
        if n == 'rz': qc.rz(p, q)
        elif n == 'sx': qc.sx(q)
        elif n == 'x': qc.x(q)

def batches(n, edges, adj):
    """greedy batches of edges whose nearest qubits are >= 3 links apart"""
    dist = {}
    for s in range(n):
        d = {s: 0}; fr = [s]
        while fr:
            nf = []
            for u in fr:
                for v in adj[u]:
                    if v not in d and d[u] < 3: d[v] = d[u] + 1; nf.append(v)
            fr = nf
        dist[s] = d
    left = list(edges); out = []
    while left:
        b = []; rest = []
        for e in left:
            if all(min(dist[x].get(y, 9) for x in e for y in f) >= 3 for f in b): b.append(e)
            else: rest.append(e)
        out.append(b); left = rest
    return out

def build(n, edges, rng):
    from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile
    from qiskit.quantum_info import Clifford
    els, ops, mul, inv, pulses = clifford_table()
    adj = {q: set() for q in range(n)}
    for a, b in edges: adj[a].add(b); adj[b].add(a)
    circs, meta = [], []
    Q = list(range(n))
    for cp in (0, 1):                      # INIT / readout, two copies
        for kind in ('I1', 'I0', 'IQ'):
            qr = QuantumRegister(n, 'q'); A = ClassicalRegister(n, 'a'); B = ClassicalRegister(n, 'b'); qc = QuantumCircuit(qr, A, B)
            qc.measure(Q, A)
            if kind in ('I1', 'IQ'): qc.x(Q)
            qc.measure(Q, B)
            if kind == 'IQ': qc.x(Q)            # IQ ends in |0>: the next shot's b reads |1> prepared after a |0>-ending shot
            circs.append(qc); meta.append(dict(kind=kind, copy=cp, shots=SHOTS_RO))
    for m in M1:                           # RB1, all qubits at once
        for s in range(S1):
            qc = QuantumCircuit(n, n); cur = np.zeros(n, int)
            for step in range(m):
                pick = rng.integers(0, 24, n)
                for q in Q: put(qc, ops[pick[q]], q); cur[q] = mul[cur[q], pick[q]]
                qc.barrier()
            for q in Q: put(qc, ops[inv[cur[q]]], q)
            qc.barrier(); qc.measure(Q, Q); circs.append(qc); meta.append(dict(kind='RB1', m=m, seed=s, shots=SHOTS_RB))
    layer = {}
    def L(a, b):
        if (a, b) not in layer:
            c2 = QuantumCircuit(2); c2.compose(els[a].to_circuit(), [0], inplace=True); c2.compose(els[b].to_circuit(), [1], inplace=True); c2.cz(0, 1)
            layer[(a, b)] = Clifford(c2)
        return layer[(a, b)]
    invcache = {}
    def inv2(c):
        k = c.tableau.tobytes()
        if k not in invcache:
            t = transpile(c.adjoint().to_circuit(), basis_gates=['rz', 'sx', 'cz'], coupling_map=[[0, 1], [1, 0]], optimization_level=1)
            invcache[k] = [(i.operation.name, [t.find_bit(x).index for x in i.qubits], snap(float(i.operation.params[0])) if i.operation.params else None) for i in t.data]
        return invcache[k]
    B2 = batches(n, edges, adj)
    for bi, bat in enumerate(B2):          # RB2, isolated batches
        for m in M2:
            for s in range(S2):
                qs = [q for e in bat for q in e]; qc = QuantumCircuit(n, len(qs))
                cur = {e: Clifford(QuantumCircuit(2)) for e in bat}
                for step in range(m):
                    for e in bat:
                        a, b = rng.integers(0, 24, 2); put(qc, ops[a], e[0]); put(qc, ops[b], e[1]); qc.cz(e[0], e[1])
                        cur[e] = cur[e].compose(L(a, b))
                    qc.barrier(qs)
                for e in bat:
                    for nm, qq, p in inv2(cur[e]):
                        tq = [e[x] for x in qq]
                        if nm == 'rz': qc.rz(p, tq[0])
                        elif nm == 'sx': qc.sx(tq[0])
                        elif nm == 'cz': qc.cz(tq[0], tq[1])
                qc.barrier(qs); qc.measure(qs, list(range(len(qs))))
                circs.append(qc); meta.append(dict(kind='RB2', batch=bi, m=m, seed=s, shots=SHOTS_RB, edges=[list(e) for e in bat]))
    return circs, meta, pulses, B2

def service():
    from qiskit_ibm_runtime import QiskitRuntimeService
    tok = open(H + '/' + KEYFILE).read().strip(); s = QiskitRuntimeService(channel='ibm_quantum_platform', token=tok)
    crn = [i['crn'] for i in s.instances() if i.get('name') == INSTANCE][0]
    return QiskitRuntimeService(channel='ibm_quantum_platform', token=tok, instance=crn)

def chip_edges(b):
    es = set()
    for a, c in b.coupling_map.get_edges():
        if ('cz' in b.target.operation_names) and b.target.instruction_supported('cz', (a, c)): es.add(tuple(sorted((a, c))))
    return sorted(es)

# ---------- analysis ----------
def rbfit(m, y, B):
    """A p^m + B with B fixed at the fully mixed survival (1/2 for one qubit, 1/4 for two)"""
    from scipy.optimize import curve_fit
    m = np.asarray(m, float); y = np.asarray(y, float)
    try:
        (A, p), _ = curve_fit(lambda x, A, p: A * p ** x + B, m, y, p0=[1 - B, 0.999], bounds=([0, 0.5], [1, 1]), maxfev=20000)
        return p, A, B
    except Exception: return np.nan, np.nan, np.nan

def med_ci(x, rng, B=4000):
    x = np.asarray([v for v in x if np.isfinite(v)])
    if len(x) < 8: return np.nan, np.nan, np.nan, len(x)
    m = float(np.median(x)); lo, hi = np.percentile([np.median(rng.choice(x, len(x))) for _ in range(B)], [2.5, 97.5]); return m, lo, hi, len(x)
def four(lo, hi, ref=1.0, bad=False):
    if bad or not np.isfinite(lo): return 'CANNOT_TELL'
    return 'ABOVE' if lo > ref else 'BELOW' if hi < ref else 'CONSISTENT'

def stated(props):
    q = {}
    for i, ps in enumerate(props['qubits']):
        q[i] = {x['name']: x['value'] for x in ps}
    g = {}
    for x in props['gates']:
        e = {y['name']: y['value'] for y in x['parameters']}
        g[(x['gate'], tuple(x['qubits']))] = e.get('gate_error')
    return q, g

def analyze(meta, data, n, pulses, st_q, st_g, label, rng=None):
    """data[i]: dict reg -> bool array (shots, bits). Returns readings; prints them. ratios are delivered / stated."""
    rng = rng or np.random.default_rng(3); res = {}
    # INIT and readout
    acc = {}
    for i, m in enumerate(meta):
        if m['kind'] in ('I1', 'I0', 'IQ'):
            for r in ('a', 'b'): acc.setdefault((m['kind'], r), []).append(data[i][r].mean(axis=0))
    P = {k: np.mean(v, axis=0) for k, v in acc.items()}
    cg = np.median(np.abs(np.array(acc[('I1', 'a')][0]) - np.array(acc[('I1', 'a')][1])))
    a1, a0 = P[('I1', 'a')], P[('I0', 'a')]; p01 = 1 - P[('IQ', 'b')]; p01_1 = 1 - P[('I1', 'b')]
    lo_init = (a1 - a0) / np.clip(1 - a0 - p01, 0.5, 1); hi_init = a1 / np.clip(1 - p01, 0.5, 1)
    sd = np.sqrt(a1 * (1 - a1) / (2 * SHOTS_RO) + a0 * (1 - a0) / (2 * SHOTS_RO))
    cnt = {'ABOVE': 0, 'BELOW': 0, 'CONSISTENT': 0, 'CANNOT_TELL': 0}; rows = []
    for q in range(n):
        s = st_q.get(q, {}).get('init_error')
        if s is None: cnt['CANNOT_TELL'] += 1; continue
        o = 'ABOVE' if lo_init[q] - 1.96 * sd[q] > s else 'BELOW' if hi_init[q] + 1.96 * sd[q] < s else 'CONSISTENT'; cnt[o] += 1; rows.append((q, s, lo_init[q], hi_init[q], o))
    res['INIT_counts'] = cnt; res['INIT_rows'] = rows
    print(f'-- {label}')
    print(f'   INIT per qubit, delivered |1> after a |1>-ending shot, bounds [{np.median(lo_init):.4f}, {np.median(hi_init):.4f}] (chip medians) vs stated init_error median {np.median([st_q[q].get("init_error", np.nan) for q in st_q]):.4f}: {cnt}')
    # IBM does not say what the shot before its readout shot ended in, so each delivered readout figure is a BRACKET:
    # previous shot ended in |0> .. ended in |1>. Per qubit: how far the stated value sits OUTSIDE the bracket (0 inside), widened by shot noise.
    sh = lambda x: 1.96 * np.sqrt(np.clip(x * (1 - x), 1e-6, None) / (2 * SHOTS_RO))
    br = {'prob_meas1_prep0': (np.minimum(a0, a1) - sh(np.minimum(a0, a1)), np.maximum(a0, a1) + sh(np.maximum(a0, a1))),
          'prob_meas0_prep1': (np.minimum(p01, p01_1) - sh(np.minimum(p01, p01_1)), np.maximum(p01, p01_1) + sh(np.maximum(p01, p01_1)))}
    br['readout_error'] = ((br['prob_meas1_prep0'][0] + br['prob_meas0_prep1'][0]) / 2, (br['prob_meas1_prep0'][1] + br['prob_meas0_prep1'][1]) / 2)
    for key, (L_, H_) in br.items():
        d = [st_q[q][key] - np.clip(st_q[q][key], L_[q], H_[q]) for q in range(n) if key in st_q.get(q, {})]
        m_, lo, hi, k = med_ci(d, rng); o = 'STATED ABOVE DELIVERED' if lo > 0 else 'STATED BELOW DELIVERED' if hi < 0 else 'CONSISTENT'; res[key] = (m_, lo, hi, o)
        nb = sum(1 for x in d if x < 0); na = sum(1 for x in d if x > 0)
        print(f'   {key}: stated outside the delivered bracket, median {m_:+.4f} [{lo:+.4f}, {hi:+.4f}] n {k} -> {o}   (qubits stated below bracket {nb}, above {na}; bracket chip medians {np.median(L_):.4f}..{np.median(H_):.4f})')
    print(f'   copy gap median |a1 copy0 - copy1| {cg:.4f}')
    # RB1
    surv = {}
    for i, m in enumerate(meta):
        if m['kind'] == 'RB1': surv.setdefault(m['m'], []).append(1 - data[i]['c'].mean(axis=0))
    ms = sorted(surv); Y = np.array([np.mean(surv[m], axis=0) for m in ms])
    epg = {}
    for q in range(n):
        p, A, Bv = rbfit(ms, Y[:, q], 0.5)
        if np.isfinite(p) and A > 0.2: epg[q] = (1 - p) / 2 / pulses
    r = [epg[q] / st_g[('sx', (q,))] for q in epg if st_g.get(('sx', (q,))) not in (None, 0)]
    m_, lo, hi, k = med_ci(r, rng); o = four(lo, hi); res['sx'] = (m_, lo, hi, o)
    print(f'   RB1 sx error per pulse delivered/stated, median {m_:.3f} [{lo:.3f}, {hi:.3f}] n {k} -> {o}   (delivered chip median {np.median(list(epg.values())):.2e}; mean pulses per Clifford {pulses:.3f}; m=1 survival {np.median(Y[0]):.3f})')
    zero_sx = sum(1 for q in range(n) if st_g.get(('sx', (q,))) == 0); one_sx = sum(1 for q in range(n) if st_g.get(('sx', (q,))) == 1)
    # RB2
    s2 = {}
    for i, m in enumerate(meta):
        if m['kind'] == 'RB2':
            bits = data[i]['c']
            for j, e in enumerate(m['edges']):
                ok = ~(bits[:, 2 * j] | bits[:, 2 * j + 1]); s2.setdefault(tuple(e), {}).setdefault(m['m'], []).append(ok.mean())
    epl = {}
    for e, d in s2.items():
        mm = sorted(d); p, A, Bv = rbfit(mm, [np.mean(d[x]) for x in mm], 0.25)
        if np.isfinite(p) and A > 0.2: epl[e] = 0.75 * (1 - p)
    def sg(e): return st_g.get(('cz', e), st_g.get(('cz', e[::-1])))
    r = [epl[e] / sg(e) for e in epl if sg(e) not in (None, 0, 1)]
    m_, lo, hi, k = med_ci(r, rng); o = four(lo, hi); res['cz'] = (m_, lo, hi, o)
    z = [(e, epl.get(e)) for e in s2 if sg(e) == 0]; one = [(e, epl.get(e)) for e in s2 if sg(e) == 1]
    print(f'   RB2 CZ layer error delivered/stated, median {m_:.3f} [{lo:.3f}, {hi:.3f}] n {k} of {len(s2)} edges -> {o}   (delivered chip median {np.median(list(epl.values())) if epl else np.nan:.2e})')
    print(f'   edges stated cz error = 0: {len(z)}, delivered ' + ', '.join(f'{e[0]}-{e[1]} {v:.2e}' if v is not None else f'{e[0]}-{e[1]} no fit' for e, v in z[:12]))
    print(f'   edges stated cz error = 1 (IBM: undefined): {len(one)}, delivered ' + ', '.join(f'{e[0]}-{e[1]} {v:.2e}' if v is not None else f'{e[0]}-{e[1]} no fit' for e, v in one[:12]))
    res['zero_cz'] = z; res['one_cz'] = one; res['epl'] = epl; res['epg'] = epg
    # tails: how many single units sit far from stated (x3 or /3) -- descriptive
    far1 = sum(1 for q in epg if st_g.get(('sx', (q,))) and (epg[q] / st_g[('sx', (q,))] > 3 or epg[q] / st_g[('sx', (q,))] < 1 / 3))
    far2 = sum(1 for e in epl if sg(e) not in (None, 0, 1) and (epl[e] / sg(e) > 3 or epl[e] / sg(e) < 1 / 3))
    print(f'   far from stated (x3 either way): sx {far1} of {len(epg)} qubits, cz {far2} of {len(epl)} edges')
    return res

# ---------- synthetic: Aer stabilizer on a heavy-hex piece, every figure planted ----------
def synth(world, seed=7):
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
    rng = np.random.default_rng(seed)
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (2, 12), (6, 13), (10, 14), (12, 15), (13, 16), (14, 17), (15, 18), (18, 19), (19, 20), (20, 16), (16, 21), (21, 22), (22, 17)]
    n = 23; f = {'match': 1, 'worse': 2.5, 'better': 0.4}.get(world, 1)
    # truth: init after a |1>-ending shot IE, after a |0>-ending shot TH, assignment errors E01 (0 read as 1), E10 (1 read as 0).
    # stated: what IBM would print if its readout shots follow a 50/50 mix of previous states (inside our bracket by construction).
    ie_s = rng.uniform(0.002, 0.02, n); E01 = rng.uniform(0.002, 0.015, n); E10 = rng.uniform(0.008, 0.03, n); TH = np.full(n, 0.0005)
    IE = ie_s * {'inithigh': 6, 'initlow': 0.3}.get(world, 1)
    if world == 'rohigh': E01 = E01 * 3; E10 = E10 * 3
    E01s = E01 / (3 if world == 'rohigh' else 1); E10s = E10 / (3 if world == 'rohigh' else 1)
    mix = 0.5 * (ie_s + TH)
    st_q = {q: {'init_error': float(ie_s[q]), 'prob_meas1_prep0': float(E01s[q] + mix[q] * (1 - E01s[q] - E10s[q])),
                'prob_meas0_prep1': float(E10s[q] + mix[q] * (1 - E01s[q] - E10s[q]))} for q in range(n)}
    for q in st_q: st_q[q]['readout_error'] = (st_q[q]['prob_meas1_prep0'] + st_q[q]['prob_meas0_prep1']) / 2
    st_g = {('sx', (q,)): float(rng.uniform(1.5e-4, 6e-4)) for q in range(n)}
    for e in edges: st_g[('cz', e)] = float(rng.uniform(2e-3, 8e-3))
    st_g[('cz', edges[3])] = 0.0                                   # an edge IBM states as 0: delivered must still read its true error
    true_cz = {e: (st_g[('cz', e)] or 4e-3) * f for e in edges}
    nm = NoiseModel()
    for q in range(n):
        nm.add_quantum_error(depolarizing_error(2 * st_g[('sx', (q,))] * f, 1), ['sx', 'x'], [q])   # EPG = 1/2 p for 1q depolarizing (d=2)
    _, _, _, _, pulses = clifford_table()
    for e in edges:   # layer error 3/4 p_layer; put (layer - 2 x 1q layer part) on CZ so the LAYER matches the planted value
        oneq = 2 * pulses * np.mean([st_g[('sx', (q,))] * f for q in e]) * 1.0
        pz = max(true_cz[e] - oneq, 1e-5) * 4 / 3
        nm.add_quantum_error(depolarizing_error(pz, 2), ['cz'], [e[0], e[1]]); nm.add_quantum_error(depolarizing_error(pz, 2), ['cz'], [e[1], e[0]])
    sim = AerSimulator(method='stabilizer', noise_model=nm)
    circs, meta, pulses, B2 = build(n, edges, np.random.default_rng(seed))
    data = []
    for c, m in zip(circs, meta):
        if m['kind'] in ('I1', 'I0', 'IQ'):   # shot-to-shot memory is not simulable in Aer: planted directly
            sh_ = m['shots']; prior1 = m['kind'] == 'I1'
            pop = rng.random((sh_, n)) < (IE if prior1 else TH)
            rd = lambda state: np.where(state, rng.random((sh_, n)) > E10, rng.random((sh_, n)) < E01)
            a = rd(pop); b = rd(~pop) if m['kind'] in ('I1', 'IQ') else rd(pop)
            data.append({'a': a, 'b': b})
        else:
            r = sim.run(c, shots=m['shots'], seed_simulator=int(rng.integers(1 << 30)), memory=True).result().get_memory()
            data.append({'c': np.array([[ch == '1' for ch in s[::-1]] for s in r])})
    st_q_used = st_q
    analyze(meta, data, n, pulses, st_q_used, st_g, f'synth {world}: planted sx and cz = {f} x stated' + (', init 0.3x' if world == 'initlow' else ', init 6x' if world == 'inithigh' else ', readout 3x' if world == 'rohigh' else ''))

if __name__ == '__main__':
    a = sys.argv
    if '--synth' in a: synth(a[a.index('--synth') + 1]); sys.exit(0)
    if '--fetch' in a:
        chip = a[a.index('--fetch') + 1]; rec = json.load(open(f'{HERE}/job_sweep_{chip}.json')); res = service().job(rec['job_id']).result(); out = {}
        for i, pr in enumerate(res):
            for r in ('a', 'b', 'c'):
                if hasattr(pr.data, r): out[f'{r}_{i}'] = getattr(pr.data, r).to_bool_array(order='little')
        np.savez_compressed(f'{HERE}/raw_sweep_{chip}.npz', **out); print('saved', chip, len(res), 'pubs'); sys.exit(0)
    if '--analyze' in a:
        chip = a[a.index('--analyze') + 1]; rec = json.load(open(f'{HERE}/job_sweep_{chip}.json')); z = np.load(f'{HERE}/raw_sweep_{chip}.npz')
        props = json.load(open(f'{HERE}/props_sweep_{chip}.json')); st_q, st_g = stated(props)
        data = [{r: z[f'{r}_{i}'] for r in ('a', 'b', 'c') if f'{r}_{i}' in z} for i in range(len(rec['meta']))]
        print(f'=== {chip} job {rec["job_id"]}; IBM record {props["last_update_date"]}; submitted {rec["submitted_utc"]}')
        analyze(rec['meta'], data, rec['n'], rec['pulses'], st_q, st_g, chip); sys.exit(0)
    CHIP = a[-1] if '--fire' not in a else a[a.index('--fire') + 1]; FIRE = '--fire' in a
    s = service(); b = s.backend(CHIP); n = b.num_qubits; edges = chip_edges(b)
    circs, meta, pulses, B2 = build(n, edges, np.random.default_rng(SEED))
    from qiskit import transpile
    tq = transpile(circs, b, optimization_level=0, initial_layout=list(range(n)), seed_transpiler=7)
    ops0 = sum(sum(c.count_ops().values()) for c in circs); ops1 = sum(sum(c.count_ops().values()) for c in tq)
    dur = np.array([c.estimate_duration(b.target, 's') if hasattr(c, 'estimate_duration') else 0 for c in tq]); shots = np.array([m['shots'] for m in meta])
    est = float(np.sum(shots * (dur + 290e-6))); rem = s.usage().get('usage_remaining_seconds')
    print(f'{CHIP}: {n} qubits, {len(edges)} CZ edges in {len(B2)} isolated batches, {len(tq)} circuits, ops {ops0} -> {ops1} after transpile (must be equal), ESTIMATED {est:.0f} s, remaining {rem}')
    if not FIRE: sys.exit(0)
    if ops0 != ops1: print('STOP: transpile changed the circuits'); sys.exit(1)
    if rem is not None and rem < est + 20: print('STOP: not enough seconds'); sys.exit(1)
    from qiskit_ibm_runtime import SamplerV2
    json.dump(b.properties().to_dict(), open(f'{HERE}/props_sweep_{CHIP}.json', 'w'), default=str)
    pubs = [(c, None, m['shots']) for c, m in zip(tq, meta)]
    job = SamplerV2(mode=b).run(pubs)
    json.dump(dict(chip=CHIP, job_id=job.job_id(), submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), meta=meta, n=n, pulses=pulses, shots_ro=SHOTS_RO, keyfile=KEYFILE, instance=INSTANCE), open(f'{HERE}/job_sweep_{CHIP}.json', 'w'))
    print('SUBMITTED', CHIP, job.job_id(), flush=True)
