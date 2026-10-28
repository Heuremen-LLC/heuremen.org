"""Fresh-clock run, 22 Oct 2026 (Heuremen LLC). PLAN_FRESHCLOCK.md, sha256 in PLAN_FRESHCLOCK.sha256. Uses runs.py's fire path unchanged.
Control-run layout: 39 measured qubits E, 117 untouched spectators S. Spectators parked in 0 (no gate) while E does a PRELOAD of
0, 27 or 72 readouts, or 72 equal idles. Then SX on S (a fresh clock) and a WINDOW of 9 steps: readout (M) or equal idle (D). S read in X, Y.
Usage: python freshclock.py CHIP [--fire]"""
import sys, json
from qiskit import QuantumCircuit, ClassicalRegister
import runs, run as base
PRELOADS = (('M', 0), ('M', 27), ('M', 72), ('D', 72))
WINDOW = 9
def circuit(n, E, S, pre_mode, P, win_mode, st):
    L = base.L_DT; qc = QuantumCircuit(n)
    regs = [ClassicalRegister(len(E), 'p%d' % k) for k in range(P if pre_mode == 'M' else 0)] + [ClassicalRegister(len(E), 'w%d' % k) for k in range(WINDOW if win_mode == 'M' else 0)]
    for r in regs: qc.add_register(r)
    fin = ClassicalRegister(len(S), 'fin'); qc.add_register(fin)
    ri = 0
    for k in range(P):
        if pre_mode == 'M': qc.measure(E, regs[ri]); ri += 1
        else: qc.delay(L, E, unit='dt')
        qc.barrier()
    qc.sx(S); qc.barrier()
    for k in range(WINDOW):
        if win_mode == 'M': qc.measure(E, regs[ri]); ri += 1
        else: qc.delay(L, E, unit='dt')
        qc.barrier()
    for q in S: base.rot(qc, q, st[0])
    qc.measure(S, fin); return qc
def build_freshclock(b, chip):
    n = b.num_qubits; E = sorted(json.load(open(base.ROLES % chip))['emitters']); S = [q for q in range(n) if q not in set(E)]
    runs.set_L(b, E[0]); pubs, meta = [], []
    for pm, P in PRELOADS:
        for wm in ('D', 'M'):
            for st in ('XX', 'YY'):
                pubs.append(circuit(n, E, S, pm, P, wm, st)); meta.append(dict(pre='%s%d' % (pm, P), win=wm, setting=st))
    return pubs, meta, 2000, dict(E=E, S=S, L_dt=base.L_DT, window=WINDOW)
runs.BUILD['freshclock'] = build_freshclock
if __name__ == '__main__':
    sys.argv = [sys.argv[0], 'freshclock'] + sys.argv[1:]
    src = open(runs.__file__).read().split("if __name__ == '__main__':")[1]
    src = src.replace("open(HERE + '/PLAN.sha256')", "open(HERE + '/PLAN_FRESHCLOCK.sha256')")
    import textwrap; exec(textwrap.dedent(src), runs.__dict__)
