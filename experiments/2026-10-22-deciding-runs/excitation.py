"""Claim 4 run, 22 Oct 2026 (Heuremen LLC). PLAN_EXCITATION.md, sha256 in PLAN_EXCITATION.sha256. Uses runs.py's fire path unchanged.
Control-run layout: 39 measured, 117 untouched spectators on the equator. Four modes: idle with measured qubits in 0 (D0), readout in 0 (M0),
idle in 1 (D1), readout in 1 (M1); N = 0, 4, 8; spectators read in X and Y. Every circuit twice per job. Usage: python excitation.py CHIP [--fire]"""
import sys, json
import runs, run as base
def build_excitation(b, chip):
    n = b.num_qubits; E = sorted(json.load(open(base.ROLES % chip))['emitters']); S = [q for q in range(n) if q not in set(E)]
    runs.set_L(b, E[0]); pubs, meta = [], []
    for e in '01':
        for N in (0, 4, 8):
            for mode in ('D', 'M'):
                if N == 0 and mode == 'M': continue
                for st in ('XX', 'YY'):
                    pubs.append(base.spectator_circuit(n, E, [], S, E, mode, N, st, e, False)); meta.append(dict(mode=mode + e, N=N, setting=st))
    return pubs, meta, 2000, dict(E=E, S=S, L_dt=base.L_DT)
runs.BUILD['excitation'] = build_excitation
if __name__ == '__main__':
    sys.argv = [sys.argv[0], 'excitation'] + sys.argv[1:]
    src = open(runs.__file__).read().split("if __name__ == '__main__':")[1]
    src = src.replace("open(HERE + '/PLAN.sha256')", "open(HERE + '/PLAN_EXCITATION.sha256')")
    import textwrap; exec(textwrap.dedent(src), runs.__dict__)
