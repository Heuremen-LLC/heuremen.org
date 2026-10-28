#!/usr/bin/env python3
"""space_time_same.py -- 2026-09-24. Spatial correlation (entangled pair) vs temporal correlation (one qubit measured twice): same cos(theta) form.
Measured on one chip, same pairs as partner_twist's balanced run (ibm_marrakesh):
  SPACE : Bell pair (|00>+|11>)/sqrt2; read a straight (Z), read b tilted by theta (RY(-theta) then Z). Correlation = cos(theta).
  TIME  : (from partner_twist, already measured) one qubit, measure, twist RX(theta), measure again. Correlation = cos(theta).
  BLOCH : a's own three directions <X>,<Y>,<Z> alone (|+>: length 1) vs entangled with b (length ~0: no net polarization of a by itself).
usage: space_time_same.py sim|real"""
import json, os, sys, math, glob
import numpy as np
from qiskit import QuantumCircuit, transpile

THETAS = [k * math.pi / 4 for k in range(8)]
PAIRS = [(1, 2), (133, 134)]
SHOTS = 2048


def space(theta):
    qc = QuantumCircuit(2, 2); qc.h(0); qc.cx(0, 1); qc.ry(-theta, 1); qc.measure([0, 1], [0, 1]); return qc


def bloch(ent, axis):
    qc = QuantumCircuit(2, 2); qc.h(0)
    if ent: qc.cx(0, 1)
    if axis == "X": qc.h(0)
    if axis == "Y": qc.sdg(0); qc.h(0)
    qc.measure([0, 1], [0, 1]); return qc


def main(mode, keyfile="~/.ibm_quantum_key", chip="ibm_marrakesh"):
    from qiskit_ibm_runtime import SamplerV2
    if mode == "sim":
        from qiskit_ibm_runtime.fake_provider import FakeMarrakesh; be = FakeMarrakesh()
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        be = QiskitRuntimeService(channel="ibm_quantum_platform", token=open(os.path.expanduser(keyfile)).read().strip()).backend(chip)
    global PAIRS
    if chip != "ibm_marrakesh":
        from entanglement_cost import pair_table, pick_pairs
        PAIRS = [(p["a"], p["b"]) for p in pick_pairs(pair_table(be), 2)]
    items = [(p, "space", t, None) for p in PAIRS for t in THETAS] + [(p, "bloch", ent, ax) for p in PAIRS for ent in (0, 1) for ax in "XYZ"]
    circs = [space(x[2]) if x[1] == "space" else bloch(x[2], x[3]) for x in items]
    pubs = [(transpile(c, backend=be, initial_layout=list(p), optimization_level=0),) for c, (p, *_ ) in zip(circs, items)]
    sm = SamplerV2(mode=be); sm.options.default_shots = SHOTS
    job = sm.run(pubs); print("job", job.job_id() if hasattr(job, "job_id") else None, flush=True); res = job.result()
    out = {"space": {}, "bloch": {}}
    for (p, kind, a, ax), r in zip(items, res):
        bits = np.array(r.data.c.to_bool_array())            # shots x 2, column order: clbit1, clbit0
        za, zb = 1 - 2 * bits[:, -1].astype(int), 1 - 2 * bits[:, -2].astype(int)
        if kind == "space": out["space"].setdefault(str(p), []).append((a, float(np.mean(za * zb))))
        else: out["bloch"].setdefault(str(p), {})[("ENT" if a else "ALONE") + "_" + ax] = float(np.mean(za))
    # time curve from partner_twist balanced run (ENT arm of those pairs, a's two-measurement correlation)
    tl = sorted(glob.glob(os.path.expanduser("./runs/partner_twist_real_ibm_marrakesh_*.json")))
    tcurve = {}
    if tl:
        d = json.load(open(tl[-1]))
        for p in PAIRS:
            tcurve[str(p)] = [(r["theta"], r["C"]) for r in sorted(d["rows"], key=lambda r: r["theta"]) if (r["a"], r["b"]) == p and r["kind"] == "ENT"]
    print("\ntheta   cos    SPACE(pair corr, b tilted)   TIME(one qubit, two measurements)")
    for i, t in enumerate(THETAS):
        sp = np.mean([dict(out["space"][str(p)])[t] for p in PAIRS])
        tm = np.mean([dict(tcurve[str(p)]).get(t, np.nan) for p in PAIRS if str(p) in tcurve]) if tcurve else float("nan")
        print(f"{t:5.3f}  {math.cos(t):+.3f}   {sp:+.3f}                        {tm:+.3f}")
    print("\nALONE vs ENTANGLED: a's own Bloch vector (length = size of the qubit's local polarization)")
    for p in PAIRS:
        b = out["bloch"][str(p)]
        la = math.sqrt(sum(b["ALONE_" + x] ** 2 for x in "XYZ")); le = math.sqrt(sum(b["ENT_" + x] ** 2 for x in "XYZ"))
        print(f"  pair {p}: alone  X{b['ALONE_X']:+.3f} Y{b['ALONE_Y']:+.3f} Z{b['ALONE_Z']:+.3f}  length {la:.3f}   |   entangled  X{b['ENT_X']:+.3f} Y{b['ENT_Y']:+.3f} Z{b['ENT_Z']:+.3f}  length {le:.3f}")
    json.dump({"backend": be.name, "pairs": PAIRS, "space": out["space"], "bloch": out["bloch"], "time_from": tl[-1] if tl else None},
              open(os.path.expanduser(f"./runs/space_time_same_{mode}_{be.name}.json"), "w"), indent=1)


if __name__ == "__main__":
    main(*(sys.argv[1:4] if len(sys.argv) > 1 else ["sim"]))
