#!/usr/bin/env python3
"""arrow_panel.py -- 2026-09-24. Reversible rotation vs mid-circuit measurement, and two-measurement correlation forward vs backward.
One qubit prepared in |+> (its X-coherence reads ~1):
  REVERSIBLE : RX(theta) then RX(-theta)  -> read X   (rotation + inverse rotation: coherence returns, ~1)
  LOOK       : mid-circuit Z measurement  -> read X   (one measurement: coherence gone, ~0; no later gate undoes it)
  LINK       : measure, RX(+theta), measure  vs  measure, RX(-theta), measure   (two-measurement correlation in both directions:
               identical, cos(theta))
Row kinds keep their stored names (BASE, REV, LOOK, LINK+, LINK-); LOOK denotes the single mid-circuit measurement.
usage: arrow_panel.py sim|real [keyfile] [chip]"""
import json, os, sys, math
import numpy as np
from qiskit import QuantumCircuit, transpile, ClassicalRegister, QuantumRegister

THETAS = [math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi]
QUBITS = [1, 2, 133]
SHOTS = 2048


def c_base():
    qc = QuantumCircuit(1, 1); qc.h(0); qc.h(0); qc.measure(0, 0); return qc          # |+>, read X (H then Z-readout)


def c_rev(t):
    qc = QuantumCircuit(1, 1); qc.h(0); qc.rx(t, 0); qc.rx(-t, 0); qc.h(0); qc.measure(0, 0); return qc


def c_look():
    q = QuantumRegister(1); m1 = ClassicalRegister(1, "m1"); m = ClassicalRegister(1, "c")
    qc = QuantumCircuit(q, m1, m); qc.h(0); qc.measure(0, m1[0]); qc.h(0); qc.measure(0, m[0]); return qc


def c_link(t):
    q = QuantumRegister(1); m1 = ClassicalRegister(1, "m1"); m = ClassicalRegister(1, "c")
    qc = QuantumCircuit(q, m1, m); qc.h(0); qc.measure(0, m1[0]); qc.rx(t, 0); qc.measure(0, m[0]); return qc


def main(mode, keyfile, chip):
    from qiskit_ibm_runtime import SamplerV2
    if mode == "sim":
        from qiskit_ibm_runtime.fake_provider import FakeMarrakesh; be = FakeMarrakesh()
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        be = QiskitRuntimeService(channel="ibm_quantum_platform", token=open(os.path.expanduser(keyfile)).read().strip()).backend(chip)
    items = []
    for q in QUBITS:
        items.append((q, "BASE", None, c_base()))
        items += [(q, "REV", t, c_rev(t)) for t in THETAS]
        items.append((q, "LOOK", None, c_look()))
        items += [(q, "LINK+", t, c_link(t)) for t in THETAS] + [(q, "LINK-", t, c_link(-t)) for t in THETAS]
    pubs = [(transpile(c, backend=be, initial_layout=[q], optimization_level=0),) for q, _, _, c in items]
    sm = SamplerV2(mode=be); sm.options.default_shots = SHOTS
    job = sm.run(pubs); jid = job.job_id() if hasattr(job, "job_id") else None; print("job", jid, flush=True)
    res = job.result(); rows = []
    for (q, kind, t, _), r in zip(items, res):
        z = 1 - 2 * np.array(r.data.c.to_bool_array()).ravel().astype(int)
        if kind.startswith("LINK"):
            z1 = 1 - 2 * np.array(r.data.m1.to_bool_array()).ravel().astype(int); v = float(np.mean(z1 * z))
        else:
            v = float(np.mean(z))                                # <X> for BASE/REV/LOOK (LOOK = one mid-circuit measurement)
        rows.append({"qubit": q, "kind": kind, "theta": t, "value": v})
    out = os.path.expanduser(f"./runs/arrow_panel_{mode}_{be.name}.json")
    json.dump({"backend": be.name, "job": jid, "shots": SHOTS, "rows": rows}, open(out, "w"), indent=1)
    m = lambda k, t=None: np.mean([r["value"] for r in rows if r["kind"] == k and (t is None or r["theta"] == t)])
    print(f"\ncoherence <X>: start {m('BASE'):+.3f} | twist+untwist " + " ".join(f"{m('REV', t):+.3f}" for t in THETAS) + f" | after ONE MID-CIRCUIT MEASUREMENT {m('LOOK'):+.3f}")
    print("link (two-measurement correlation)  theta: forward / backward / cos")
    for t in THETAS: print(f"   {t:5.3f}: {m('LINK+', t):+.3f} / {m('LINK-', t):+.3f} / {math.cos(t):+.3f}")
    print("saved", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim", sys.argv[2] if len(sys.argv) > 2 else "~/.ibm_quantum_key",
         sys.argv[3] if len(sys.argv) > 3 else "ibm_marrakesh")
