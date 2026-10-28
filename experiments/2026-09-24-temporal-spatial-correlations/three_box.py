#!/usr/bin/env python3
"""three_box.py -- 2026-09-24. Three-box paradox measured with a weak probe on IBM hardware.
The Aharonov-Vaidman THREE-BOX PARADOX with a weak probe:
  boxes on two qubits: A=|00>, B=|01>, C=|10>; the 4th corner |11> is never populated.
  PRE-select   |psi> = (A + B + C)/sqrt3
  WEAK PROBE   ancilla p rotated RY(2 eps) only if the system is in the box being checked (one box per circuit)
  POST-select  |phi> = (A + B - C)/sqrt3          (keep only shots whose final readout says phi)
  weak value of "the particle is in box X" ~ <X_probe>/(2 eps) on kept shots.  TEXTBOOK: A=+1, B=+1, C=-1, corner 11 = 0.
usage: three_box.py sim|real [keyfile] [chip]"""
import json, os, sys, math
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import RYGate

EPS = 0.25
SHOTS = 8192
BOXES = {"A": (0, 0), "B": (0, 1), "C": (1, 0), "corner_11": (1, 1)}      # (s1, s0)
OUTNAME = "three_box"
TEXTBOOK = {"A": 1, "B": 1, "C": -1, "corner_11": 0}


def prep(qc, s1, s0):
    qc.ry(2 * math.asin(1 / math.sqrt(3)), s1)          # s1=1 w.p. 1/3
    qc.x(s1); qc.ch(s1, s0); qc.x(s1)                    # if s1=0: s0 -> |+>


def unprep(qc, s1, s0):
    qc.x(s1); qc.ch(s1, s0); qc.x(s1)
    qc.ry(-2 * math.asin(1 / math.sqrt(3)), s1)


def circuit(box):
    qc = QuantumCircuit(3, 3)                            # q0 = s0, q1 = s1, q2 = probe
    s0, s1, p = 0, 1, 2
    prep(qc, s1, s0)
    b1, b0 = BOXES[box]
    if b1 == 0: qc.x(s1)
    if b0 == 0: qc.x(s0)
    qc.append(RYGate(2 * EPS).control(2), [s1, s0, p])  # rotate probe only if system is in this box
    if b1 == 0: qc.x(s1)
    if b0 == 0: qc.x(s0)
    qc.z(s1)                                             # |phi> = Z_s1 |psi>, so undo: Z_s1 then prep^dagger
    unprep(qc, s1, s0)
    qc.h(p)                                              # probe read in X
    qc.measure([s0, s1, p], [0, 1, 2])
    return qc


def weak_values(counts_by_box):
    out = {}
    for box, bits in counts_by_box.items():
        keep = (bits[:, 0] == 0) & (bits[:, 1] == 0)     # post-selection: system back to |00> after unprep = was |phi>
        xp = 1 - 2 * bits[keep, 2].astype(int)
        out[box] = {"kept": int(keep.sum()), "frac": float(keep.mean()), "Xp": float(xp.mean()) if keep.any() else float("nan"),
                    "weak": float(xp.mean() / (2 * EPS)) if keep.any() else float("nan"),
                    "se": float(xp.std(ddof=1) / math.sqrt(keep.sum()) / (2 * EPS)) if keep.sum() > 1 else float("nan")}
    return out


def run_on(be, sm_cls):
    pubs = [(transpile(circuit(b), backend=be, optimization_level=1),) for b in BOXES]
    sm = sm_cls(mode=be); sm.options.default_shots = SHOTS
    job = sm.run(pubs); jid = job.job_id() if hasattr(job, "job_id") else None
    res = job.result(); cb = {}
    for b, r in zip(BOXES, res):
        raw = np.array(r.data.c.to_bool_array())          # columns: clbit2, clbit1, clbit0
        cb[b] = np.column_stack([raw[:, 2], raw[:, 1], raw[:, 0]]).astype(int)   # -> [s0, s1, probe]
    return jid, weak_values(cb)


def main(mode, keyfile, chip):
    from qiskit_ibm_runtime import SamplerV2
    from qiskit_aer import AerSimulator
    ideal = run_on(AerSimulator(), SamplerV2)[1]          # finite-eps ideal reference (what perfect hardware gives at this eps)
    if mode == "sim":
        from qiskit_ibm_runtime.fake_provider import FakeMarrakesh; jid, hw = run_on(FakeMarrakesh(), SamplerV2); name = "fake_marrakesh"
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        be = QiskitRuntimeService(channel="ibm_quantum_platform", token=open(os.path.expanduser(keyfile)).read().strip()).backend(chip)
        jid, hw = run_on(be, SamplerV2); name = be.name
    print(f"job {jid}   eps={EPS}   (post-selection keeps ~1/9 of shots)")
    print("box        textbook   ideal(eps)   HARDWARE weak value      kept")
    tb = TEXTBOOK
    for b in BOXES:
        print(f"{b:10s} {tb[b]:+d}         {ideal[b]['weak']:+.3f}       {hw[b]['weak']:+.3f} +/- {hw[b]['se']:.3f}      {hw[b]['kept']}")
    s = sum(hw[b]["weak"] for b in BOXES if TEXTBOOK[b] != 0)    # the three real boxes: should sum to one particle
    print(f"sum of the three real boxes (one particle) = {s:+.3f}")
    json.dump({"backend": name, "job": jid, "eps": EPS, "shots": SHOTS, "ideal": ideal, "hardware": hw},
              open(os.path.expanduser(f"./runs/{OUTNAME}_{mode}_{name}.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sim", sys.argv[2] if len(sys.argv) > 2 else "~/.ibm_quantum_key",
         sys.argv[3] if len(sys.argv) > 3 else "ibm_marrakesh")
