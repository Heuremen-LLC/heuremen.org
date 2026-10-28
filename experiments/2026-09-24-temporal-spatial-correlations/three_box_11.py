#!/usr/bin/env python3
"""three_box_11.py -- 2026-09-24. Same three-box paradox, with box C placed at the |11> corner of the two-qubit register.
  boxes: A=|00>, B=|01>, C=|11> (the corner); |10> is the unused slot.
  PRE  |psi> = (A + B + C)/sqrt3      WEAK PROBE (3rd qubit = probe)      POST |phi> = (A + B - C)/sqrt3
  TEXTBOOK weak values: A = +1, B = +1, corner 11 = -1, unused 10 = 0.
usage: three_box_11.py sim|real [keyfile] [chip]"""
import os, sys, math, json
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import RYGate
import three_box as tb

BOXES = {"A_00": (0, 0), "B_01": (0, 1), "corner_11": (1, 1), "unused_10": (1, 0)}
ANG = 2 * math.asin(1 / math.sqrt(3))

def prep(qc, s1, s0):
    qc.ry(ANG, s1); qc.x(s1); qc.ch(s1, s0); qc.x(s1); qc.cx(s1, s0)      # (|00>+|01>+|11>)/sqrt3

def unprep(qc, s1, s0):
    qc.cx(s1, s0); qc.x(s1); qc.ch(s1, s0); qc.x(s1); qc.ry(-ANG, s1)

def circuit(box):
    qc = QuantumCircuit(3, 3); s0, s1, p = 0, 1, 2
    prep(qc, s1, s0)
    b1, b0 = BOXES[box]
    if b1 == 0: qc.x(s1)
    if b0 == 0: qc.x(s0)
    qc.append(RYGate(2 * tb.EPS).control(2), [s1, s0, p])
    if b1 == 0: qc.x(s1)
    if b0 == 0: qc.x(s0)
    qc.cz(s1, s0)                        # |phi> = CZ|psi> (flips only the 11 corner); undo: CZ then prep^dagger
    unprep(qc, s1, s0)
    qc.h(p); qc.measure([s0, s1, p], [0, 1, 2]); return qc

tb.BOXES = BOXES; tb.circuit = circuit; tb.TEXTBOOK = {"A_00": 1, "B_01": 1, "corner_11": -1, "unused_10": 0}; tb.OUTNAME = "three_box_11"
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sim"
    import io, contextlib
    tb.main(mode, sys.argv[2] if len(sys.argv) > 2 else "~/.ibm_quantum_key", sys.argv[3] if len(sys.argv) > 3 else "ibm_marrakesh")
