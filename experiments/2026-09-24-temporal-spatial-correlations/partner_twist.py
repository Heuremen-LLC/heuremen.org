#!/usr/bin/env python3
"""partner_twist.py -- 2026-09-24. Does a qubit's two-time correlation depend on its entanglement with a partner (FORWARD),
and does the partner's statistics depend on the rotation applied to a between measurements (INVERSE)? The textbook says neither.

CIRCUIT (qubits a, b; one pinned pair; mid-circuit measurement = dynamic circuit on IBM Heron)
  prep   : ENT  -> H(a) CX(a,b)                  (a entangled with b)
           PROD -> H(a)                          (same marginal for a: 50/50; b untouched in |0>)
  meas 1 : measure a -> m1
  twist  : RX(theta) on a                        (the interval between measurements; theta swept)
  meas 2 : measure a -> m2 ; measure b -> mb    (b is never touched after prep)

TEXTBOOK PREDICTIONS (exact, for projective measurements):
  FORWARD  space -> time : C_a(theta) = <z1 z2> = cos(theta) for BOTH ENT and PROD. The first measurement collapses a; whether a
                           was entangled cannot change its two-time correlation. ENT - PROD = 0 at every theta.
  INVERSE  time -> space : P(mb=1) cannot depend on theta (no-signaling: what is done to a after measurement 1 cannot reach b's own
                           statistics). For ENT, P(mb=1) = 1/2 at every theta; for PROD, 0.
A LINK = ENT-PROD difference in C_a that isn't zero, or a theta-dependence of P(mb) in ENT.

Systematics designed out: circuit ORDER randomized (lesson of entanglement_cost), ENT and PROD share a's marginal, so asymmetric
readout errors hit both equally; PROD's b sits in |0> (not driven), and the ENT-PROD readout crosstalk difference is checked
at theta=0 and theta=pi, where the textbook values are +/-1 for both.
usage: partner_twist.py sim|real [pairs] [chip]      analyze <json>"""
import json, sys, os, math, random
from datetime import datetime, timezone
import numpy as np
from qiskit import QuantumCircuit, transpile, ClassicalRegister, QuantumRegister
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entanglement_cost import pair_table, pick_pairs

THETAS = [k * math.pi / 4 for k in range(8)]     # 0..7pi/4
SHOTS = 4096
OUT_DIR = os.path.expanduser("./runs")


def circuit(kind, theta):
    q = QuantumRegister(2, "q"); m1 = ClassicalRegister(1, "m1"); m2 = ClassicalRegister(1, "m2"); mb = ClassicalRegister(1, "mb")
    qc = QuantumCircuit(q, m1, m2, mb)
    # v2 BALANCED (9/24 ~20:25Z): every arm has exactly one CX and b is excited in half the shots on average.
    #   ENT  : H(a) CX(a->b)          -> (|00>+|11>)/sqrt2, b excited half the time, entangled
    #   PROD0: H(a) CX(b->a), b=|0>   -> no-op CX, a=|+>, b at rest,   not entangled
    #   PROD1: X(b) H(a) CX(b->a)     -> X on a leaves |+> unchanged, b excited, not entangled
    #   PROD = mean(PROD0, PROD1): same gate count, same length, same b-excitation pattern as ENT.
    if kind == "PROD1": qc.x(1)
    qc.h(0)
    if kind == "ENT": qc.cx(0, 1)
    elif kind in ("PROD0", "PROD1"): qc.cx(1, 0)
    qc.measure(0, m1[0])
    qc.rx(theta, 0)
    qc.measure(0, m2[0]); qc.measure(1, mb[0])
    return qc


def run(mode, n_pairs, chip):
    from qiskit_ibm_runtime import SamplerV2
    if mode == "sim":
        from qiskit_ibm_runtime.fake_provider import FakeKingston, FakeMarrakesh
        be = {"ibm_marrakesh": FakeMarrakesh, "ibm_kingston": FakeKingston}.get(chip, FakeKingston)()
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        key = open(os.environ.get("DF_KEY") or os.path.expanduser("~/.ibm_quantum_key")).read().strip()
        be = QiskitRuntimeService(channel="ibm_quantum_platform", token=key, instance=os.environ.get("DF_INSTANCE") or None).backend(chip)
    pairs = pick_pairs(pair_table(be), n_pairs)
    items = [(p, kind, th) for p in pairs for kind in ("ENT", "PROD0", "PROD1") for th in THETAS]
    random.Random(212).shuffle(items)                              # ORDER RANDOMIZED
    pubs = [(transpile(circuit(kind, th), backend=be, initial_layout=[p["a"], p["b"]], optimization_level=0),) for p, kind, th in items]
    print(f"[{be.name}] {len(pairs)} pairs x 3 kinds x {len(THETAS)} thetas = {len(pubs)} circuits, randomized order", flush=True)
    # chunks of 16 in a Batch (v2 lesson: one oversized job -> IBM 6073). Chunks keep the randomized order.
    chunks = [pubs[i:i + 16] for i in range(0, len(pubs), 16)]
    if mode == "real":
        from qiskit_ibm_runtime import Batch
        with Batch(backend=be) as batch:
            sampler = SamplerV2(mode=batch); sampler.options.default_shots = SHOTS
            jobs = [sampler.run(c) for c in chunks]
    else:
        sampler = SamplerV2(mode=be); sampler.options.default_shots = SHOTS
        jobs = [sampler.run(c) for c in chunks]
    print(f"[{be.name}] {len(jobs)} jobs: {[j.job_id() if hasattr(j, 'job_id') else None for j in jobs]}", flush=True)
    res = [r for j in jobs for r in j.result()]
    rows = []
    for (p, kind, th), r in zip(items, res):
        b1 = np.array(r.data.m1.to_bool_array()).ravel(); b2 = np.array(r.data.m2.to_bool_array()).ravel(); bb = np.array(r.data.mb.to_bool_array()).ravel()
        z1, z2 = 1 - 2 * b1.astype(int), 1 - 2 * b2.astype(int)
        rows.append({"a": p["a"], "b": p["b"], "kind": kind, "theta": th, "C": float(np.mean(z1 * z2)), "Pb": float(np.mean(bb)),
                     "P1": float(np.mean(b1)), "n": int(len(b1))})
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, f"partner_twist_{mode}_{be.name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
    json.dump({"mode": mode, "backend": be.name, "pairs": pairs, "rows": rows}, open(path, "w"), indent=1)
    print("saved", path, flush=True); analyze(path)


def analyze(path):
    d = json.load(open(path)); rows = d["rows"]
    print(f"\n=== partner_twist ({d['mode']}, {d['backend']}) ===")
    fwd, inv, ctl = [], [], []
    for p in d["pairs"]:
        k = (p["a"], p["b"]); R = {(r["kind"], round(r["theta"], 6)): r for r in rows if (r["a"], r["b"]) == k}
        if ("PROD0", 0.0) in R:   # balanced v2: PROD = mean of the two product arms
            for t in THETAS:
                a0, a1 = R[("PROD0", round(t, 6))], R[("PROD1", round(t, 6))]
                R[("PROD", round(t, 6))] = {"C": (a0["C"] + a1["C"]) / 2, "Pb": (a0["Pb"] + a1["Pb"]) / 2}
        dif = [R[("ENT", round(t, 6))]["C"] - R[("PROD", round(t, 6))]["C"] for t in THETAS]
        th = np.array(THETAS); X = np.vstack([np.ones_like(th), np.cos(th), np.sin(th)]).T
        swing = {}
        for kind in ("ENT", "PROD"):
            pb = np.array([R[(kind, round(t, 6))]["Pb"] for t in THETAS])
            (c0, c1, s1), *_ = np.linalg.lstsq(X, pb, rcond=None); swing[kind] = (math.hypot(c1, s1), c0)
        fwd.append(np.mean(dif)); inv.append(swing["ENT"][0]); ctl.append(swing["PROD"][0])
        print(f"  {k[0]:>3}-{k[1]:<3} FORWARD mean(C_ENT - C_PROD) {np.mean(dif):+.4f}   INVERSE P(mb) theta-swing ENT {swing['ENT'][0]:.4f} vs PROD {swing['PROD'][0]:.4f} (control)")
    se_f = np.std(fwd, ddof=1) / math.sqrt(len(fwd)) if len(fwd) > 1 else float("nan")
    print(f"  FORWARD  space->time : ENT-PROD = {np.mean(fwd):+.5f} +/- {se_f:.5f}   (textbook: 0)")
    print(f"  INVERSE  time->space : P(mb) theta-swing ENT median {np.median(inv):.4f}  vs PROD-control median {np.median(ctl):.4f}")
    print(f"           (textbook: both are noise. PROD's b sits near |0> so its binomial floor is SMALLER -- compare ENT to the")
    print(f"            ENT floor 0.5/sqrt(shots)*sqrt(2/n_theta) = {0.5 / math.sqrt(SHOTS) * math.sqrt(2 / len(THETAS)):.4f} per component)")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "analyze": analyze(sys.argv[2])
    else: run(sys.argv[1] if len(sys.argv) > 1 else "sim", int(sys.argv[2]) if len(sys.argv) > 2 else 6, sys.argv[3] if len(sys.argv) > 3 else "ibm_kingston")
