#!/usr/bin/env python3
"""entanglement_cost_v2.py -- 2026-09-24, follow-up designed to test an ordering / TLS-bath explanation of run 1.

Run 1 found: entangled circuits lose ~0.0076 more than balanced product nulls after a 30 us idle (z~3, 2 chips),
nothing at t=0. The leading ordinary explanation: bell ALWAYS ran first in each block, so it saw the coldest
TLS bath; the nulls excite one qubit every shot (TLS saturation -> longer effective T1). Coherence DOES survive 30 us
(XX=YY up to +/-0.7), and collective decay / exchange are dead (sign, secularity). So two questions:

  Q1 ORDER (signed, free): half the pairs run bell-FIRST, half bell-LAST. Ordering/TLS predicts D flips sign or
     vanishes in the bell-last half. An effect of the state itself doesn't care about order.
  Q2 PHASE: bell with RZ(theta) on b BETWEEN the CX and the RY(pi) -> (|01> - e^{i theta}|10>)/sqrt2. The RZ is
     virtual but rides on a real pi pulse, so every theta has identical pulses, schedule, populations and per-shot
     excitation statistics. Any theta-dependence of the loss is coherence-driven -- which no textbook channel allows
     at these detunings (secular approximation). The theta values run in a shuffled binding order.

Observables: ZZ (gap), XX, YY, XY, YX (full single-excitation coherence, not a lower bound).
usage: entanglement_cost_v2.py sim|real [pairs] [chip]   (real needs ~/.ibm_quantum_key or DF_KEY=path)
       entanglement_cost_v2.py analyze <json>
"""
import json, sys, os, math, random
from datetime import datetime, timezone
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import Parameter
from qiskit.quantum_info import SparsePauliOp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entanglement_cost import pair_table, pick_pairs      # same pair selection + impossible-calibration filter as run 1

T_US = 30
N_THETA = 12
SHOTS_THETA = 2048   # per binding: pairs beat shots (per-pair D is shot-noise-limited anyway)
SHOTS_NULL = 4096
OUT_DIR = os.path.expanduser("./runs")
OBS = ["ZZ", "XX", "YY", "XY", "YX"]


def circuits():
    th = Parameter("theta")
    bell = QuantumCircuit(2); bell.h(0); bell.cx(0, 1); bell.rz(th, 1); bell.ry(math.pi, 1)
    nb = QuantumCircuit(2); nb.h(0); nb.h(0); nb.cx(0, 1); nb.ry(math.pi, 1)
    na = QuantumCircuit(2); na.x(0); na.cx(0, 1); na.ry(math.pi, 1)
    for qc in (bell, nb, na):
        qc.delay(T_US, 0, unit="us"); qc.delay(T_US, 1, unit="us")
    return bell, nb, na


def build(backend, pairs, seed=212):
    rng = random.Random(seed)
    bell, nb, na = circuits()
    thetas = [2 * math.pi * k / N_THETA for k in range(N_THETA)]
    pubs, index = [], []
    for i, p in enumerate(pairs):
        order = "bell_first" if i % 2 == 0 else "bell_last"
        lay = [p["a"], p["b"]]
        tb = transpile(bell, backend=backend, initial_layout=lay, optimization_level=0)
        tn = [transpile(c, backend=backend, initial_layout=lay, optimization_level=0) for c in (nb, na)]
        obs_b = [SparsePauliOp(o).apply_layout(tb.layout) for o in OBS]
        th_order = thetas[:]; rng.shuffle(th_order)
        bell_pub = (tb, [[o] for o in obs_b], np.array(th_order).reshape(-1, 1), 1 / math.sqrt(SHOTS_THETA))   # 4th slot is PRECISION (1/sqrt(shots)), not shots -- rehearsal bug 9/24
        null_pubs = [(t, [SparsePauliOp(o).apply_layout(t.layout) for o in OBS], None, 1 / math.sqrt(SHOTS_NULL)) for t in tn]
        seq = [("bell", bell_pub)] + list(zip(("null_b", "null_a"), null_pubs))
        if order == "bell_last":
            seq = seq[1:] + seq[:1]
        for kind, pub in seq:
            pubs.append(pub)
            index.append({"a": p["a"], "b": p["b"], "kind": kind, "order": order,
                          "thetas": th_order if kind == "bell" else None})
    return pubs, index


def run(mode, n_pairs, chip):
    from qiskit_ibm_runtime import EstimatorV2
    if mode == "sim":
        from qiskit_ibm_runtime.fake_provider import FakeKingston, FakeMarrakesh
        be = {"ibm_marrakesh": FakeMarrakesh, "ibm_kingston": FakeKingston}.get(chip, FakeKingston)()
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        key = open(os.environ.get("DF_KEY") or os.path.expanduser("~/.ibm_quantum_key")).read().strip()
        be = QiskitRuntimeService(channel="ibm_quantum_platform", token=key, instance=os.environ.get("DF_INSTANCE") or None).backend(chip)
    pairs = pick_pairs(pair_table(be), n_pairs)
    pubs, index = build(be, pairs)
    shots = sum(round((len(p[2]) if p[2] is not None else 1) / p[3] ** 2) for p in pubs)
    print(f"[{be.name}] {len(pairs)} pairs, {len(pubs)} pubs, {shots} shots", flush=True)
    # 9/24: one job for all 30 pubs failed with IBM error 6073 (exceeds classical control memory: 12 thetas x 32 readout
    # randomizations x 10 pairs). Now ONE JOB PER PAIR (its 3 pubs, in that pair's pre-registered order), all submitted
    # inside a Batch so they queue together.
    groups = [list(range(i, i + 3)) for i in range(0, len(pubs), 3)]
    def make_est(m):
        est = EstimatorV2(mode=m)
        if mode == "real":
            est.options.resilience_level = 1
            est.options.dynamical_decoupling.enable = False
        return est
    jobs = []
    if mode == "real":
        from qiskit_ibm_runtime import Batch
        with Batch(backend=be) as batch:
            est = make_est(batch)
            for g in groups:
                jobs.append(est.run([pubs[i] for i in g]))
    else:
        est = make_est(be)
        jobs = [est.run([pubs[i] for i in g]) for g in groups]
    jid = [j.job_id() if hasattr(j, "job_id") else None for j in jobs]
    print(f"[{be.name}] {len(jobs)} jobs: {jid}", flush=True)
    res = [r for j in jobs for r in j.result()]
    rows = []
    for meta, r in zip(index, res):
        ev = np.asarray(r.data.evs)
        if meta["kind"] == "bell":
            ev = ev.reshape(len(OBS), -1)               # obs x theta
            for j, th in enumerate(meta["thetas"]):
                rows.append({**{k: meta[k] for k in ("a", "b", "kind", "order")}, "theta": th,
                             **{o: float(ev[i, j]) for i, o in enumerate(OBS)}})
        else:
            rows.append({**{k: meta[k] for k in ("a", "b", "kind", "order")}, "theta": None,
                         **{o: float(v) for o, v in zip(OBS, ev.ravel())}})
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, f"entanglement_cost_v2_{mode}_{be.name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
    json.dump({"mode": mode, "backend": be.name, "job_id": jid, "t_us": T_US, "pairs": pairs, "rows": rows,
               "started": datetime.now(timezone.utc).isoformat()}, open(path, "w"), indent=1)
    print("saved", path, flush=True)
    analyze(path)


def analyze(path):
    """Pre-registered statistics (9/24):
    Q1  D_invariant = (D_forward + D_reversed)/2  <- the ONLY physics number;  D_ordering = (D_forward - D_reversed)/2 reported always.
    Q2  per pair: deficit(theta) = D0 + A1 cos(theta+phi1) + A2 cos(2 theta+phi2). A1 is the test statistic. phi1 is a NUISANCE
        (absorbs accumulated drive-detuning phase), never interpreted. A2 must be ~0: a nonzero A2 is a FALSIFIER (binding systematic),
        not a bonus. Within-pair fit: pair-to-pair variance cancels out of A1."""
    d = json.load(open(path))
    rows = d["rows"]
    print(f"\n=== entanglement_cost v2 ({d['mode']}, {d['backend']}, t={d['t_us']}us) ===")
    per = []
    for p in d["pairs"]:
        k = (p["a"], p["b"])
        R = [r for r in rows if (r["a"], r["b"]) == k]
        bell = [r for r in R if r["kind"] == "bell"]
        nulls = [r for r in R if r["kind"].startswith("null")]
        if not bell or len(nulls) < 2: continue
        null_def = 1 + sum(r["ZZ"] for r in nulls) / len(nulls)          # deficit = 1 - gap = 1 + ZZ
        th = np.array([r["theta"] for r in bell]); defi = np.array([1 + r["ZZ"] for r in bell])
        X = np.vstack([np.ones_like(th), np.cos(th), np.sin(th), np.cos(2 * th), np.sin(2 * th)]).T
        beta, res, *_ = np.linalg.lstsq(X, defi, rcond=None)
        sigma = math.sqrt(res[0] / (len(th) - 5)) if len(res) and len(th) > 5 else float("nan")
        cov = sigma ** 2 * np.linalg.inv(X.T @ X)
        a1 = math.hypot(beta[1], beta[2]); a2 = math.hypot(beta[3], beta[4])
        coh = np.mean([math.hypot((r["XX"] + r["YY"]) / 2, (r["XY"] - r["YX"]) / 2) for r in bell])
        per.append({"k": k, "order": bell[0]["order"], "D": beta[0] - null_def, "c1": beta[1], "s1": beta[2],
                    "a1": a1, "a2": a2, "se1": math.sqrt((cov[1, 1] + cov[2, 2]) / 2), "coh": coh})
        print(f"  {k[0]:>3}-{k[1]:<3} {bell[0]['order']:10s} D {beta[0] - null_def:+.4f}  A1 {a1:.4f}  A2 {a2:.4f}  (fit sd per point {sigma:.4f})  |coh| {coh:.2f}")
    f = [x["D"] for x in per if x["order"] == "bell_first"]; r = [x["D"] for x in per if x["order"] == "bell_last"]
    if len(f) >= 2 and len(r) >= 2:
        mf, mr = np.mean(f), np.mean(r)
        sf, sr = np.std(f, ddof=1) / math.sqrt(len(f)), np.std(r, ddof=1) / math.sqrt(len(r))
        se = 0.5 * math.hypot(sf, sr)
        print(f"  Q1 D_forward {mf:+.5f}+/-{sf:.5f}  D_reversed {mr:+.5f}+/-{sr:.5f}")
        print(f"  Q1 D_invariant = {(mf + mr) / 2:+.5f} +/- {se:.5f}  (95% [{(mf + mr) / 2 - 1.96 * se:+.5f}, {(mf + mr) / 2 + 1.96 * se:+.5f}])   <- headline")
        print(f"  Q1 D_ordering  = {(mf - mr) / 2:+.5f} +/- {se:.5f}")
    if per:
        # A1 is positive-biased by noise; the unbiased pooled test is on the (c1, s1) components: under "no modulation"
        # sum over pairs of (a1/se1)^2 ~ chi2 with 2*npairs dof.
        chi = sum((x["a1"] / x["se1"]) ** 2 for x in per if x["se1"] > 0); dof = 2 * len(per)
        print(f"  Q2 A1: median {np.median([x['a1'] for x in per]):.4f};  chi2 = {chi:.1f} on {dof} dof (no-modulation expectation ~{dof})")
        print(f"  Q2 A2 (falsifier, must be ~noise): median {np.median([x['a2'] for x in per]):.4f}")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "analyze":
        analyze(sys.argv[2])
    else:
        run(sys.argv[1] if len(sys.argv) > 1 else "sim", int(sys.argv[2]) if len(sys.argv) > 2 else 10,
            sys.argv[3] if len(sys.argv) > 3 else "ibm_kingston")
