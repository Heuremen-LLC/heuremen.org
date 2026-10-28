#!/usr/bin/env python3
"""entanglement_cost.py -- 2026-09-24. Does a Bell pair lose more fidelity than an unentangled product-state control?

BACKGROUND
  A precision run (2026-08-13, ibm_marrakesh, 200k shots, readout-mitigated) measured the Bell-pair endpoint gap at
  0.98983 +/- 0.00034 where perfect is 1.0. This script asks whether the 1.02% shortfall is ordinary hardware noise
  or specific to the entangled state. The reference size 1/pi^4 = 0.010266 is used only as a comparison scale for the
  bell-minus-null difference.

HOW HARDWARE NOISE IS SEPARATED FROM AN ENTANGLEMENT-SPECIFIC LOSS (decided before any real run)
  * Across qubit pairs: noise TRACKS each pair's two-qubit (CZ) gate error; a pair-independent floor would not.
  * Against a null: the PRODUCT state |0>|1> reads the same ideal gap (1.0) with the same gates and no entanglement.
    Noise hits both about the same. An entanglement-specific loss would show in the Bell pair only.
  * Time: an idle delay after the pair is made. Ordinary decay is PREDICTED from
    each qubit's calibrated T1 (gap loses p0+p1 for the pair, 2*p1 for the null, p = 1-exp(-t/T1)). What is compared is
    whatever is left after that prediction, and whether it changes with time.

WHAT IS MEASURED (EstimatorV2, readout mitigation resilience_level=1 = same as the August run)
  gap = -<ZZ>        (the endpoint quantity; 1.0 ideal for BOTH circuits -- not entanglement-specific)
  coh = -(<XX>+<YY>)/2  (1.0 ideal for the entangled pair, 0 for the product null -- only entanglement lights it)

CIRCUITS (physical pair pinned by initial_layout; optimization_level=0 so no gates are optimized away)
  bell: H(a) CX(a,b) RY(pi,b) DELAY(t) measure    -> (|01> - |10>)/sqrt2 : one excitation per branch during delay
  null: H(a) H(a) CX(a,b) RY(pi,b) DELAY(t)       -> |0>|1>              : one excitation during delay, same CX

usage:  entanglement_cost.py sim  [pairs_per_chip]     local dress rehearsal on IBM's fake-chip noise models (free)
        entanglement_cost.py real [pairs_per_chip]     real hardware (needs ~/.ibm_quantum_key)
        entanglement_cost.py analyze <results.json>"""
import json, sys, os, math, random
from datetime import datetime, timezone
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import SparsePauliOp

DELAYS_US = [int(x) for x in os.environ["DF_DELAYS"].split(",")] if os.environ.get("DF_DELAYS") else [0, 30]   # two points: now, and after time passes (QPU budget; T1 ~100-300us makes 30us a clear, predicted decay)
SHOTS = 8192   # per pub; rehearsal showed bell-minus-null sd ~0.013 at 4096 -- need tighter to see ~0.01
OUT_DIR = os.path.expanduser("./runs")


KINDS = ("bell", "null_b", "null_a")

def circuits(t_us):
    """bell  : (|01>-|10>)/sqrt2 -- the excitation is shared by both qubits
       null_b: |0>|1>  (b excited)  -- same CX, no entanglement
       null_a: |1>|0>  (a excited)  -- same CX, no entanglement
       The two nulls AVERAGED carry exactly the bell pair's excitation load (p_a + p_b), so unequal T1s can't fake a
       bell-minus-null difference. All three read gap = 1.0 ideally."""
    bell = QuantumCircuit(2); bell.h(0); bell.cx(0, 1); bell.ry(math.pi, 1)
    nb = QuantumCircuit(2); nb.h(0); nb.h(0); nb.cx(0, 1); nb.ry(math.pi, 1)
    na = QuantumCircuit(2); na.x(0); na.cx(0, 1); na.ry(math.pi, 1)
    for qc in (bell, nb, na):
        if t_us > 0:
            qc.delay(t_us, 0, unit="us"); qc.delay(t_us, 1, unit="us")
    return bell, nb, na


def pair_table(backend):
    """Every coupled pair with CZ error, both qubits' T1 and readout error, from the backend's own calibration."""
    t = backend.target
    props = {}
    for q in range(backend.num_qubits):
        qp = t.qubit_properties[q] if t.qubit_properties else None
        ro = t["measure"][(q,)].error if ("measure" in t and (q,) in t["measure"] and t["measure"][(q,)]) else None
        props[q] = {"T1_us": (qp.t1 * 1e6) if (qp and qp.t1) else None, "T2_us": (qp.t2 * 1e6) if (qp and qp.t2) else None, "ro": ro}
    rows = []
    for (a, b), ip in t["cz"].items():
        if a > b or ip is None or ip.error is None: continue
        pa, pb = props[a], props[b]
        if None in (pa["T1_us"], pb["T1_us"], pa["ro"], pb["ro"]): continue
        if pa["ro"] > 0.05 or pb["ro"] > 0.05 or pa["T1_us"] < 40 or pb["T1_us"] < 40: continue   # dead/dying qubits out
        if any(q["T2_us"] is None or q["T2_us"] > 2 * q["T1_us"] for q in (pa, pb)): continue       # impossible calibration (T2 > 2*T1) = untrustworthy numbers, skip
        if ip.error >= 0.5: continue
        rows.append({"a": a, "b": b, "cz_err": ip.error, "T1a": pa["T1_us"], "T1b": pb["T1_us"], "roa": pa["ro"], "rob": pb["ro"]})
    return rows


def pick_pairs(rows, n):
    """Spread by VALUE across the working CZ-error range (<= 2%: a pair worse than that is broken, not 'high noise',
    and one broken pair dominates the slope fit -- rehearsal lesson 9/24). Evenly spaced targets from best to 2%-cap
    worst, nearest pair to each target, no shared qubits."""
    rows = sorted([r for r in rows if r["cz_err"] <= 0.02], key=lambda r: r["cz_err"])
    if not rows: return []
    lo, hi = rows[0]["cz_err"], rows[-1]["cz_err"]
    ranks = [min(range(len(rows)), key=lambda i: abs(rows[i]["cz_err"] - (lo + k * (hi - lo) / max(1, n - 1)))) for k in range(n)]
    chosen, used = [], set()
    for k in range(n):
        target = ranks[k]
        for off in range(len(rows)):
            for i in (target - off, target + off):
                if 0 <= i < len(rows) and rows[i]["a"] not in used and rows[i]["b"] not in used:
                    chosen.append(rows[i]); used |= {rows[i]["a"], rows[i]["b"]}; break
            else:
                continue
            break
    return chosen


def build_pubs(backend, pairs):
    obs = SparsePauliOp.from_list([("ZZ", 1)]), SparsePauliOp.from_list([("XX", 1)]), SparsePauliOp.from_list([("YY", 1)])
    pubs, index = [], []
    for p in pairs:
        for t in DELAYS_US:
            for kind, qc in zip(KINDS, circuits(t)):
                # no scheduling_method: ALAP pads EVERY idle qubit on the chip with delays (incl. broken ones);
                # the backend schedules its own timing. Our explicit delays stay on our two qubits only.
                tqc = transpile(qc, backend=backend, initial_layout=[p["a"], p["b"]], optimization_level=0)
                pubs.append((tqc, [o.apply_layout(tqc.layout) for o in obs]))
                index.append({"a": p["a"], "b": p["b"], "t_us": t, "kind": kind})
    return pubs, index


def run(mode, n_pairs):
    from qiskit_ibm_runtime import EstimatorV2
    if mode == "sim":
        from qiskit_ibm_runtime.fake_provider import FakeKingston, FakeMarrakesh
        backends = [FakeKingston(), FakeMarrakesh()]
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        key = open(os.path.expanduser("~/.ibm_quantum_key")).read().strip()
        svc = QiskitRuntimeService(channel="ibm_quantum_platform", token=key)
        names = [b.name for b in svc.backends(operational=True)]
        req = [x for x in (os.environ.get("DF_CHIPS") or "").split(",") if x]
        want = [n for n in (req or ["ibm_kingston", "ibm_marrakesh"]) if n in names] or names[:2]
        backends = [svc.backend(n) for n in want]
    os.makedirs(OUT_DIR, exist_ok=True)
    record = {"mode": mode, "started": datetime.now(timezone.utc).isoformat(), "delays_us": DELAYS_US, "shots": SHOTS, "chips": []}
    for be in backends:
        pairs = pick_pairs(pair_table(be), n_pairs)
        pubs, index = build_pubs(be, pairs)
        est = EstimatorV2(mode=be)
        est.options.default_shots = SHOTS
        if mode == "real":
            est.options.resilience_level = 1
            est.options.dynamical_decoupling.enable = False   # delays must be real idle time, not DD-protected
        print(f"[{be.name}] {len(pairs)} pairs x {len(DELAYS_US)} delays x {len(KINDS)} kinds = {len(pubs)} pubs", flush=True)
        job = est.run(pubs)
        jid = job.job_id() if hasattr(job, "job_id") else None
        print(f"[{be.name}] job {jid}", flush=True)
        res = job.result()
        rows = []
        for meta, r in zip(index, res):
            zz, xx, yy = (float(v) for v in r.data.evs)
            sd = [float(v) for v in r.data.stds] if hasattr(r.data, "stds") else [None] * 3
            rows.append({**meta, "zz": zz, "xx": xx, "yy": yy, "gap": -zz, "coh": -(xx + yy) / 2, "gap_std": sd[0]})
        record["chips"].append({"backend": be.name, "job_id": jid, "pairs": pairs, "rows": rows})
    path = os.path.join(OUT_DIR, f"entanglement_cost_{mode}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
    json.dump(record, open(path, "w"), indent=1)
    print("saved", path, flush=True)
    analyze(path)


def fit_line(x, y):
    x, y = np.asarray(x), np.asarray(y)
    A = np.vstack([x, np.ones_like(x)]).T
    (m, c), *_ = np.linalg.lstsq(A, y, rcond=None)
    # bootstrap over pairs
    rng = np.random.default_rng(212)
    bs = []
    for _ in range(2000):
        i = rng.integers(0, len(x), len(x))
        if len(set(x[i])) < 2: continue
        (bm, bc), *_ = np.linalg.lstsq(A[i], y[i], rcond=None); bs.append((bm, bc))
    bs = np.array(bs)
    return m, c, np.percentile(bs[:, 0], [2.5, 97.5]), np.percentile(bs[:, 1], [2.5, 97.5])


def analyze(path):
    rec = json.load(open(path))
    print(f"\n=== entanglement_cost analysis ({rec['mode']}) ===")
    for chip in rec["chips"]:
        P = {(p["a"], p["b"]): p for p in chip["pairs"]}
        print(f"\n[{chip['backend']}] job {chip['job_id']}")
        print("  pair     cz_err   t_us  bell_def  null_def  bell-null  T1pred_b  T1pred_n  bell_resid null_resid  coh")
        d0 = []
        for r in sorted(chip["rows"], key=lambda r: (P[(r['a'], r['b'])]['cz_err'], r["t_us"], r["kind"])):
            if r["kind"] != "bell": continue
            p = P[(r["a"], r["b"])]
            nn = [x for x in chip["rows"] if x["a"] == r["a"] and x["b"] == r["b"] and x["t_us"] == r["t_us"] and x["kind"].startswith("null")]
            t = r["t_us"]; pa = 1 - math.exp(-t / p["T1a"]); pb = 1 - math.exp(-t / p["T1b"])
            tb = tn = pa + pb          # bell shares one excitation; the two nulls averaged carry the same load
            bd, nd = 1 - r["gap"], 1 - sum(x["gap"] for x in nn) / len(nn)
            print(f"  {r['a']:>3}-{r['b']:<3} {p['cz_err']:.4f}  {t:>4}  {bd:8.4f}  {nd:8.4f}  {bd-nd:+8.4f}  {tb:8.4f}  {tn:8.4f}  {bd-tb:+9.4f} {nd-tn:+9.4f}  {r['coh']:.3f}")
            if t == min(rec["delays_us"]): d0.append((p["cz_err"], bd, nd))
        if len(d0) >= 3:
            x = [d[0] for d in d0]
            for label, idx in (("bell", 1), ("null", 2)):
                m, c, mci, cci = fit_line(x, [d[idx] for d in d0])
                print(f"  t={min(rec['delays_us'])} {label} deficit vs CZ error: slope {m:+.3f} [{mci[0]:+.3f},{mci[1]:+.3f}]  intercept {c:+.5f} [{cci[0]:+.5f},{cci[1]:+.5f}]")
            diff = [d[1] - d[2] for d in d0]
            print(f"  t={min(rec['delays_us'])} bell-minus-null deficit: mean {np.mean(diff):+.5f}  (sd {np.std(diff, ddof=1):.5f}, n={len(diff)})")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "analyze":
        analyze(sys.argv[2])
    else:
        run(sys.argv[1] if len(sys.argv) > 1 else "sim", int(sys.argv[2]) if len(sys.argv) > 2 else 6)
