#!/usr/bin/env python3
"""verify_september.py -- recompute every number on heuremen.org/one-idea.html from the stored JSON in data/.
Standard library only (math, json, os). No network. Reads data/*.json; writes nothing.
Usage:  python3 verify_september.py
Status per claim:
  OK                |recomputed - claimed| <= half the last printed digit of the claim
  DIFFERS           recomputed from stored data, outside that tolerance
  NOT REPRODUCIBLE  the number cannot be derived from the files in data/ (reason given)
Recipes follow the scripts that produced the numbers (space_time_same.py, arrow_panel.py, three_box.py,
three_box_11.py, entanglement_cost.py, entanglement_cost_v2.py, partner_twist.py). Statistics (SE = sd(ddof=1)/sqrt(n) across pairs)
follow the analyze() functions in those scripts.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
def load(name): return json.load(open(os.path.join(DATA, name)))

COUNTS = {"OK": 0, "DIFFERS": 0, "NOT REPRODUCIBLE": 0}
def decimals(s):
    s = s.replace("−", "-").strip()
    return len(s.split(".")[1]) if "." in s else 0
def num(s): return float(s.replace("−", "-").replace("+", ""))

def check(claim, claimed, got, note=""):
    """claimed: string exactly as printed on the page."""
    d = decimals(claimed); tol = 0.5 * 10 ** (-d) + 1e-9
    ok = abs(got - num(claimed)) <= tol
    st = "OK" if ok else "DIFFERS"
    COUNTS[st] += 1
    print(f"[{st:<16}] {claim} | page {claimed} | recomputed {got:+.{d + 3}f} (tol {tol - 1e-9:.4g}){' | ' + note if note else ''}")
    return ok

def check_int(claim, claimed, got, note=""):
    ok = abs(got - claimed) <= 0.5
    st = "OK" if ok else "DIFFERS"; COUNTS[st] += 1
    print(f"[{st:<16}] {claim} | page {claimed} | recomputed {got:.4f} (integer rounding){' | ' + note if note else ''}")

def nr(claim, claimed, reason):
    COUNTS["NOT REPRODUCIBLE"] += 1
    print(f"[NOT REPRODUCIBLE] {claim} | page {claimed} | {reason}")

def mean(x): x = list(x); return sum(x) / len(x)
def sd(x):
    x = list(x); m = mean(x); return math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))
def se(x): return sd(x) / math.sqrt(len(list(x)))
def section(t): print(f"\n=== {t} ===")

THETAS = [k * math.pi / 4 for k in range(8)]

# ---------------------------------------------------------------- 1. SPACE AND TIME
section("1 Space and time (space_time_same_real.json; time curve = partner_twist_..._162216.json ENT arm)")
st = load("space_time_same_real.json")
dl2 = load("partner_twist_real_ibm_marrakesh_20260924_162216.json")
check_pairs = sorted(st["space"].keys())
print(f"    pairs in file: {check_pairs}; time_from: {st['time_from']}; backend of time file: {dl2['backend']}")
claim_pairs_ok = check_pairs == ["(1, 2)", "(133, 134)"]
COUNTS["OK" if claim_pairs_ok else "DIFFERS"] += 1
print(f"[{'OK' if claim_pairs_ok else 'DIFFERS':<16}] 'ibm_marrakesh, pairs 1-2 and 133-134' | file pairs {check_pairs}, time file backend {dl2['backend']}")
PAIRS = [(1, 2), (133, 134)]
space = lambda t: mean(dict((round(a, 6), v) for a, v in st["space"][str(p)])[round(t, 6)] for p in PAIRS)
def time_c(t):
    vals = []
    for p in PAIRS:
        r = [r for r in dl2["rows"] if (r["a"], r["b"]) == p and r["kind"] == "ENT" and abs(r["theta"] - t) < 1e-9]
        vals.append(r[0]["C"])
    return mean(vals)
table = [("0.000", "+1.000", "+0.929", "+0.930"), ("0.785", "+0.707", "+0.666", "+0.664"),
         ("1.571", "+0.000", "-0.008", "+0.016"), ("2.356", "-0.707", "-0.672", "-0.664"),
         ("3.142", "-1.000", "-0.937", "-0.936"), ("3.927", "-0.707", "-0.680", "-0.657"),
         ("4.712", "+0.000", "-0.003", "-0.004"), ("5.498", "+0.707", "+0.683", "+0.653")]
for t, (th, cs, sp, tm) in zip(THETAS, table):
    check(f"table theta column ({th})", th, t)
    check(f"table textbook cos at theta={th}", cs, math.cos(t))
    check(f"table SPACE at theta={th}", sp, space(t))
    check(f"table TIME  at theta={th}", tm, time_c(t))

section("1 Replications at pi (space_time_same_real_ibm_kingston / _fez; arrow_panel kingston)")
kg = load("space_time_same_real_ibm_kingston.json"); fz = load("space_time_same_real_ibm_fez.json")
def at_pi(d, p): return dict((round(a, 6), v) for a, v in d["space"][str(p)])[round(math.pi, 6)]
check("kingston space at pi (pair 144-145)", "-0.813", at_pi(kg, (144, 145)),
      "page uses ONLY pair 144-145; the other kingston pair 16-23 is excluded in the lab log (dead readout, space at pi = %+.3f); mean of both pairs would be %+.3f" % (at_pi(kg, (16, 23)), mean([at_pi(kg, (144, 145)), at_pi(kg, (16, 23))])))
check("fez space at pi (mean of pairs 22-23, 81-82)", "-0.891", mean([at_pi(fz, (22, 23)), at_pi(fz, (81, 82))]),
      "pair values %+.4f, %+.4f" % (at_pi(fz, (22, 23)), at_pi(fz, (81, 82))))
ak = load("arrow_panel_real_ibm_kingston.json")
def arrow_link(d, kind, t=math.pi):
    return mean(r["value"] for r in d["rows"] if r["kind"] == kind and r["theta"] is not None and abs(r["theta"] - t) < 1e-9)
check("kingston 'time' at pi = one-qubit two-measurement link (arrow_panel kingston LINK+ at pi, mean of qubits 1,2,133)", "-0.840",
      arrow_link(ak, "LINK+"),
      "source is arrow_panel_real_ibm_kingston.json, NOT the partner_twist file named in space_time_same 'time_from' (no kingston partner_twist run is in data/); not pair-specific to 144-145")

# ---------------------------------------------------------------- 2. THE ARROW
section("2 The arrow (arrow_panel_real_ibm_marrakesh.json, qubits 1,2,133, 2048 shots)")
am = load("arrow_panel_real_ibm_marrakesh.json"); af = load("arrow_panel_real_ibm_fez.json")
def km(d, kind): return mean(r["value"] for r in d["rows"] if r["kind"] == kind)
check("coherence: start", "0.98", km(am, "BASE"), "mean <X> over qubits 1,2,133")
check("coherence: twist and untwist", "0.98", km(am, "REV"), "mean over 3 qubits x 4 angles (pi/4..pi)")
check("coherence: after one mid-circuit measurement", "+0.02", km(am, "LOOK"), "mean over 3 qubits")
check("fez link forward at pi", "-0.925", arrow_link(af, "LINK+"), "mean qubits 1,2,133")
check("fez link backward at pi", "-0.919", arrow_link(af, "LINK-"), "mean qubits 1,2,133")
print(f"    shots field arrow_panel marrakesh = {am['shots']}; backend {am['backend']}; qubits {sorted({r['qubit'] for r in am['rows']})}")
print(f"    (marrakesh link fwd/bwd at pi for reference: {arrow_link(am,'LINK+'):+.4f} / {arrow_link(am,'LINK-'):+.4f}; the page does not quote these)")

# ---------------------------------------------------------------- 3. ONE QUBIT OF AN ENTANGLED PAIR
section("3 One qubit of an entangled pair, read alone (space_time_same bloch: length = sqrt(X^2+Y^2+Z^2) of qubit a)")
def blen(d, p, arm): b = d["bloch"][str(p)]; return math.sqrt(sum(b[arm + "_" + x] ** 2 for x in "XYZ"))
for p in ["(1, 2)", "(133, 134)"]:
    print(f"    marrakesh pair {p}: alone {blen(st,p,'ALONE'):.4f}  entangled {blen(st,p,'ENT'):.4f}")
check("marrakesh alone length (mean of 2 pairs)", "0.99", mean(blen(st, p, "ALONE") for p in ["(1, 2)", "(133, 134)"]))
check("marrakesh entangled length, pair (1,2)", "0.06", blen(st, "(1, 2)", "ENT"))
check("marrakesh entangled length, pair (133,134)", "0.05", blen(st, "(133, 134)", "ENT"))
check("kingston entangled (pair 144-145)", "0.061", blen(kg, "(144, 145)", "ENT"), "pair 16-23 excluded as dead readout (alone length %.3f, impossible >1)" % blen(kg, "(16, 23)", "ALONE"))
check("fez alone (pair 22-23)", "0.995", blen(fz, "(22, 23)", "ALONE"))
check("fez entangled (pair 22-23)", "0.064", blen(fz, "(22, 23)", "ENT"))
print(f"    fez other pair (81-82): alone {blen(fz,'(81, 82)','ALONE'):.4f} entangled {blen(fz,'(81, 82)','ENT'):.4f} (not quoted on page)")

# ---------------------------------------------------------------- 4. THREE BOX
section("4 Three-box weak values (three_box_real_*.json; weak = Xp/(2 eps); SE from kept + Xp as in three_box.py)")
def se_from(xp, n, eps): return math.sqrt((1 - xp * xp) * n / (n - 1)) / math.sqrt(n) / (2 * eps)
chips = {"marrakesh": load("three_box_real_ibm_marrakesh.json"), "kingston": load("three_box_real_ibm_kingston.json"), "fez": load("three_box_real_ibm_fez.json")}
claims = {"marrakesh": {"A": ("+0.58", "0.06"), "B": ("+0.58", "0.06"), "C": ("-0.69", "0.06"), "corner_11": ("+0.05", "0.06")},
          "kingston": {"A": ("+0.80", "0.06"), "B": ("+0.84", "0.05"), "C": ("-0.59", "0.06"), "corner_11": ("+0.17", "0.06")},
          "fez": {"A": ("+0.75", "0.06"), "B": ("+0.91", "0.05"), "C": ("-0.50", "0.06"), "corner_11": ("+0.09", "0.06")}}
sig = {}
for chip, d in chips.items():
    check(f"{chip}: eps", "0.25", d["eps"]); check_int(f"{chip}: shots (8192)", 8192, d["shots"])
    for box, (w, s) in claims[chip].items():
        h = d["hardware"][box]; eps = d["eps"]
        weak = h["Xp"] / (2 * eps); s_re = se_from(h["Xp"], h["kept"], eps)
        check(f"{chip} {box} weak value", w, weak, f"kept {h['kept']}, stored weak {h['weak']:+.4f}")
        check(f"{chip} {box} +/-", s, s_re, f"stored se {h['se']:.4f}")
        if box == "C": sig[chip] = abs(weak) / s_re
print("    C-box significance (|weak|/se): " + ", ".join(f"{k} {v:.1f}" for k, v in sig.items()))
check_int("intro: 'nine error bars below zero' (weakest chip C-box, fez)", 9, min(sig.values()), "marrakesh %.1f, kingston %.1f are larger" % (sig["marrakesh"], sig["kingston"]))

section("4 Second encoding: corner 11 (three_box_11_real_ibm_marrakesh.json)")
d11 = load("three_box_11_real_ibm_marrakesh.json"); h = d11["hardware"]; eps = d11["eps"]
W = {k: h[k]["Xp"] / (2 * eps) for k in h}; S = {k: se_from(h[k]["Xp"], h[k]["kept"], eps) for k in h}
check("corner 11 weak value", "-0.51", W["corner_11"]); check("corner 11 +/-", "0.06", S["corner_11"])
check("box 00 weak value", "+0.77", W["A_00"]); check("box 01 weak value", "+0.76", W["B_01"])
check("sum of 00 + 01 + 11", "1.03", W["A_00"] + W["B_01"] + W["corner_11"])
check("unused slot 10 weak value", "+0.17", W["unused_10"]); check("unused slot 10 +/-", "0.06", S["unused_10"])
check_int("'nine standard errors below zero' (corner 11 |weak|/se)", 9, abs(W["corner_11"]) / S["corner_11"])
check_int("11-corner run shots", 8192, d11["shots"])

# ---------------------------------------------------------------- FIVE-ROW TABLE
section("The five that came back textbook")

# --- row 1/2: entanglement_cost.py analyze(): D = (1-bell.gap) - (1 - mean(null gaps)) per pair, SE across pairs
def floor_D(chipdict, t):
    out = []
    for p in chipdict["pairs"]:
        R = [r for r in chipdict["rows"] if r["a"] == p["a"] and r["b"] == p["b"] and r["t_us"] == t]
        bell = [r for r in R if r["kind"] == "bell"]; nul = [r for r in R if r["kind"].startswith("null")]
        if bell and nul: out.append((1 - bell[0]["gap"]) - (1 - mean(x["gap"] for x in nul)))
    return out
f1 = load("entanglement_cost_real_20260924_131028.json")["chips"][0]; f2 = load("entanglement_cost_real_20260924_135201.json")["chips"][0]
D0 = floor_D(f1, 0)
check("row 1: floor D at t=0 (marrakesh, %d pairs)" % len(D0), "+0.0012", mean(D0))
check("row 1: floor D at t=0 +/-", "0.0024", se(D0))
D30m, D30k = floor_D(f1, 30), floor_D(f2, 30)
print(f"    t=30us marrakesh {mean(D30m):+.5f} +/- {se(D30m):.5f} (n={len(D30m)}); kingston {mean(D30k):+.5f} +/- {se(D30k):.5f} (n={len(D30k)})")
w = [1 / se(D30m) ** 2, 1 / se(D30k) ** 2]; comb = (w[0] * mean(D30m) + w[1] * mean(D30k)) / sum(w); comb_se = math.sqrt(1 / sum(w))
check("row 2: combined D at 30us (inverse-variance, 2 chips)", "+0.0076", comb)
check("row 2: combined D at 30us +/-", "0.0025", comb_se)
check_int("row 2: '(3 sigma)'", 3, comb / comb_se, "z = %.2f" % (comb / comb_se))

# --- row 3/4: entanglement_cost_v2.py analyze() reimplemented without numpy
def solve(A, b):
    n = len(A); M = [r[:] + [b[i]] for i, r in enumerate(A)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(M[r][i])); M[i], M[piv] = M[piv], M[i]
        for r in range(n):
            if r != i:
                f = M[r][i] / M[i][i]; M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [M[i][n] / M[i][i] for i in range(n)]
def invert(A):
    n = len(A); cols = [solve(A, [1.0 if i == j else 0.0 for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]
v2 = load("entanglement_cost_v2_real_ibm_marrakesh_20260924_145627.json")
per = []
for p in v2["pairs"]:
    k = (p["a"], p["b"]); R = [r for r in v2["rows"] if (r["a"], r["b"]) == k]
    bell = [r for r in R if r["kind"] == "bell"]; nulls = [r for r in R if r["kind"].startswith("null")]
    if not bell or len(nulls) < 2: continue
    null_def = 1 + mean(r["ZZ"] for r in nulls)
    th = [r["theta"] for r in bell]; y = [1 + r["ZZ"] for r in bell]
    X = [[1, math.cos(t), math.sin(t), math.cos(2 * t), math.sin(2 * t)] for t in th]
    XtX = [[sum(X[n][i] * X[n][j] for n in range(len(X))) for j in range(5)] for i in range(5)]
    Xty = [sum(X[n][i] * y[n] for n in range(len(X))) for i in range(5)]
    beta = solve(XtX, Xty)
    rss = sum((y[n] - sum(X[n][i] * beta[i] for i in range(5))) ** 2 for n in range(len(X)))
    sigma2 = rss / (len(th) - 5); inv = invert(XtX)
    se1 = math.sqrt(sigma2 * (inv[1][1] + inv[2][2]) / 2)
    per.append({"order": bell[0]["order"], "D": beta[0] - null_def, "a1": math.hypot(beta[1], beta[2]), "se1": se1})
f = [x["D"] for x in per if x["order"] == "bell_first"]; r_ = [x["D"] for x in per if x["order"] == "bell_last"]
check("row 3: D_forward (bell first)", "+0.0124", mean(f)); check("row 3: D_forward +/-", "0.0046", se(f))
check("row 3: D_reversed (bell last)", "-0.0013", mean(r_)); check("row 3: D_reversed +/-", "0.0154", se(r_))
chi = sum((x["a1"] / x["se1"]) ** 2 for x in per)
check("row 4: chi2 of A1 (phase sweep)", "19.4", chi, "dof = 2 x %d pairs = %d" % (len(per), 2 * len(per)))
check_int("row 4: dof", 20, 2 * len(per))

# --- row 5: partner_twist.py analyze(): FORWARD = mean over pairs of mean_theta(C_ENT - C_PROD)
def link_forward(d):
    fwd = []
    for p in d["pairs"]:
        k = (p["a"], p["b"]); R = {(r["kind"], round(r["theta"], 6)): r for r in d["rows"] if (r["a"], r["b"]) == k}
        if ("PROD0", 0.0) in R:
            for t in THETAS: R[("PROD", round(t, 6))] = {"C": (R[("PROD0", round(t, 6))]["C"] + R[("PROD1", round(t, 6))]["C"]) / 2}
        fwd.append(mean(R[("ENT", round(t, 6))]["C"] - R[("PROD", round(t, 6))]["C"] for t in THETAS))
    return fwd
dl1 = load("partner_twist_real_ibm_marrakesh_20260924_161242.json")
f1_, f2_ = link_forward(dl1), link_forward(dl2)
check("row 5: first run ENT-PROD (unbalanced), %d pairs" % len(f1_), "-0.0076", mean(f1_)); check("row 5: first run +/-", "0.0036", se(f1_))
check("row 5: balanced rerun ENT-PROD, %d pairs" % len(f2_), "-0.0014", mean(f2_)); check("row 5: balanced rerun +/-", "0.0026", se(f2_))
print("\nSUMMARY: " + ", ".join(f"{k} {v}" for k, v in COUNTS.items()))
