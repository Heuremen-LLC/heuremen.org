"""Derive the per-measured-qubit flip-rule failure on 21 Oct calibration-1 data (the numbers quoted in PLAN_FLIP.md, flip-readable guard).
Heuremen LLC, Wayfinder and Bones, 22 October 2026. Read-only on control_run/iqjob_<chip>.json + iqraw_<chip>.npz.
DATA: every CAL state '1' pub (both copies; register c, 156 columns), 2 x 2000 = 4000 shots per qubit; CAL state '0' likewise for the 0 cluster.
INDEXING: column = PHYSICAL qubit number (the calibration circuit measures qubit q into clbit q); measured qubits = rec['E'] (physical numbers).
CLUSTERS (as analyze_iq.py): mu0, mu1 = mean of ALL cal-0 / cal-1 shots; t = projection on mu1 - mu0; cov0 from cal-0 shots with t < 0.5,
cov1 from cal-1 shots with t > 0.5; 'inside' = Mahalanobis distance < 3.5.
RULES on each prepared-1 shot:
  PLAN  (the rule PLAN_FLIP.md quoted, and analyze_flip.py's flip rule): FAIL = inside the 0 cluster OR not inside the 1 cluster
  A     not inside the 1 cluster
  B     inside the 0 cluster (whether or not also inside 1)
  BOTH  inside both clusters
  P     reads 0 by projection (t < 0.5)"""
import json, os, numpy as np
CTRL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'control_run'); D2 = 3.5 ** 2
def cov2(x): return np.cov(np.stack([x.real, x.imag])) + 1e-12 * np.eye(2)
def m2(x, mu, cv):
    d = np.stack([x.real - mu.real, x.imag - mu.imag], -1); return np.einsum('...i,ij,...j->...', d, np.linalg.inv(cv), d)
print(__doc__)
for chip in ('ibm_kingston', 'ibm_marrakesh', 'ibm_fez'):
    rec = json.load(open(os.path.join(CTRL, 'iqjob_%s.json' % chip))); z = np.load(os.path.join(CTRL, 'iqraw_%s.npz' % chip)); E = rec['E']
    i0 = [i for i, m in enumerate(rec['meta']) if m['arm'] == 'CAL' and m['state'] == '0']; i1 = [i for i, m in enumerate(rec['meta']) if m['arm'] == 'CAL' and m['state'] == '1']
    C0 = np.concatenate([z['c_%d' % i] for i in i0]); C1 = np.concatenate([z['c_%d' % i] for i in i1])
    mu0 = C0.mean(0); mu1 = C1.mean(0); u = mu1 - mu0; proj = lambda x, q: np.real((x - mu0[q]) * np.conj(u[q])) / np.abs(u[q]) ** 2
    rows = []
    for q in E:
        c0 = cov2(C0[:, q][proj(C0[:, q], q) < 0.5]); c1 = cov2(C1[:, q][proj(C1[:, q], q) > 0.5]); x = C1[:, q]
        in0 = m2(x, mu0[q], c0) < D2; in1 = m2(x, mu1[q], c1) < D2
        rows.append((q, np.mean(in0 | ~in1), np.mean(~in1), np.mean(in0), np.mean(in0 & in1), np.mean(proj(x, q) < 0.5)))
    r = np.array([x[1] for x in rows])
    print('=' * 92); print('%s: cal-0 pubs %s, cal-1 pubs %s, %d prepared-1 shots per qubit, %d measured qubits' % (chip, i0, i1, len(C1), len(E)))
    print('  PLAN rule over measured qubits: mean %.1f%%, median %.1f%%, max %.1f%%, > 20%%: %d' % (100 * r.mean(), 100 * np.median(r), 100 * r.max(), (r > 0.2).sum()))
    for lab, k in (('PLAN', 1), ('A not-in-1', 2), ('P reads-0', 5)):
        top = sorted(rows, key=lambda x: -x[k])[:6]
        print('  top 6 by %-11s: ' % lab + ', '.join('q%d %.1f%%' % (x[0], 100 * x[k]) for x in top))
    print('  %5s %7s %7s %7s %7s %7s   (top 6 by PLAN rule, all columns)' % ('qubit', 'PLAN', 'A', 'B', 'BOTH', 'P'))
    for x in sorted(rows, key=lambda x: -x[1])[:6]: print('  %5d %6.1f%% %6.1f%% %6.1f%% %6.1f%% %6.1f%%' % (x[0], *[100 * v for v in x[1:]]))
