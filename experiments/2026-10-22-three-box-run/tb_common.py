"""tb_common.py -- three-box weak-value run, 22 Oct 2026. Shared circuit builder and analysis (used by run, synth and analyze)."""
import math, numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import RYGate, StatePreparation
EPS = (0.10, 0.15, 0.20, 0.25, 0.30, 0.40)
BOXES = ('00', '01', '10', '11')                     # basis state of (s1 s0); every circuit probes one box
# variant -> (preselected support, postselection signs); weak value textbook w_b = sign_b for b in support, 0 otherwise
VARIANTS = {
    'V11': {'support': ('00', '01', '11'), 'sign': {'00': 1, '01': 1, '11': -1}},   # negative box at the corner (as 24 Sep / 11 Oct)
    'V10': {'support': ('00', '01', '10'), 'sign': {'00': 1, '01': 1, '10': -1}},   # standard Aharonov-Vaidman: negative box at |10>
    'VA':  {'support': ('00', '01', '11'), 'sign': {'00': -1, '01': 1, '11': 1}},   # same support as V11, negative box moved to |00>
}
def vec(support, sign=None):
    v = np.zeros(4, complex)
    for b in support: v[int(b, 2)] = (sign or {}).get(b, 1)
    return v / np.linalg.norm(v)
def textbook(var):
    V = VARIANTS[var]; psi = vec(V['support']); phi = vec(V['support'], V['sign']); den = np.vdot(phi, psi)
    return {b: float(np.real(np.conj(phi[int(b, 2)]) * psi[int(b, 2)] / den)) for b in BOXES}
def circuit(var, box, eps):
    """q0 = s0, q1 = s1, q2 = probe. Pre-select psi, rotate probe RY(2 eps) only if the system is in `box`, rotate phi -> |00>, read probe in X.
    Kept shots: s0 = s1 = 0. Weak value of the box ~ <X_probe>/(2 eps) on kept shots, in the weak limit."""
    V = VARIANTS[var]; qc = QuantumCircuit(3, 3, name=f'{var}_{box}_{eps:.2f}')
    qc.append(StatePreparation(vec(V['support'])), [0, 1])
    b1, b0 = int(box[0]), int(box[1])
    if b1 == 0: qc.x(1)
    if b0 == 0: qc.x(0)
    qc.append(RYGate(2 * eps).control(2), [1, 0, 2])
    if b1 == 0: qc.x(1)
    if b0 == 0: qc.x(0)
    qc.append(StatePreparation(vec(V['support'], V['sign'])).inverse(), [0, 1])
    qc.h(2); qc.measure([0, 1, 2], [0, 1, 2]); return qc
def keys():
    return [(v, b, e) for v in VARIANTS for e in EPS for b in BOXES]
def xp_from_bits(bits):
    """bits: (shots,) packed ints, bit0 = s0, bit1 = s1, bit2 = probe. Returns +-1 probe X values on kept shots."""
    k = ((bits & 1) == 0) & (((bits >> 1) & 1) == 0)
    return 1.0 - 2.0 * ((bits >> 2) & 1)[k]
def fit_SP(eps, y, s):
    """weighted LS: <Xp> = S + P sin(2 eps). Returns S, P, chi2/dof."""
    e = np.asarray(eps); X = np.column_stack([np.ones_like(e), np.sin(2 * e)]); W = 1 / np.asarray(s) ** 2
    A = X.T @ (X * W[:, None]); b = np.linalg.solve(A, X.T @ (W * np.asarray(y))); r = np.asarray(y) - X @ b
    return float(b[0]), float(b[1]), float((W * r ** 2).sum() / max(1, len(y) - 2))
