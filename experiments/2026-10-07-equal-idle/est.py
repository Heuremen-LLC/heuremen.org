"""Shared estimator pieces for the 7 Oct re-read. The statistical functions (wrap, boot, boot2, outcome, reading, bfs) are IMPORTED UNCHANGED from
deciding_runs/analyze.py. total_gap() is copied verbatim from the nested total() in deciding_runs/analyze.py count() (phase gap M-D unwrapped along N by
adding wrapped increments). The only new code is data loading and the hand -> (<X>,<Y>) mapping."""
import sys, os
import numpy as np
ROOT = '/Users/heuremen/bones-forge/arxiv_one_readout'
sys.path.insert(0, os.path.join(ROOT, 'deciding_runs'))
_argv = sys.argv; sys.argv = sys.argv[:1]
import analyze as A          # deciding_runs/analyze.py, unchanged
sys.argv = _argv

def paper_xy_from_hands(eXh, eYh):
    """7 Oct hands -> the paper's <X>, <Y>. Derived in sign_check.py (exact gate algebra): hand 'X' (final sx) reads Bloch <Y>; hand 'Y' (rz(pi/2), sx) reads Bloch <X>."""
    return eYh, eXh          # (<X>, <Y>)

def phase(x, y): return np.arctan2(y, x)       # paper: phi = atan2(<Y>, <X>)

def total_gap(ph_M, ph_D, Ns):
    """verbatim logic of deciding_runs/analyze.py count().total: ph_M, ph_D dicts N -> phase array. Returns dict N -> unwrapped gap in degrees."""
    dphi = {N: A.wrap(ph_M[N] - ph_D[N]) for N in Ns}; tot = {Ns[0]: dphi[Ns[0]]}
    for a, b_ in zip(Ns[:-1], Ns[1:]): tot[b_] = tot[a] + A.wrap(dphi[b_] - dphi[a])
    return {N: np.degrees(v) for N, v in tot.items()}

def slope0(tot, Ns):
    """per-qubit slope of the unwrapped gap against N, line through the origin (the gap is 0 at N = 0 by construction); with one N it is gap/N, i.e. R1 / R2.1."""
    Nn = np.array(Ns, float); return sum(N * tot[N] for N in Ns) / np.sum(Nn ** 2)

def slope_icpt(tot, Ns):
    """per-qubit least-squares slope WITH intercept, as deciding_runs R2.3 fits (lstsq on [N, 1])."""
    A_ = np.stack([np.array(Ns, float), np.ones(len(Ns))], 1); return np.linalg.lstsq(A_, np.stack([tot[N] for N in Ns], 1).T, rcond=None)[0][0]
