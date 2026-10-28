"""Sign convention check. (1) exact gate algebra of the 7 Oct hand circuits vs the paper's X / Y settings, Qiskit gate definitions, numpy only.
(2) the 21 Oct control-run kingston kick re-computed from its raw bits through the SAME estimator functions as the 7 Oct data (est.py), compared with the
committed results_prereg.txt line. Prints; imported by analyze_oct07.py."""
import os, sys, json, re
import numpy as np
from est import *
def run():
    I2 = np.eye(2, dtype=complex); Z = np.diag([1, -1]).astype(complex)
    H = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2); SDG = np.diag([1, -1j]); SX = 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]])
    RZ = lambda a: np.diag([np.exp(-1j * a / 2), np.exp(1j * a / 2)])
    zero = np.array([1, 0], complex)
    ez = lambda U, psi: float(np.real(np.vdot(U @ psi, Z @ (U @ psi))))      # mean(1 - 2 bit)
    print('SIGN CHECK 1: exact gate algebra, state = SX|0> precessed by alpha about z (7 Oct clock: sx, wait, hand), Qiskit gate definitions')
    print('  Paper settings (control_run/run.py lines 23-24): XX = h then measure; YY = sdg, h then measure; phase = atan2(<Y>, <X>) = atan2(e_YY, e_XX).')
    print('  7 Oct hands (look_floor_*.py lines 130-131): hand X = sx then measure; hand Y = rz(pi/2), sx then measure.')
    print('  %8s | %12s %12s | %12s %12s | %14s %14s' % ('alpha', 'e_hand_X', 'e_hand_Y', 'e_XX(paper)', 'e_YY(paper)', 'paper phase', 'hands->phase'))
    worst = 0
    for al in (0, 20, 45, 90, 135, -30, -100, 170):
        a = np.radians(al); psi = RZ(a) @ SX @ zero
        ehX, ehY = ez(SX, psi), ez(SX @ RZ(np.pi / 2), psi)           # circuit order: rz first then sx => matrix SX @ RZ
        eXX, eYY = ez(H, psi), ez(H @ SDG, psi)                     # sdg first then h => matrix H @ SDG
        pp = np.degrees(np.arctan2(eYY, eXX)); xs, ys = paper_xy_from_hands(ehX, ehY); ph = np.degrees(np.arctan2(ys, xs))
        worst = max(worst, abs(((pp - ph + 180) % 360) - 180))
        print('  %8.1f | %+12.5f %+12.5f | %+12.5f %+12.5f | %+14.4f %+14.4f' % (al, ehX, ehY, eXX, eYY, pp, ph))
    print('  max |paper phase - hands->phase| = %.2e deg.  Precession by +alpha moves the paper phase by +alpha (=-90+alpha here), and the hand mapping reproduces it exactly: SAME sign convention.' % worst)
    print('  Consequence: the history_inventory reader used H = (2 lvX - 1) - i (2 lvY - 1) with lv ~ P(1): H = -<Y> + i<X> = i exp(i phi), so its "kick" = -arg(H_look/H_dark) = -(phi_M - phi_D) = MINUS the paper-convention gap.')
    print('  A positive inventory kick is therefore a NEGATIVE kick in the paper convention, i.e. the same direction as the 21-22 Oct kick (all negative there).')
    # ---- check 2
    print('SIGN CHECK 2: 21 Oct control run kingston, N = 8, arm S, pooled copies, through est.py (hands mapping applied to the 21 Oct <X>,<Y> and back)')
    CR = os.path.join(ROOT, 'control_run'); _a = sys.argv; sys.argv = sys.argv[:1] + [CR]
    import importlib.util
    sp = importlib.util.spec_from_file_location('ctl_analyze', os.path.join(CR, 'analyze.py')); C = importlib.util.module_from_spec(sp); sp.loader.exec_module(C); sys.argv = _a
    rec = json.load(open(os.path.join(CR, 'job_ibm_kingston.json'))); z = np.load(os.path.join(CR, 'raw_ibm_kingston.npz')); meta = rec['meta']; nS = len(rec['S'])
    def expv(mode, N, st):
        ii = [i for i, m in enumerate(meta) if m.get('arm') == 'S' and m.get('mode') == mode and m.get('N') == N and m.get('setting') == st]
        a = np.concatenate([z['fin_%d' % i] for i in ii]); return np.array([(1 - 2 * C.bits(a, k)).mean() for k in range(nS)])
    out = {}
    for mode in ('M', 'D'):
        x, y = expv(mode, 8, 'XX'), expv(mode, 8, 'YY')
        eXh, eYh = y, x                                        # what the 7 Oct hands would have read (inverse of the mapping)
        xs, ys = paper_xy_from_hands(eXh, eYh); out[mode] = phase(xs, ys)
    tot = total_gap({8: out['M']}, {8: out['D']}, [8]); k = slope0(tot, [8])
    print('  21 Oct kingston, recomputed through est.py: median kick = %+.4f deg per readout, n=%d' % (np.median(k), len(k)))
    txt = open(os.path.join(CR, 'results_prereg.txt')).read().split('=' * 100)
    blk = [b for b in txt if b.strip().startswith('ibm_kingston')][0]
    line = [l for l in blk.split('\n') if 'R1 kick' in l and 'pooled' in l][0]; print('  committed results_prereg.txt:  ' + line.strip())
    cm = float(re.search(r'pooled\s+([+-][0-9.]+)', line).group(1)); print('  difference %.2e deg -> %s' % (abs(np.median(k) - cm), 'IDENTICAL (code path and sign agree)' if abs(np.median(k) - cm) < 5e-4 else 'MISMATCH'))
    return abs(np.median(k) - cm) < 5e-4
if __name__ == '__main__': run()
