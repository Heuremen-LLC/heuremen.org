"""How much coherence the count run's phases were read from (after the data; reported with the cannot-tell readings).
Per chip and N: median coherence (idle, readout; both copies, 2000 shots), spectators whose readout coherence is below 5 x the shot noise
(1 / sqrt(2000)), and the phase noise at the median coherence in degrees."""
import numpy as np, analyze as A
for chip in ('ibm_fez', 'ibm_kingston', 'ibm_marrakesh'):
    J = A.Job('count', chip); J.props_n = 156; nS = len(J.info['S'])
    co = lambda mo, N: np.hypot(J.mean_pm((0, 1), 'fin', nS, mode=mo, N=N, setting='XX'), J.mean_pm((0, 1), 'fin', nS, mode=mo, N=N, setting='YY'))
    for N in (1, 9, 27, 54, 81):
        d, m = co('D', N), co('M', N); noise = 1 / np.sqrt(2000)
        print('%s N=%2d median coherence idle %.3f readout %.3f | below 5x shot noise: %d of %d | phase noise at median %.1f deg'
              % (chip, N, np.median(d), np.median(m), np.sum(m < 5 * noise), nS, np.degrees(noise / np.median(m))))
