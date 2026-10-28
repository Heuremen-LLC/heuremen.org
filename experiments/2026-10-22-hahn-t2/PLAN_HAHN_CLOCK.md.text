# Hahn T2 decay clock: delivered T2 on ibm_fez against hours since IBM's T2 calibration, one unchanged record (22 Oct 2026 night)
Written and pushed BEFORE the first clock point. Question (Wayfinder): does delivered T2 fall after calibration while the stated value stays?
Tonight's chips cannot say: chip and stated-T2 age are confounded (fez ~1 h ~1.0, marrakesh ~3 h, kingston ~14 h ~0.91-0.95), and within a chip
the ages barely spread (age_test.text, defined after the data). The clock holds the chip and the stated values fixed and moves only time.
Design: hahn.py ibm_fez --fire, unchanged, open instance of .ibm_quantum_key_3, at roughly 2, 4, 8, 16 h after fez's T2 calibration
(record 2026-10-22 18:17:23 EDT, T2 dates ~18:16), as seconds and presence allow. Each point keeps its own job/props/raw files (_clockN).
Runs 1 and 2 (0.7 h and 0.9 h) are points 0a and 0b. If IBM's record changes, the clock stops and that point is reported separately.
Reading: per point, hahn_pool.py (FIXED and FREE fits), median delivered/stated; then the slope of log(median ratio) against hours since
calibration, by least squares over points, with a bootstrap over qubits within each point (2000). Outcome on the slope: BELOW 0 (decays),
ABOVE 0, CONSISTENT, or CANNOT_TELL if the two fits disagree in sign with intervals excluding 0.
KNOWN-GOOD: median P(1) at tau = 0 at least 0.90 per point. Synthetic test: synth_hahn.text (unchanged code).
