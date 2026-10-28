# Hahn T2 repeat, kingston and marrakesh (22 Oct 2026, night)
Written and pushed BEFORE the repeat jobs. Code unchanged (CODE_HAHN.sha256). Run 1 files kept as *_run1.*.
Why: run 1 read kingston delivered/stated T2 0.910 [0.871, 0.962] BELOW, marrakesh 0.907 [0.873, 1.012] CONSISTENT. The vendor's
answers to a single job are "one job" and "T2 drifts"; a second job on the same chips with the same code answers the first and,
if IBM's record has updated between jobs, part of the second.
Design: hahn.py CHIP --fire, unchanged (open instance of .ibm_quantum_key_3), IBM record saved at firing.
Reading: hahn.py --analyze CHIP on run 2, outcome by the same rule. Then run 1 and run 2 combined per chip by value and width
(sigma = 95% width / 3.92), quoted only if chi2/dof < 2; otherwise both runs printed side by side and the chip reads CANNOT_TELL.
The kingston bill line stands as a deviation only if run 2 alone also reads BELOW.
KNOWN-GOOD: the median P(1) at tau = 0 must be at least 0.90, per chip, as in run 1.
Synthetic test: synth_hahn.text (unchanged code, unchanged worlds).
