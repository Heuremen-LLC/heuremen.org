# Dose repeat on ibm_fez: 10 of 39 read out (22 Oct 2026, night)
Written, hashed and pushed BEFORE the job. Why: in the 21 Oct control run, fez's kick with 10 of the 39 measured qubits read out was
CANNOT_TELL by the committed rule (copies disagree), pooled -0.133 [-0.226, -0.030] deg per readout. A repeat with the committed code decides it.
Design: the UNCHANGED control-run run.py (21 Oct, PREREG.sha256 unchanged) fired through fire_open.py with --chip ibm_fez --arms S,C on the
open instance of .ibm_quantum_key_3 (the 21 Oct job used Darkwave Disco; instance name is the only change). Arm S: 39 read, idle/readout at
N = 2, 4, 8; arm C: K = 1 and 10 of the 39 read (the same seeded K10 set), N = 8; settings XX, YY; every pub twice (copy 0, copy 1), shuffled.
44 pubs, 44000 shots, estimated 13 s. IBM record saved at firing.
Reading: dose_read.py, the committed R1/R2 code copied unchanged (only arms S and C read): R2 kick with K = 10 measured, copies 0/1/pooled,
CANNOT_TELL where copies disagree; plus the per-qubit copy gap.
Rule for the paper and bill: this run alone decides fez K = 10. If its committed outcome is BELOW (copies agree), fez joins kingston and
marrakesh as a resolved dose reading; if CANNOT_TELL again, it stays so and both runs are printed; if CONSISTENT, the bound is printed.
KNOWN-GOOD: the 39-read kick (R1, arm S) on fez must read BELOW zero, as in the 21 Oct control run (-1.62 [-2.08, -1.18] deg per readout);
if not, the run is not decided.
Synthetic test: synth_dose.text (planted none / proportional / global kicks read back by dose_read.py).
