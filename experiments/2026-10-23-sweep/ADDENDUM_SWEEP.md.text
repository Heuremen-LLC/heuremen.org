# Addendum to PLAN_SWEEP.md (23 Oct 2026, after fez/kingston/marrakesh run 1 was read, before any repeat)
1. BILL SCOPE CORRECTED. The plan's line "Figures where IBM is right or conservative ... do not go on the bill" is withdrawn. A published
   figure that does not match the hardware is a bill line in either direction, with its direction stated. In scope from run 1:
   init_error stated ABOVE delivered (fez 12, kingston 7, marrakesh 30 qubits); init_error ABSENT from IBM's record (fez 40, kingston 31,
   marrakesh 9 qubits). Readout: every stated figure inside the delivered bracket on all three chips -- no line.
2. REPEATS. R2 = sweep.py unchanged, fired now through fire_sweep_rep.py (file names suffixed _r2) on fez, kingston, marrakesh: a
   same-record repeat (IBM re-dated 0-8 of 352 cz and 0-4 of 156 sx values since run 1). It tests our own reproducibility: a line needs
   run 1 and R2 to read the same outcome. R3 = the fresh-record repeat, fired when at least half of a chip's sx and cz gate_error dates
   (and, for init lines, init_error dates) are newer than run 1's record; the bill rule's "fresh record" is R3.
3. Per-qubit/per-edge lines (init above delivered, absent values) need the same unit to read the same in run 1 and R2; counts printed.
4. Staleness test run after run 1 (zero QPU, defined after the data): Spearman(age of IBM's value, delivered/stated) includes 0 on every
   chip and gate; the younger half of each chip already reads 1.22 to 1.71.
