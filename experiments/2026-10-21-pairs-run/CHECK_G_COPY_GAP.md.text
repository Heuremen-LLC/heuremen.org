# Check G: where the committed analyses printed "copies disagree", is there a real shift between the copies?
Defined after the raw-readout and pair results were seen, before this was computed. Heuremen LLC, 22 Oct 2026. Script check_g_copy_gap.py.
The committed rule compared the two copies' outcome LABELS, so a copy that is merely less precise counted as disagreeing. The copy gap itself is
the measurement: per qubit (or pair), value(copy 0) - value(copy 1), median over qubits (pairs), 95 % bootstrap.
- Interval contains 0 -> no shift between the copies is resolved; the pooled reading (as printed by the committed analysis) is the measurement.
- Interval excludes 0 -> a real copy-to-copy shift; the reading stays CANNOT TELL.
Applied to: raw-readout run R1 (fez N=2, 4; kingston N=2; marrakesh N=2); pair run R1 (fez N=2; kingston N=2, 4; marrakesh N=2) and R2 (fez).
The committed outcomes are not changed; this is reported beside them.
