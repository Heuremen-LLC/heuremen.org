# Addendum to PLAN_LENGTH.md -- written after the three jobs were submitted and BEFORE any result was read (22 Oct 2026)

1. The two instructions are not known to differ only in length. IBM's calibration record saved at firing (backend_properties_<chip>.json)
   gives measure and measure_2 different discriminator thresholds (kingston -3.79e6 vs -2.78e6; marrakesh -9.07e7 vs -2.04e7; fez
   -9.43e6 vs -3.78e6) and, on marrakesh and fez, different stated errors (0.0160 vs 0.0220; 0.0101 vs 0.0153); kingston states the same
   error for both (0.00928). So the integrated signal, and possibly the pulse amplitude, differ between the two. R2 and R3 therefore read
   how the kick depends on WHICH readout instruction is used, of which the stated length is one known difference, not the only one.
2. The plan's mechanism wording is narrowed. A ratio of 1 (PER-READOUT) means the kick does not depend on the instruction; it is what a
   phase-frame update applied once per readout would give, and also what any effect fixed per readout would give -- for example the
   resonator's ring-down after the pulse, whose duration is set by the resonator, not the pulse. A ratio at L2/L1 (PROPORTIONAL) is what a
   frequency shift during the pulse at equal amplitude would give, and could also come from an amplitude difference. NEITHER or CANNOT_TELL
   are read as they fall. No reading of this run alone names the mechanism; it narrows which ones remain.
3. The readings, code and outcomes are unchanged (analyze_len.py as hashed).
