# Readout-instruction repeat on ibm_kingston (22 Oct 2026, night), pushed BEFORE the job
Why: the 22 Oct length run on kingston read the measure_2 / measure kick ratio at 4 readouts as PER-READOUT (0.920 [0.774, 1.169]) but at
8 readouts as CANNOT_TELL because its two copies disagreed (copy gap [+0.003, +0.091]). The bill's phase-kick line says "kingston not decided".
Design: the UNCHANGED run_len.py / fetch_len.py / analyze_len.py of the 22 Oct length run (CODE hashes below match the published ones),
fired again on ibm_kingston with the default 1000 shots on the open instance of .ibm_quantum_key, from this folder (the published 22 Oct
files are not overwritten). Estimated 15 s.
Reading: analyze_len.py, its committed outcomes (R1 kick per length; R2 ratio k(L2)/k(L1) at N = 4 and 8 against 1 and against L2/L1 =
0.633, copies compared). This run alone decides kingston; both runs are printed side by side.
fez is not repeated tonight: its proportional ratio (0.807) is too close to 1 for 1000 shots to separate; that needs ~4000 shots (58 s).
KNOWN-GOOD: as in PLAN_LENGTH.md, the measure (L1) kick at 8 readouts must read below zero on kingston (22 Oct: -2.755 [-3.623, -2.067]);
if not, the run is not decided.
Synthetic test: synth_len_results.text (the 22 Oct synthetic test of the same code).
