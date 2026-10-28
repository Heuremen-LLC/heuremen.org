# Addendum: second readout-instruction repeat (22 Oct 2026, ~22:20 EDT), pushed BEFORE these jobs
Why: after the first repeat, kingston is decided at 8 readouts (ratio 1.01 [0.80, 1.37], length ratio 0.63 excluded) but CANNOT_TELL at 4
(copies disagree); fez is CANNOT_TELL at 4 and 8 in the 22 Oct run because its length ratio (0.807) is close to 1 at 1000 shots.
Design: the UNCHANGED run_len.py / fetch_len.py / analyze_len.py (hashes as CODE_LENGTH_REPEAT.sha256), only the shot count changed:
  kingston: 2000 shots, open instance of .ibm_quantum_key (as the 22 Oct run), folder length_repeat/r2 ... (estimate ~29 s)
  fez: 4000 shots, open instance of .ibm_quantum_key_3 through fire_key3.py (KEYFILE swapped only) (estimate ~58 s)
Reading: analyze_len.py's committed outcomes. Each of these runs alone decides its processor; all runs are printed side by side; no pooling.
KNOWN-GOOD: as PLAN_LENGTH.md: the L1 (measure) kick at N = 8 must read BELOW 0 with an interval overlapping the control run's committed
kick (fez -1.62 [-2.08, -1.18], kingston -2.86 [-3.55, -2.09]); if not, that processor is not decided.
Synthetic test: synth_len_results.text (unchanged code).
