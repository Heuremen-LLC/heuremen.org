# Fresh-clock run: does the readout kick fade after many readouts? Read with clocks that are alive.
Written, SHA-256 hashed and pushed to GitHub BEFORE any job is submitted. Heuremen LLC (Wayfinder and Bones), 22 October 2026.
Why: the count run read the kick over readouts 27 to 81 from untouched qubits whose coherence had fallen to the shot-noise floor (67 to 97 of
117 below 5 x noise at 54 to 81 readouts); a pure-noise phase step has median 0, so that reading cannot tell a fade from dead clocks.
Design (freshclock.py; ibm_kingston, ibm_marrakesh, ibm_fez): control-run layout. Untouched qubits parked in 0, no gate, while the 39
measured qubits do a PRELOAD of 0, 27 or 72 readouts, or 72 equal idles. Then SX on the untouched qubits (a fresh clock) and a WINDOW of
9 readouts or 9 equal idles; untouched qubits read in X and Y. 2000 shots, every circuit twice per job. Dry run 24-28 s per chip.
Analysis (analyze_freshclock.py, analyze.py rulers unchanged): per untouched qubit, kick per readout = wrap(phase(window readout) -
phase(window idle)) / 9; size = kick x sign of its median after no preload (sign printed). The window's idle coherence is printed per arm.
RF.1 size after none, 27 readouts, 72 readouts, 72 idles. RF.2 size after 72 readouts minus after none: BELOW -> the kick fades;
CONSISTENT -> no fade resolved, interval is the bound; ABOVE -> grows. RF.3 the same after 27. RF.4 after 72 readouts minus after 72 idles
(the readouts, not the elapsed time). Four outcomes; copies by their gap; chips by value and width.
Synthetic test (synth_freshclock.py): no fade -> RF.2 +0.02738 [-0.03526, +0.09026]; planted fade of 1.2 -> RF.2 -1.23 BELOW;
no kick -> every reading CONSISTENT.
