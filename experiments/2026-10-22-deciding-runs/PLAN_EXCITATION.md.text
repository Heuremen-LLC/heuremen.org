# Claim 4 run: is the readout kick bigger when the measured qubits are in 1?
Written, SHA-256 hashed and pushed to GitHub BEFORE any job is submitted. Heuremen LLC (Wayfinder and Bones), 22 October 2026.
Why: claim 4 rested on two chips, one copy per job, an analysis written after the data, and a plan with no outside timestamp. On that data the
gap (check J, committed before computing) reads kingston +0.204 deg/step [-0.092, +0.320], marrakesh -0.673 [-0.982, -0.403].
Circuits: excitation.py (control-run layout and builder; runs.py fire path). Analysis: analyze_excitation.py (analyze.py's rulers unchanged),
tested on synthetic data (synth_excitation.py): planted gap -0.6 read back -0.5764; null: no false ABOVE/BELOW (one copy flag, CANNOT TELL).
Design: ibm_fez, ibm_kingston, ibm_marrakesh; 39 measured, 117 untouched spectators on the equator; idle D and readout M with the measured
qubits in 0 or in 1; N = 0, 4, 8; X and Y; 2000 shots; every circuit twice per job, shuffled. Dry run: 24 s per chip.
Readings, per chip then combined by value and width (analyze.py summary): phase slope per step over N = 0, 4, 8, per spectator.
RJ.1 kick in 0 (M0 - D0); RJ.2 kick in 1 (M1 - D1); RJ.3 G = RJ.2 - RJ.1: ABOVE -> bigger in 1, BELOW -> smaller in 1, CONSISTENT -> no
difference resolved, the interval is the bound; RJ.4 measured qubits in 1 without readout (D1 - D0). Four outcomes; copies by their gap.
