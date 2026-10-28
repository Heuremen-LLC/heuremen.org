# QIA Portfolio Run: does choosing qubits from IBM's published error figures cost a customer-style workload?
Pre-registered. Written and SHA-256 hashed (PLAN_DEMO.sha256, CODE_DEMO.sha256) BEFORE any job. Heuremen LLC, Wayfinder and Bones,
23 October 2026. Code: demo.py; planted-world test: synth_demo.py (synth_demo.text). Chips: ibm_kingston, ibm_fez (Darkwave Disco).

Workload: 10-qubit QAOA-shaped circuit (MaxCut on a line, two layers, fixed angles), the shape of the portfolio and routing pilots run on
IBM hardware, followed by its exact inverse behind a barrier: the correct output is all zeros, so success is known without a classical model.
Three placements of the same circuit, compiled with identical settings (optimization level 1, fixed layout), 4000 shots x 2 copies each:
  QISKIT  Qiskit's own choice (preset pass manager, optimization level 2, the default) from IBM's published figures
  IBM     the best 10-qubit line by IBM's published sx, cz and readout figures (record at firing), one fixed product formula
  OURS    the best 10-qubit line by our measured figures (sweep run 1 and repeat, about 1 h old at firing), the same formula
  KG_*    the same placements with measurement only.

KNOWN-GOOD: readout-only success at least 0.80 on all three placements; otherwise every reading is CANNOT_TELL.
Readings (success ratio, 4000-draw bootstrap over shots, four outcomes):
  R1 OURS / IBM  -- primary: same search, same formula, only the numbers differ.
  R2 OURS / QISKIT;  R3 IBM / QISKIT (Qiskit's search against the plain formula on the same IBM numbers).
  R4 per placement: delivered success against the success predicted from IBM's figures and from ours (descriptive).
Disclosed before the run: our cz figure is a layer error that already contains the single-qubit layer, so our predictions double-count
single-qubit error (they read low); both picks are the maximum over 1674 (kingston) candidate lines, so both predictions carry selection
optimism; a mirror circuit can cancel part of a coherent error. Our figures are about 1 h old at firing; IBM's are its record at firing.
Use: a demonstration for QIA (indahl.ai), reported with every outcome; not a bill line.
Synthetic (synth_demo.text): IBM-true world -> R1 CONSISTENT; ours-true world -> R1 ABOVE; equal world -> 5 ABOVE, 3 BELOW in 200 draws.
