# Two checks on the control run, defined AFTER its results were seen and BEFORE these two were computed
Heuremen LLC, 21 October 2026. Rules fixed here; script checks_after_data.py; both committed before it is run.
Wiring: the processor's two-qubit gate list (cz) in backend_properties_<chip>.json, saved at firing.
Four outcomes as in PREREG.md; 95 % bootstrap over qubits (or pairs), 2000 resamples; arms S (coherence) and P (concurrence).

## Check A: can measured qubits that are stuck (excited or leaked) account for the superposition and entanglement costs?
- For every measured qubit e and every readout circuit with N >= 2, from its own N back-to-back mid-circuit readouts:
  f_e = fraction of shots in which e reads 1 on two consecutive readouts (conservative: counts misreads that repeat by chance).
  q_e = mean over consecutive pairs of [P(1 then 1) - P(1)^2] (the persistent part; misreads are independent between readouts).
- Bound: a neighbour that is stuck in a fraction f of shots and acts on the untouched qubit through a phase (a diagonal coupling)
  can lower its coherence, or a pair's concurrence, by at most 2f, whatever the size of that phase. For an untouched qubit (or a pair),
  B = 2 * sum of f_e over the measured qubits wired directly to it (to either member of a pair). B_est is the same with q_e.
- Reading, per processor, per N, for coherence (arm S, N = 2, 4, 8) and concurrence (arm P, N = 2, 4):
  median over qubits (pairs) of (loss - B), where loss = idle value minus readout value (positive = lost more during readouts).
  ABOVE 0 -> stuck neighbours cannot account for the loss even on the conservative bound (EXCLUDED).
  Otherwise, the same with B_est: ABOVE 0 -> excluded on the persistent estimate only (CANNOT TELL on the conservative bound).
  Otherwise -> NOT EXCLUDED.
- Untouched qubits with no measured qubit wired to them (bound 0 under the diagonal-coupling assumption): their median loss vs 0, reported apart.
- Assumption stated: a stuck neighbour acts through a phase; energy exchange with it is not bounded here.

## Check B: where does the kick land?
- Split the 117 untouched qubits into NEIGHBOURS (wired directly to at least one qubit that was read out in that arm) and the rest.
- For K = 1, 10 (arm C) and K = 39 (arm S), N = 8, pooled copies: median kick per readout in each group, vs 0.
- Outcome: NEIGHBOUR-LOCAL if the non-neighbour median is CONSISTENT with 0 while the neighbour median is clear of 0;
  CHIP-WIDE if the non-neighbour median is clear of 0; CANNOT TELL otherwise. Reported per K and processor.
- Also printed: the neighbour median for K = 10 and K = 39 side by side (a local effect predicts similar values per neighbour).
