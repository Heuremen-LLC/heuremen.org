# The sweep: every per-qubit and per-edge figure IBM publishes, measured the way IBM documents measuring it
Pre-registered. Written and SHA-256 hashed (PLAN_SWEEP.sha256, CODE_SWEEP.sha256) BEFORE any job is submitted. Heuremen LLC, Wayfinder
and Bones, 23 October 2026. Circuits, fire path and analysis: sweep.py (one file). Synthetic test: synth_sweep.text.
Chips: ibm_fez, ibm_kingston, ibm_marrakesh, ibm_boston, ibm_pittsburgh, ibm_phoenix (Darkwave Disco). One job each, IBM's record saved at
submission. T1 and T2 are measured in ../ladder_run and ../hahn_t2, zz in the direct-zz run; this job covers every other figure.

## What IBM states, and how this job reads it (IBM documentation rows from experiments/ibm-documentation/claims.md)
  init_error (row 25: |1> population after the default rep-delay when the prior shot prepared |1>): circuits I1 [measure a; X; measure b]
    and I0 [measure a; measure b]; 'a' of I1 reads each qubit after a shot that ended in |1>, 'a' of I0 after one that ended in |0>.
    Delivered init population is bounded [(a_I1 - a_I0)/(1 - a_I0 - p01), a_I1/(1 - p01)]; per qubit ABOVE / BELOW / CONSISTENT with
    1.96 shot widths; counts per chip.
  prob_meas1_prep0, prob_meas0_prep1, readout_error (rows 19, 20, 26, 27): IBM does not say what the shot before its readout shot ended in,
    so each delivered figure is a BRACKET over the two cases (a third circuit IQ [measure a; X; measure b; X] reads |1> prepared after a
    |0>-ending shot), widened by 1.96 shot widths. Reading: per qubit, how far the stated value sits outside the bracket (0 inside); median
    over qubits with a 4000-draw bootstrap; STATED BELOW DELIVERED / STATED ABOVE DELIVERED / CONSISTENT. 4000 shots, two copies.
  sx gate_error (rows 55, 60: RB measured simultaneously on all qubits; SX, ID, X errors assumed equal): simultaneous single-qubit Clifford
    RB on every qubit, independent random Cliffords, lengths 1, 150, 500, 1000, three sequences each, 300 shots; fit A p^m + 1/2; error
    per pulse = (1-p)/2 / mean pulses per Clifford (1.000 in this rz/sx synthesis). Reading: median delivered/stated over qubits, four outcomes.
  cz gate_error (row 59: RB in isolation, batches with a minimum separation between edges, alternating layers of single-qubit Cliffords
    and two-qubit gates, the value including the single-qubit layer): edges in batches whose nearest qubits are at least 3 links apart,
    layers [random Clifford on each qubit, CZ], lengths 1, 10, 30, 60, 100, three sequences, 300 shots, inverse two-qubit Clifford; fit
    A p^m + 1/4; error per layer 3/4 (1-p). Reading: median delivered/stated over edges, four outcomes; edges stated 0 and stated 1 listed
    with their delivered error.
  Not measured, and why: layer fidelity / EPLG (rows 62, 63: IBM reports the best of six pre-selected 100-qubit chains, which are not
  published, so no like-for-like chain exists); rzz error (row 57: a copy of the cz value in the record); gate and readout lengths
  (configuration, not performance).

KNOWN-GOOD: real-hardware readings whose answer is known, in each job: (1) single-qubit RB survival at m = 1 median >= 0.95 and the
two-qubit layer survival at m = 1 median >= 0.85 (state preparation, one Clifford and readout only); (2) median P(1) for a qubit read
after a |0>-ending shot <= 0.05. A chip failing either has its readings printed and not used.

## Bill rule (written now)
A figure goes on the bill only where IBM's stated value is BETTER than delivered (sx or cz delivered/stated ABOVE; readout STATED BELOW
DELIVERED; init ABOVE on a majority of qubits), on a chip with known-good met, AND a repeat against a fresh record reads the same.
An edge or qubit IBM states as exactly 0 error whose delivered error is resolved above 0 is a line on its own (a 0 is not a measurement).
Figures where IBM is right or conservative are reported on heuremen.org and do not go on the bill.

## Synthetic test (synth_sweep.text): Aer stabilizer simulation of a 23-qubit heavy-hex piece, every figure planted
match (all = stated) -> all CONSISTENT; worse (sx, cz 2.5x) -> sx 2.48, cz 2.56 ABOVE; better (0.4x) -> 0.40, 0.40 BELOW; inithigh (6x) ->
23 of 23 ABOVE; initlow (0.3x) -> 9 BELOW, 14 CONSISTENT, 0 ABOVE; rohigh (readout 3x) -> all three readout figures STATED BELOW DELIVERED.
An edge planted at true error 4e-3 and stated 0 reads 4.1e-3 (match world).
