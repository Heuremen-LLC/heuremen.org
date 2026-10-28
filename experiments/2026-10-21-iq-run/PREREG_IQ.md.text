# Raw-readout run: does the superposition (T2) cost of untouched qubits survive when every measured neighbour is verified in 0?
Written, SHA-256 hashed and pushed to GitHub BEFORE any job is submitted. Heuremen LLC, 21 October 2026. Circuits: run_iq.py (uses run.py's
circuit builder unchanged). Analysis: analyze_iq.py, tested on synthetic data and committed before submission.

## Why this run
The control run (same day) found that untouched qubits keep less coherence after mid-circuit readouts than after an equal idle. After-data
checks excluded, as its explanation: measured neighbours stuck in 1 (the loss survives in shots where every neighbour read 0), extra elapsed
time per readout step and any chip-wide cause (the loss is absent with 1 qubit read out and present with 39, on the same schedule; a direct
bound on kingston). The one ordinary cause left: a measured neighbour leaked to a higher level that the standard readout reports as 0.
This run returns every readout as a raw IQ point, so a neighbour can be required to sit inside its own calibrated 0 cluster.

## Design (one job per processor, order ibm_marrakesh, ibm_kingston, ibm_fez; a processor is skipped if the remaining time is short)
- Roles as in the control run: 39 measured qubits E, 117 untouched spectators S, all spectators unpaired on the equator.
- Calibration: every qubit prepared in 0 and read once; every qubit prepared in 1 (x) and read once.
- Arm S: idle N = 0; idle (D) and readout (M) at N = 2, 4, 8; spectators read in X and Y.
- Every circuit twice (copy 0, copy 1), all shuffled together; 2000 shots per circuit; meas_type = kerneled (raw IQ for every readout).

## Classification (fixed now)
- Per qubit q: mu0, mu1 = mean IQ point of calibration-0 and calibration-1 shots. The axis u = mu1 - mu0; a point's projection
  t = Re((x - mu0) * conj(u)) / |u|^2 (0 at mu0, 1 at mu1).
- Spectator final readouts: bit = 1 if t > 0.5, else 0.
- Measured-qubit mid readouts: "verified 0" if the point lies inside the 0 cluster: Mahalanobis distance to mu0 < 3.5 using the covariance of
  the calibration-0 shots with t < 0.5 (the 0 side). "Outside both" if the distance to the 0 cluster and to the 1 cluster (calibration-1 shots
  with t > 0.5) both exceed 3.5.

## Readings (four outcomes; 95 % bootstrap over spectators, 2000 resamples; copy 0, copy 1 and pooled; copies that disagree -> CANNOT_TELL)
Coherence = |<X> + i<Y>| from the classified spectator bits. Neighbours of a spectator = measured qubits wired to it (cz gates in IBM's
calibration saved at firing). Only spectators with at least one wired measured qubit enter R1 and R2.
R1 (the reading this run exists for): per processor and N = 2, 4, 8, median over spectators of
   coherence(idle, all shots) - coherence(readout, only shots in which every neighbour was verified 0 on every readout).
   ABOVE 0 -> the superposition cost survives with every measured neighbour verified in its 0 cluster: a leaked neighbour is EXCLUDED as its
   cause on that processor. CONSISTENT -> cannot tell (the loss may be leakage). BELOW -> reported as such.
R2: the same with all shots (no selection), for comparison.
R3: fraction of measured-qubit mid readouts "outside both" clusters, per processor and N (a direct count of leak candidates).
R4: fraction of shots kept by the R1 selection, per processor and N.
