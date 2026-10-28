# Pair run with raw readout: does the entanglement cost survive (1) with every measured neighbour verified in 0 and (2) against extra time?
Written, SHA-256 hashed and pushed to GitHub BEFORE any job is submitted. Heuremen LLC, 21 October 2026. Circuits: run_pairs_iq.py (run.py's
builder unchanged). Analysis: analyze_pairs_iq.py, tested on synthetic data and committed before submission.

## Why this run
The control run found Bell pairs among the untouched qubits keep less concurrence after mid-circuit readouts than after an equal idle. Two
ordinary causes are not excluded for the pairs: a measured neighbour leaked to a level that reads as 0, and extra elapsed time in the readout steps.

## Design (per processor: ibm_kingston, ibm_marrakesh, ibm_fez; two jobs, copy 0 and copy 1; a job is skipped if the remaining time is short)
- Roles, pairs (48 Bell pairs on coupled untouched qubits) and singles as in the control run; meas_type = kerneled (raw IQ for every readout).
- Calibration: every qubit in 0, and every qubit in 1, read once (2000 shots each).
- Pairs: idle N = 0; idle (D) and readout of all 39 (M, K = 39) at N = 2, 4; readout of ONE measured qubit (K = 1, the same qubit as the
  control run's dose arm; the other 38 idle exactly L) at N = 4. Settings ZZ, XX, YY, XY, YX. 4000 shots per circuit; order shuffled per job.

## Classification (identical to PREREG_IQ.md)
mu0, mu1 per qubit from calibration; projection t onto mu1 - mu0; final readouts bit = 1 if t > 0.5. A measured-qubit mid readout is
"verified 0" if its Mahalanobis distance to the 0 cluster (calibration-0 shots with t < 0.5) is below 3.5.

## Readings (four outcomes; 95 % bootstrap over pairs, 2000 resamples; copy 0, copy 1 and pooled; copies that disagree -> CANNOT_TELL)
Concurrence: X-state form from the five settings. Neighbours of a pair = measured qubits wired to either member (cz in IBM's calibration
saved at firing). Pairs whose idle concurrence is 0 at that N are excluded (count printed).
R1 per processor, N = 2, 4: median over pairs with at least one wired measured qubit of
   concurrence(idle, all shots) - concurrence(K = 39 readout, only shots where every neighbour was verified 0 on every readout).
   ABOVE 0 -> the entanglement cost survives with every neighbour verified in 0: a leaked neighbour is EXCLUDED as its cause.
R2 per processor, N = 4: median over pairs of loss(K = 39) - loss(K = 1), each loss = idle concurrence minus readout concurrence (all shots).
   Same schedule in both, so ABOVE 0 -> extra elapsed time per step (and anything tied to the schedule) is EXCLUDED as the cause.
R3 the same as R1 with all shots, for comparison. R4 fraction of shots kept by the R1 selection. R5 fraction of mid readouts outside both clusters.
