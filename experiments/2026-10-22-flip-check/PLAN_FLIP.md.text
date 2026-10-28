# The X-flip check: is the coherence loss caused by a measured neighbour stuck in a level that reads as 0?
Heuremen LLC, Wayfinder and Bones, 22 October 2026. Written and SHA-256 hashed (PLAN_FLIP.sha256; code in CODE_FLIP.sha256) BEFORE any job
is submitted. Circuits: run_flip.py (control_run/run.py's `spectator_circuit`, unchanged, plus a thin wrapper). Analysis: analyze_flip.py,
tested on synthetic data (synth_flip.py) before submission. Instance: Darkwave Disco.

## Why this run
Claim 6: after mid-circuit readouts, untouched qubits keep less coherence than after an equal idle (0.004-0.021 less). The raw-readout run
(control_run/run_iq.py, PREREG_IQ.md, analyze_iq.py; site page experiments/2026-10-21-iq-run) kept only the shots in which every measured
neighbour's raw IQ point sat inside its own calibrated 0 cluster (Mahalanobis < 3.5) on every mid readout, and the loss survived.
One ordinary cause is still open. A measured neighbour may have leaked out of the 0-1 subspace into a level whose IQ point happens to fall
inside the 0 cluster. The selection cannot see that.

The X-flip check: after the last step, apply X to every measured qubit and read it once more. A qubit truly in 0 moves to the 1 cluster. A qubit
stuck in a level outside 0-1 does not move, because X is a resonant 0-1 pulse and is detuned from 1-2 by the anharmonicity (about 300 MHz).
That qubit reads where it was, inside the 0 cluster, or outside both clusters.

## Design (one job per processor; order ibm_kingston, ibm_marrakesh, then ibm_fez if the remaining time allows)
- Roles as in run_iq.py: 39 measured qubits E, 117 untouched spectators S, all spectators unpaired on the equator, read in X and Y.
- Spectator arm: idle N = 0; idle (D) and readout (M) at N = 4 and 8. Each step either idles E for one readout duration L or reads E.
- NEW, the flip: right after the barrier that closes the last step (or closes preparation when N = 0), X goes on every measured qubit. A barrier
  on all qubits fences it. Then the spectators' basis rotation, then one final readout layer: the spectators (register `fin`) and the measured
  qubits (register `efin`, NEW). This happens in BOTH arms, so the idle arm gives the baseline flip failure.
- The flip is verified in the transpiled circuit (`verify()` in run_flip.py; firing stops if it fails). In every spectator circuit, each measured
  qubit carries exactly one `x`, after its last mid readout and directly before its final measure. The negative test strips the x after
  transpiling: 22 of 22 circuits flagged. Transpiled sequence of a measured qubit at N = 4: `measure x4, x, delay, measure`. The delay is the
  spectator rotation, one single-qubit gate.
- Calibration: every qubit in 0, read once; every qubit in 1, read once (both as run_iq.py). NEW `0X`: every qubit in 0, read (register `c0`),
  X fenced by barriers, read again (`c1`). This places the post-flip cluster after one readout, as in the readout arm.
- Every circuit twice (copy 0, copy 1), all 26 shuffled together (seed 20261023). 3000 shots per circuit. meas_type = kerneled, so every
  readout comes back as a raw IQ point. On 21 Oct this mode already ran on all three processors with this circuit family (control_run/fire_iq.log).

### What changed vs run_iq.py
N = 0, 4, 8 (N = 2 dropped to buy shots); the X flip and the `efin` final readout of the measured qubits; the `0X` calibration; 3000 shots
(was 2000); a transpile-time verifier for the flip; output files flipjob_/flipraw_/flip_backend_properties_<chip>; the plan hash stored as
`plan_sha256`. The final readout layer now also reads the 39 measured qubits at the same time as the spectators. This is the same in both arms.

## Classification (fixed now; as analyze_iq.py where it overlaps)
- Per qubit q: mu0 and mu1 are the mean IQ points of calibration-0 and calibration-1. The axis is u = mu1 - mu0, and a point's projection is
  t = Re((x - mu0) conj(u)) / |u|^2. Spectator bit = 1 if t > 0.5.
- 0 cluster: mu0 with the covariance of calibration-0 shots that have t < 0.5. Verified 0 means Mahalanobis distance < 3.5.
- Post-flip cluster (primary): from `0X`, take the second-read points of shots whose first read was verified 0. Its mean is mu1x, and its
  covariance comes from those points with t > 0.5. A sensitivity line uses the calibration-1 cluster instead.
- Flipped = inside the post-flip cluster AND NOT inside the 0 cluster. Failed = inside the 0 cluster (even if also inside the 1 cluster), or
  outside both.
- Flip-readable guard: a measured qubit whose `0X` second read fails the flip rule in more than 20 % of verified-0 shots cannot report a flip.
  It is left out of F1 and F2, and every spectator wired to it is left out of L1, L2 and L3. On 21 Oct, calibration-1 points failed the rule
  as follows: kingston median 2.6 % (max 12 %); fez median 3.4 % (one qubit 24 %); marrakesh median 3.1 % (three qubits > 20 %, one 99 %).
- Neighbours of a spectator = measured qubits wired to it by a cz gate in IBM's calibration saved at firing.

## Readings (ruler as analyze_iq.py: four outcomes ABOVE / CONSISTENT / BELOW / CANNOT_TELL; 95 % bootstrap, 2000 resamples; copy 0,
## copy 1 and pooled; copies that disagree -> CANNOT_TELL; per processor and per N)
Coherence = |<X> + i<Y>| from the classified spectator bits. Per-spectator loss = coherence(idle) - coherence(readout).
- F1: fraction of (shot, measured qubit) readings that failed to flip, as the mean over flip-readable measured qubits. It is reported per
  processor, per N, for the idle arm and the readout arm, split into inside 0 and outside both, and also against the CAL1 cluster. For the
  readout arm it is also reported among shots in which that qubit was verified 0 on every mid readout.
- F2, the readout-induced hidden-leak rate: per measured qubit, F1(readout, shots verified 0 on every mid read) - F1(idle), same N; mean over
  qubits, bootstrap over qubits. F2u is the same difference over all shots. It also counts the visible excitations the mid reads already see.
- L1: loss with the previous selection. Idle arm = all shots; readout arm = shots in which every wired neighbour was verified 0 on every mid
  readout. Median over spectators (wired, all neighbours flip-readable). L1w is the same over all wired spectators, repeating run_iq's R1.
- L2: loss with BOTH that selection AND every wired measured neighbour flipping correctly at the end. The flip selection is applied to the idle
  arm too (symmetric), so the baseline flip failure, about 4 % per neighbour, cancels.
- L3 = L2 - L1, per spectator, then the median. This is a paired difference on the same spectators and the same data, so it is the sharpest
  reading here (null half-width about 0.0005 against about 0.003 for L1 and L2).
- Shots kept by the L1 and L2 selections (median over the spectators used).

## Decision rule (stated before data; per processor and N, on the FINAL state after the copy check)
- L2 ABOVE 0 -> a neighbour stuck in a level that reads as 0 is EXCLUDED as the cause of the loss, to the extent the flip test catches it.
- L2 CONSISTENT with 0 while L1 ABOVE -> the loss MAY come from such neighbours.
- Anything else is reported as is.
- Companion readings, stated now: L3 BELOW means the shots removed by the flip selection carried loss, and -L3 is the size of that part.
  L3 CONSISTENT means the flip selection removed no loss down to its CI. So if L2 reads CONSISTENT while L3 reads CONSISTENT, the L2 reading
  is a shortfall of power (fewer shots), not leak removal, and is reported that way. F2 ABOVE is direct evidence of readout-induced hidden leakage.

## Known-good checks (real hardware, the answer already known; analyze_flip.py prints each with PASS/FAIL)
KNOWN-GOOD: (1) the idle arm, which has no readouts, so no readout-induced leak is possible: at every N (0, 4, 8), each flip-readable measured
qubit's flip-failure rate minus its own CAL 0X failure rate, averaged over qubits with a 2000-resample bootstrap, must NOT read ABOVE in the
pooled FINAL. The idle arm must sit at the calibration's own floor. (2) The 'in 0 then X' calibration (CAL 0X): verified-0 shots must land in the
post-flip (1) cluster in nearly all shots, median over flip-readable measured qubits >= 95 %. On 21 Oct, calibration-1 points failed the rule at
a median of 2.6-3.4 %. (3) L1 must reproduce the raw-readout run's resolved clean-shot loss within its width. Wherever that run's R1 FINAL read
ABOVE (control_run/results_iq.txt, pooled), this run's L1w (the same spectator set and selection as R1) must fall inside 21 Oct's 95 % interval:
- ibm_kingston N=4: +0.0056 [+0.0013, +0.0105]
- ibm_kingston N=8: +0.0062 [+0.0013, +0.0155]
- ibm_marrakesh N=4: +0.0087 [+0.0032, +0.0160]
- ibm_marrakesh N=8: +0.0205 [+0.0112, +0.0284]
- ibm_fez is not resolved on 21 Oct (N=4 copies disagree, N=8 CONSISTENT), so it has no L1 known-good.
A FAIL on (1) or (2) means the flip instrument is not working on that processor, and its F2 and L2 are reported as CANNOT_TELL. A FAIL on (3)
is reported as a failure to reproduce, next to the readings.

## What the flip test cannot catch
1. A leaked state that X does move: anything inside the 0-1 subspace (an excited or rotated neighbour; the mid selection sees most of these),
   or a level the 0-1 pulse happens to drive.
2. Leakage that relaxes before the X. A neighbour can leak at mid readout k, cost coherence, and decay 2 -> 1 -> 0 before the flip. If it
   passes through 1 at a later mid read, L1 already removes the shot. If it decays within one step, nothing sees it. A leak at the LAST readout
   is caught, because the X follows within one gate.
3. A leaked level whose IQ point lands inside the post-flip (1) cluster and outside the 0 cluster. Then it reads as flipped. But it would also
   have read in the 1 cluster at the mid readouts, which L1 already removes.
4. Ordinary flip failures: X error, T1 decay during the final readout, misassignment. These are the same in both arms, cancel in F2, and are
   cut symmetrically in L2.
5. A coupling that is not by cz wiring: neighbours are cz-wired measured qubits only, as in run_iq.py.
6. Leakage of qubits not measured (spectators). This run does not test that.

## Synthetic test (synth_flip.py, results in synth_results_full.txt)
Everything is in the real raw-IQ schema: same pubs, register names, shapes and complex64 dtype as run_flip.py writes. The roles, cz wiring and
every qubit's 0/1 cluster centres and spreads are the real ibm_kingston ones from 21 Oct. Common to all scenarios: 0.5 % visible readout
excitation, 1.5 % baseline flip failure in both arms. The synthetic baseline F1 came out 4.1-4.3 %; the real kingston calibration-1 figure is 4.1 %.
- (a) planted loss 0.008, no hidden leak (1 seed): EXCLUDED at N = 4 and 8.
  L1 +0.0076 / +0.0091 ABOVE; L2 +0.0076 / +0.0100 ABOVE; L3 0.0000 / +0.0001 CONSISTENT; F2 CONSISTENT.
- (b) loss caused only by hidden-leaked neighbours (6 seeds x 2 N = 12 readings). Each measured qubit leaks at its first mid readout in 1 % of
  readout shots, its IQ point is then 1 sigma off mu0 (inside the 0 cluster) on every later read and after the X, and each wired spectator loses
  0.8 coherence in that shot.
  - L1 ABOVE in 10 of 12 (2 copies-disagree); medians +0.004 to +0.011.
  - L2 CONSISTENT in 10 of 12 (2 copies-disagree); medians -0.0035 to +0.0025; never ABOVE.
  - L3 BELOW in 12 of 12; -0.0060 to -0.0080, half-width 0.0007.
  - F2 ABOVE in 12 of 12; +0.82 to +0.99 percentage points, against about 0.96 planted net of visible excitation.
  - Decision: MAY in 8 of 12, CANNOT_TELL in 4 (a copy split on L1 or L2), EXCLUDED in 0.
- (c) null, 6 seeds x 2 N = 12 readings:
  - L2 ABOVE 0 of 12, so no false EXCLUDED.
  - L1 BELOW 2 of 12 (null2 N = 4, -0.0022; null3 N = 8); L2 BELOW 1 of 12 (null2 N = 4, -0.0023); copies-disagree 2 of 12 each.
  - L3 nonzero 0 of 12 (2 copies-disagree).
  - F2 ABOVE 1 of 12 (null3 N = 4, +0.08 points, at the floor); 1 copies-disagree.
  - F2u is ABOVE in every null, as it should be, because 0.5 % visible excitation is planted.
  - Decision branches fired: 0 of 12.
- Known-good checks on the synthetic data (all 13 runs):
  - CAL 0X flip success: median 97 %, PASS 13 of 13.
  - Idle arm at the calibration floor: PASS 39 of 39.
  - L1w within 21 Oct's interval: PASS 26 of 26 in (a) and (b), where a loss is planted. In the null runs it FAILs in 10 of 12 and passes in
    2, which is correct: no loss is planted there, so failing to reproduce a real loss is the right reading.

## Floor (null half-widths, 3000 shots x 2 copies, kingston geometry)
- L1 and L2: about +/-0.003. A loss at or above about 0.004-0.005 reads ABOVE. The low end of claim 6 (0.004) is at the edge, so L3 carries
  the weight there.
- L3: about +/-0.0005.
- F2: about +/-0.08 percentage points (a hidden-leak rate of about 0.1 % per measured qubit per shot is resolvable).

## Cost (dry-run, 22 Oct about 08:07, nothing submitted; dryrun_flip.log)
| processor | pubs | circuits | ESTIMATED |
|---|---|---|---|
| ibm_kingston | 26 | 2.2-19.7 us | 23.4 s |
| ibm_marrakesh | 26 | 2.7-24.3 us | 23.6 s |
| ibm_fez | 26 | 1.7-15.0 us | 23.2 s |

The estimate model is the one from run_iq.py: shots x (duration + 290 us). All three total 70.2 s against 797 s remaining at the dry-run, so
ibm_fez is included. Raw IQ is about 132 MB per processor, counted at 5.43 bytes per value as on 21 Oct (that run was 90 MB for 64000 shots).

## Addendum before firing (22 Oct, after reconciling two flip-failure counts; no data from this run exists yet)
The flip-failure rule, named: a calibration-1 shot FAILS if it is inside the 0 cluster OR not inside the 1 cluster (inside = Mahalanobis < 3.5;
inside BOTH counts as a failure). Per-qubit rates on 21 Oct data, both rules side by side: derive_flip_guard.py, derive_flip_guard.txt.
Under this rule, on 21 Oct data, marrakesh q82 (98.9 %), q100 (76.8 %), q119 (71.4 %) and fez q149 (24.0 %) exceed 20 %. The committed code
does NOT use that list: it drops, with their wired untouched qubits, whichever measured qubits exceed 20 % on THIS run's own 'in 0 then X'
calibration, so the list on the day may differ; it is printed. Re-reading the 21 Oct marrakesh clean-shot loss without the 5 untouched
qubits wired to q82/q100/q119 (done after the data, 22 Oct, before this run): N=2 +0.0063 [+0.0018, +0.0108], N=4 copies disagree,
N=8 +0.0166 [+0.0112, +0.0259] (81 of 86 qubits). q82 and q100 fail mostly because their 0 and 1 clusters OVERLAP (94.5 % and 74.0 % of calibration-1 shots inside both):
on those qubits "verified in its 0 cluster" cannot tell 0 from 1. The 21 Oct marrakesh clean-shot selection used that verification, so this
run also prints, for marrakesh, L1 and L2 with and without the untouched qubits wired to q82, q100 and q119 (both printed; neither replaces
the committed reading). The L1 known-good comparison for marrakesh is read on the version WITHOUT those qubits as well as with.
