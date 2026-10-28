# Pair run REPEAT: does the entanglement cost of claim 7 come back on new jobs, read exactly as on 21 October?
Written and SHA-256 hashed (PREREG_PAIRS.sha256, first line) BEFORE any repeat job is submitted. Heuremen LLC, Wayfinder and Bones, 22 October 2026.
Original: control_run/PREREG_PAIRS.md (sha256 3270c4b2...133b), control_run/results_pairs.txt (21 Oct). Nothing in control_run is modified.

## Why this repeat
Claim 7 (heuremen.org): Bell pairs among untouched qubits keep 0.007-0.010 less concurrence after 4 readouts than after the same idle, on fez,
kingston and marrakesh, surviving clean-shot selection (every measured neighbour verified in its 0 IQ cluster) and 39-vs-1 readouts on the same
schedule. It rests on ONE run per processor (two jobs, copy 0 and copy 1, minutes apart in one session, 21 Oct 21:15-21:32 EDT). A repeat on a
new day, new calibration and new jobs asks whether it comes back.

## What is identical
- Runner: pairs_repeat/run_pairs_iq.py is a BYTE-IDENTICAL copy of control_run/run_pairs_iq.py (sha256 b4c08604...e828, checked with shasum
  and cmp). It is executed by pairs_repeat/repeat_pairs.py, which only puts control_run on sys.path (so 'import run as base' loads
  control_run/run.py unchanged, sha256 88fd7f71...bd96e) and runs the copy with the same argv. Because the copy lives here, its HERE is
  pairs_repeat: every output (pairsjob_*, pairsraw_*, pairs_backend_properties_*) is written here, and the prereg hash it stamps into each job
  record is the first line of pairs_repeat/PREREG_PAIRS.sha256, i.e. THIS plan's hash.
- Circuits builder (run.py spectator_circuit, pairs_for), roles (countwide_roles_<chip>.json), 48 pairs, 21 singles, K1 = [0] (the same
  single measured qubit), L per chip, settings ZZ XX YY XY YX, shots (CAL 2000, P 4000), arms (idle N = 0; idle and readout of all 39 at
  N = 2, 4; readout of K1 only at N = 4), meas_type kerneled, Darkwave Disco instance. The pub shuffle seeds (20261023 + copy) are part of the
  byte-identical runner, so the pub order inside each copy is the SAME as 21 Oct's.
- Checked 22 Oct 08:28 EDT (read-only, no job): today's build gives E, S, pairs and singles identical to the 21 Oct job records on all three
  processors, and the cz edge set in today's properties equals the 21 Oct saved properties (176 edges each, no difference).
- Analysis: control_run/analyze_pairs_iq.py UNCHANGED (sha256 74e846a9...99c3), run on this folder:
  `venv_quantum/bin/python control_run/analyze_pairs_iq.py pairs_repeat > pairs_repeat/results_pairs_repeat.txt`

## What changed
New jobs (new day, new calibration, IBM's properties saved again at firing), new folder (pairs_repeat), the launcher, this plan and its hash.
Firing order as 21 Oct: ibm_kingston, ibm_marrakesh, ibm_fez: `venv_quantum/bin/python pairs_repeat/repeat_pairs.py --chip <chip> --fire`.
The runner skips a copy when remaining seconds < estimate + 3 (unchanged rule).

## Cost (dry run 22 Oct 08:28 EDT, no --fire, pairs_repeat/dryrun_pairs.log)
ibm_fez 37 s per copy, 73 s per processor; ibm_kingston 37 s / 74 s; ibm_marrakesh 37 s / 74 s; total ESTIMATED 221 s for six jobs.
21 Oct charged 35-37 s per job (fire_pairs.log: kingston 35 + 35, marrakesh 37 + 37, fez 35 + the last job's reading was overtaken by an
allowance refresh). Remaining at the dry run: 578 s.

## Readings (exactly as the original's R1-R5; four outcomes; 95 % bootstrap over pairs, 2000 resamples; copy 0, copy 1 and pooled; copies that disagree -> CANNOT_TELL)
Concurrence: X-state form from the five settings. Neighbours of a pair = measured qubits wired to either member (cz in IBM's calibration
saved at firing). Pairs whose idle concurrence is 0 at that N are excluded (count printed).
R1 per processor, N = 2, 4: median over pairs with at least one wired measured qubit of
   concurrence(idle, all shots) - concurrence(K = 39 readout, only shots where every neighbour was verified 0 on every readout).
   ABOVE 0 -> the entanglement cost survives with every neighbour verified in 0: a leaked neighbour is EXCLUDED as its cause.
R2 per processor, N = 4: median over pairs of loss(K = 39) - loss(K = 1), each loss = idle concurrence minus readout concurrence (all shots).
   Same schedule in both, so ABOVE 0 -> extra elapsed time per step (and anything tied to the schedule) is EXCLUDED as the cause.
R3 the same as R1 with all shots, for comparison. R4 fraction of shots kept by the R1 selection. R5 fraction of mid readouts outside both clusters.
(Classification identical to PREREG_IQ.md: mu0, mu1 per qubit from calibration; projection onto mu1 - mu0, bit = 1 above 0.5; a mid readout is
verified 0 if its Mahalanobis distance to the 0 cluster is below 3.5.)

## Decision
The repeat SUPPORTS claim 7 on a processor if R1 N=4 FINAL is ABOVE. FINAL CONSISTENT, BELOW or CANNOT_TELL: the repeat does not support
claim 7 on that processor (CANNOT_TELL is its own bin, not a refutation). R2 N=4 FINAL is read the same way for the 39-vs-1 half of the claim.
The repeat is reported SIDE BY SIDE with 21 Oct, per processor and per reading, and NEVER pooled with it (same qubits, same pairs, same circuits).
Route, printed beside each FINAL ABOVE (decision unchanged): "both copies ABOVE" or "copies both CONSISTENT, pooled ABOVE". In the synthetic
nulls below every false FINAL ABOVE came by the second route.

KNOWN-GOOD: on real hardware the idle arm at N = 0 (no readouts, no idle steps) must show the wired pairs entangled: pairs_repeat/known_good_pairs.py
reads the median N = 0 concurrence over wired pairs per processor and copy with the analysis's own ruler (95 % bootstrap, 2000, four
outcomes); known-good holds if ABOVE 0, and the value is printed beside 21 Oct's (known_good_20261021.txt). No new threshold. A reading that is not ABOVE 0 voids
that processor's repeat readings (reported, not read). The 21 Oct resolved readings below are the reference the repeat is read against.

## 21 Oct reference (control_run/results_pairs.txt, FINAL; pooled median)
| processor | R1 N=4 (clean) | R2 N=4 (39 vs 1) | R3 N=4 (all shots) | R1 N=2 | R3 N=2 |
|---|---|---|---|---|---|
| ibm_fez | ABOVE +0.0070 (both copies ABOVE) | CANNOT_TELL (copy 0 CONSISTENT) | ABOVE +0.0093 | CANNOT_TELL | CANNOT_TELL |
| ibm_kingston | CANNOT_TELL (copy 0 CONSISTENT +0.0098) | ABOVE +0.0104 | ABOVE +0.0106 | CANNOT_TELL | CANNOT_TELL |
| ibm_marrakesh | ABOVE +0.0093 (both copies ABOVE) | ABOVE +0.0088 | ABOVE +0.0124 | CANNOT_TELL | ABOVE +0.0089 |
Note for the referee: on 21 Oct the clean-shot reading R1 N=4 resolved ABOVE on fez and marrakesh, CANNOT_TELL on kingston; the 39-vs-1 reading
R2 resolved ABOVE on kingston and marrakesh, CANNOT_TELL on fez. Only marrakesh had both.

## SIDE reading, ibm_marrakesh (printed beside, does not replace)
Measured qubits q82, q100, q119 on marrakesh have 0/1 IQ clusters that overlap or miss (flipcheck_run/derive_flip_guard.txt, 21 Oct calibration-1:
flip-rule failure 98.9 %, 76.8 %, 71.4 %; q119 lands outside the 1 cluster in 45.9 % of prepared-1 shots). "Verified 0" is weak for them.
pairs_repeat/side_marrakesh.py runs the UNCHANGED analysis on the marrakesh records with the pairs wired (cz) to 82/100/119 removed (their fin
columns dropped, every other register as is). On 21 Oct this removes 5 pairs: [76,81] [83,84] [101,102] [112,113] [132,133]
(side_marrakesh_20261021.txt: R1 N=4 FINAL ABOVE +0.0093, n=38; R2 FINAL ABOVE +0.0088; medians unchanged because the removed pairs split
two above, two below; two of them, [112,113] and [132,133], lose 0.26 and 0.23 concurrence all-shot at N = 4, copy 0).
R5 in the side output still counts all 39 measured qubits.

## Synthetic tests
- The ORIGINAL 21 Oct test exists: its generator (make.py) and its synthetic pairsjob/pairsraw/properties files were written 21 Oct
  21:15:16-18 EDT, 13 s before the 21 Oct plan commit (21:15:29). Its analysis output was not saved and the generator was not published.
  Both are preserved here unchanged (pairs_repeat/synth_orig_20261021/, file times kept) and the unchanged analysis is RE-RUN on the same
  files, output labelled re-run 22 Oct: synth_orig_rerun_20261022.txt. make.py plants a toy (fake 156-qubit chip, 800 shots, loss 0.03,
  leaked neighbours off-cluster with concurrence 0 in those shots, no null). Re-run: R1 N=4 FINAL ABOVE +0.0260 (planted 0.03), R1 N=2 ABOVE
  +0.0278, R3 N=4 +0.0446 (carries the leaked shots), R2 ABOVE +0.0448, R4 kept 93.0 % at N = 4.
- NEW for the repeat: pairs_repeat/synth_pairs.py, the REAL schema at the real scale (real 21 Oct roles, pairs, singles, K1, pub order and meta
  per copy, real cz wiring, every qubit's IQ centre and spread from the real 21 Oct calibration pubs; 4000 shots; D-mode mid registers zeros as
  on hardware; 0.25 % visible excitation per measured qubit per mid readout). Scenarios: a = loss 0.009 (N=4, K=39), 0.0045 (N=2), 0.001
  (K=1) in every shot of wired pairs plus a 0.5 kick in shots with an excited wired neighbour; leak = only the kick; null = nothing (six seeds).
  Results (three processors each; synth_*.txt):
  a (3 seeds, 9 processor runs): R1 N=4 FINAL ABOVE 9/9, pooled 0.0077-0.0115 (planted 0.009 x readout contrast); R2 ABOVE 9/9 (0.014-0.017,
  includes kicks); R1 N=2 ABOVE 7/9, CANNOT_TELL 2/9.
  leak (2 seeds, 6 runs): R1 N=4 FINAL never ABOVE (CONSISTENT 4, CANNOT_TELL 2); R3 N=4 ABOVE 6/6: the clean-shot selection removes a loss
  carried only by visibly excited neighbours.
  null (6 seeds, 18 runs): R1 N=4 FINAL ABOVE 0/18 (pooled medians -0.0035 to +0.0027). False FINAL ABOVE elsewhere: R1 N=2 3/18, R3 N=2 2/18,
  R2 2/18, R3 N=4 0/18; every one by the route "copies both CONSISTENT, pooled ABOVE", pooled medians +0.0013 to +0.0025.
  Floor: null pooled R1 N=4 within +-0.0035; the 21 Oct R1 N=4 pooled values (+0.0070 to +0.0100) sit above it.
  What the synthetic test cannot see: a hidden leak that reads INSIDE the 0 cluster (flipcheck_run addresses that), drift between copies.

## Files and hashes
CODE_PAIRS.sha256: run_pairs_iq.py (copy), repeat_pairs.py, control_run/analyze_pairs_iq.py (original path), control_run/run.py, synth_pairs.py,
side_marrakesh.py, known_good_pairs.py, synth_orig_20261021/make.py.
