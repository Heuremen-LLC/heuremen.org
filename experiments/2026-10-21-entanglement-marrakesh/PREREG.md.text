# Entanglement run: does a mid-circuit readout destroy entanglement in untouched qubits beyond an equal idle?
Written and SHA-256 hashed BEFORE any job is submitted. Heuremen LLC, 21 Oct 2026.

## Why
The 21 Oct deciding run showed untouched qubits turn 17-20 % further in phase during a readout than during an equal idle.
A phase shift that is a plain frequency shift is a LOCAL UNITARY and cannot change entanglement. Extra decoherence can.
So concurrence of Bell pairs among the untouched qubits separates "readout shifts frequency" from "readout destroys entanglement".

## Design (one job per chip; primary chip ibm_marrakesh, L = published measure duration; replicate ibm_kingston if time allows)
- Emitters E: the 39 qubits in countwide_roles_<chip>.json "emitters" (same as the deciding run). Spectators S: all others.
- Pairs: greedy matching over coupling-map edges with both ends in S, in sorted edge order; each spectator in at most one pair.
- Every pair (a,b): h(a); cx(a,b) -> |Phi+>. Barrier on all qubits.
- N steps of MODE, barrier after each step:
  - D (reference): emitters idle for exactly L. No measurement.
  - M0: emitters measured (held in |0>).
- N in {0, 8, 16}; N=0 is shared (no steps). Pubs: D at 0, 8, 16; M0 at 8, 16 -> 5 conditions.
- Pair readout settings (same on every pair): ZZ, XX, YY, XY, YX (X: h; Y: sdg,h). 5 conditions x 5 settings = 25 pubs.
- 1200 shots per pub. Free open-plan account first; the paid instance only for the replicate.

## Analysis (fixed now)
- Per pair, per condition: from ZZ the four populations p00,p01,p10,p11; from the X/Y settings <XX>,<YY>,<XY>,<YX>.
  rho03 = (XX - YY - i(XY + YX))/4,  rho12 = (XX + YY + i(XY - YX))/4.
  Concurrence (X-state form): C = 2*max(0, |rho03| - sqrt(p01 p10), |rho12| - sqrt(p00 p11)).
  (|Phi+> under amplitude damping, dephasing and local Z rotations stays an X state; local Z rotations leave C unchanged.)
- Per pair: loss per step k_mode = (C(0) - C(16 of mode)) / 16, and also the least-squares slope over N = 0, 8, 16.
- Ratio Q = k_M0 / k_D per pair. Keep pairs with C(0) >= 0.5 and k_D >= 0.002 per step (declared now; results also printed with NO cut).
- Per chip: median Q, 68 % and 95 % intervals from 2000 bootstrap resamples of pairs (one job: within-job intervals).
- Also: median over pairs of C_M0(16) - C_D(16), bootstrap 95 %.

## Readings, each in one of four states (ABOVE / BELOW / CONSISTENT / CANNOT_TELL), never two
1. Q vs 1: does a readout cost a neighbouring Bell pair more entanglement than an equal-length idle?
   - Q CONSISTENT with 1: the deciding run's phase excess is a frequency shift (unitary), not entanglement loss.
   - Q ABOVE 1 and its interval includes ~1.2: consistent with extra elapsed time (cannot separate from dephasing at that size).
   - Q ABOVE 1 and its interval excludes 1.2 from above: readout destroys entanglement beyond time (measurement-induced dephasing).
   - Interval too wide to exclude both 1 and 1.2: CANNOT_TELL.
2. C_M0(16) - C_D(16) vs 0.
Exploratory only (labelled as such): pairs adjacent to an emitter vs pairs not adjacent.
