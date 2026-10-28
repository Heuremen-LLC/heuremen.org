# One-neighbour revival: does an untouched qubit's coherence come back on the schedule its own coupling to ONE neighbour predicts, and where is it while it is gone?
Written, SHA-256 hashed and pushed to GitHub BEFORE any job is submitted. Heuremen LLC (Wayfinder and Bones), 22 October 2026.
Circuits: zz.py (dry-run against the live chips, nothing submitted). Analysis: analyze_zz.py, tested on synthetic data (synth_zz.py) with
the ZZ scale planted at 0.5 (two seeds) and 1.0, and with no coupling at all (six seeds); outputs synth_*.txt, committed with this plan.

## Why
On ibm_phoenix (deciding runs, 22 Oct) untouched qubits whose neighbours were parked in superposition lost their coherence by about 26 us
and got some of it back at 58-78 us. A referee found the timing and shape follow IBM's calibrated per-edge ZZ scaled by a fitted factor of
about 0.5, while the size of the return was about four times what that model gives. Each qubit there had 3-4 neighbours in superposition,
so their beats did not line up. Here each untouched qubit has exactly one or two DESIGNATED neighbours, and every other neighbour is idle in 0.

## Design (ibm_phoenix and ibm_fez; one job each; every circuit twice, copy 0 and copy 1, shuffled together; 1000 shots)
Layout (zz.py layout(), seeded, from IBM's calibration at firing): triples (B1, A, B2) with B1, B2 coupled to A; A has readout error < 0.05
and T2 > 40 us, B1 and B2 readout error < 0.05; no qubit of a triple is coupled to the A of another triple, and A is coupled to no qubit of
another triple. Dry run: 21 triples on phoenix, 33 on fez. A: SX, wait tau, read in X or Y. B1, B2 read in Z at the end (never mid-circuit).
Arms: Z (B1, B2 in 0); P1 (B1 in +); P2 (B2 in +); P12 (both in +); O1 (B1 in 1); EZ and EP1 (as Z and P1, Hahn echo X on A at tau/2).
tau: phoenix 2-202 us in steps of 5; fez 3-303 us in steps of 7.5 (41 values each). Dry-run estimate 452 s (phoenix) + 511 s (fez).

## Readings (unit = triple; median over triples; 95 % bootstrap, 2000; four outcomes; copies as in the deciding runs)
Coherence C = |<X> + i<Y>| of A. zz = IBM's calibrated zz for the edge (props at firing, GHz -> Hz).
R1 SCALE. Per triple, fit C_P1(tau) = C_Z(tau) |cos(pi s zz tau)| over s in 0.10-2.50 (step 0.0025). Printed: s, s - 1.0, s - 0.5, and the
   same for P2. The two candidates are 1.0 (IBM's zz is the beat frequency) and 0.5 (the referee's fit on the 22 Oct phoenix data).
   Synthetic: planted 0.5 -> 0.4875 [0.48, 0.49] twice; planted 1.0 -> 0.995 [0.985, 1.0]: the fit reads 1-2.5 % LOW because it ignores B1's
   T1. So "s - 0.5 BELOW" by up to about 0.015 is the fit's own bias, not a miss. No coupling at all -> s pinned at the grid floor 0.10.
R2 RETURN HEIGHT (read only if R1 is inside the grid). Mean C_P1/C_Z at the 3 taus nearest the first revival t* = 1/(s zz), minus
   (1 + exp(-t*/T1_B1))/2. That prediction is an UPPER BOUND for static ZZ with B1 relaxing: synthetic ZZ-only gives -0.083 to -0.004.
   ABOVE -> more return than static ZZ with a relaxing neighbour can give: the extra return seen on 22 Oct survives with one neighbour.
   Inside or below the synthetic band -> ZZ accounts for the return with one neighbour.
R3 WHERE IT IS WHILE IT IS GONE. At the dip t_d = 1/(2 s zz): R3a = C_P1/C_Z (unconditional); R3b = the coherence of A computed separately in
   the shots where B1 ended in 0 and in 1, weighted by branch size, divided by the same statistic on Z with the shots split into two seeded
   random halves (bias matched); R3c = R3b - R3a. ABOVE -> the coherence that disappears from A alone is held in its correlation with B1.
   Synthetic ZZ-only: R3a 0.12-0.16, R3b 0.85-0.91, R3c +0.70 to +0.78. No coupling: R3c consistent with 0 in 6 of 6 seeds.
R4 ECHO. EP1/EZ - 1 at t_d. A static ZZ is refocused by the echo; synthetic ZZ-only with B1 relaxing: -0.014 and -0.048 (third seed: copies shifted). Well below that band
   -> the coupling is not static over the window.
R5 PRODUCT OR SUM, both in +. Over taus where C_Z > 5x its noise floor: SSE(C_P12/C_Z = x1 x2) - SSE(1 - C_P12/C_Z = (1-x1) + (1-x2)).
   BELOW -> the two neighbours combine as a product. Synthetic ZZ: -2.5 to -2.9; no coupling: consistent in 5 of 6 seeds (1 above).
R6 B1 in 1 (reference): mean C_O1/C_Z over live taus. Static ZZ shifts the frequency only (ratio near 1, lowered by B1's T1 jumps).
R7 the scale on phoenix minus on fez (value and width).

## Floor
Six seeds with no coupling: 1 of 18 readings that can resolve (R3c, R4, R5 x 6 seeds) came out resolved (R5, seed 16). With no coupling R1 sits at the grid floor, so R2 is not read (its printed ABOVE there is the upper-bound prediction against a flat curve). On one ZZ seed the echo reading's copies shifted (CANNOT TELL), on the others it read -0.014 and -0.048. Ratios
(R1, R2 as measured, R3a, R3b, R6) are LEVELS: their interval is the reading; their ABOVE/BELOW against 0 means nothing and is not read.
