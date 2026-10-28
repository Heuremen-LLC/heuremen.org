# Control run: does a mid-circuit readout disturb qubits it never touches, and do the ordinary causes explain it?
Written, SHA-256 hashed, and the hash pushed to the public repository github.com/claude-wayfinder/Heuremen.org BEFORE any job is
submitted. Heuremen LLC, 21 October 2026. Circuits: run.py in this folder (same hash commit). Analysis: analyze.py, written before
any result file is opened; its hash is committed before the first result is read.

## Design (one job per processor: ibm_kingston, then ibm_marrakesh, then ibm_fez; a processor is skipped if the remaining time is short)
- Roles: the 39 measured qubits E and the 117 untouched spectators S are the `emitters` in countwide_roles_<chip>.json (same as 21 Oct).
- A readout step measures the chosen measured qubits; the matched idle step holds every measured qubit idle for exactly the stated
  readout duration L of that processor. Barriers on all qubits after every step. Untouched qubits get no gates during the steps.
- Every circuit appears twice in the job (copy 0 and copy 1); all circuits are shuffled together with a fixed seed.
- Arms (shots per copy):
  S (1000): all 117 spectators prepared on the equator (sx); idle N=0, and idle (D) / readout (M) at N = 2, 4, 8; spectators read in X and Y.
  C (1000): as S, N = 8, but only 1 (K1) or 10 (K10) of the 39 are read out each step; the rest idle exactly L. K10 = 10 drawn from E with
            seed 20261021; K1 = the first of those.
  P (4000): 48 Bell pairs |Phi+> on coupled spectator pairs (+21 singles); idle N=0, and D / M at N = 2, 4; settings ZZ, XX, YY, XY, YX.
  L (1000): leakage. E prepared |0>, 8 readouts, then X on E and a final readout of E. Reference: E prepared |1>, one readout, X, final readout.
- IBM's stated gate/readout errors and T1/T2 are saved at the moment of firing (backend_properties_<chip>.json).

## Readings (fixed now). Every reading has four outcomes: ABOVE / BELOW / CONSISTENT (with the named value) / CANNOT_TELL,
## from 95 % bootstrap intervals over qubits (or pairs), 2000 resamples. Every reading is computed on copy 0, copy 1 and pooled;
## if copy 0 and copy 1 give different outcomes, the pooled reading is reported as CANNOT_TELL with both copies printed.
Phase per spectator: phi = atan2(<Y>, <X>) from the X and Y settings. Coherence: |<X> + i<Y>|. Concurrence: X-state form from the five settings.

R1 KICK. Per spectator, kick = (phi_M(8) - phi_D(8)) / 8 in degrees per readout (arm S). Median over spectators vs 0.
R2 DOSE. Kick with 1, 10 and 39 measured (arms C and S, N = 8). Two alternatives, both tested:
   PER-QUBIT SUM: kick(K) = kick(39) * K / 39.     GLOBAL: kick(K) = kick(39).
   For K = 1 and K = 10: interval of the median kick vs each prediction. Outcome GLOBAL if it contains kick(39) and excludes K/39 * kick(39);
   PER-QUBIT if the reverse; BOTH (cannot tell) if it contains both; NEITHER if it excludes both.
R3 STRETCH. Rate per spectator s = (phi_D(8) - phi_D(0)) / 8 (deg per step). Regress kick on s using s from copy 0 and kick from copy 1, and
   again with the copies swapped (independent noise in x and y). Slope a: a pure extra-time stretch gives a > 0. Reported vs 0 and vs
   12/11 - 1 = 0.0909 (both distances printed). Spectators precessing backward (s < 0): median kick vs 0.
R4 SUPERPOSITION COST. Per spectator, coherence_M(N) - coherence_D(N) at N = 2, 4, 8 (arm S). Median vs 0 at each N.
   Growth: slope of the per-spectator difference against N, median vs 0.
R5 ENTANGLEMENT COST. Per pair, concurrence_M(N) - concurrence_D(N) at N = 2, 4 (arm P). Median vs 0. A pair whose idle concurrence is 0
   at that N is excluded and the count excluded is printed (a 0-vs-0 pair cannot tell).
R6 LEAKAGE. Per measured qubit: P(final = 1 | last mid readout = 1) in the leak arm minus P(final = 1 | readout = 1) in the reference arm.
   Median over the 39 vs 0. A leaked qubit reads 1 before and after the X; a qubit genuinely in 1 reads 0 after the X.
R7 EXCITATION. Fraction p of 1 outcomes in the mid readouts of the S-arm readout circuits (N = 8), per chip. Printed; used only to bound
   the excitation route: if R4 at N = 8 is ABOVE 0 in size while p is under 1 %, excitation cannot be the whole of it (reported as a bound,
   not a reading).

## What would count against each claim (stated now)
- Kick: R1 CONSISTENT with 0 on a processor.
- Superposition cost: R4 at N = 8 CONSISTENT with 0, or R6 ABOVE 0 by enough that leakage at the measured rate covers it.
- Entanglement cost: R5 CONSISTENT with 0 at both N.
- Any reading whose two copies disagree is CANNOT_TELL, not a result.
