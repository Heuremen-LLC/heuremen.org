# Headline run: does a mid-circuit readout destroy superposition AND entanglement in qubits it never touched, beyond an equal idle?
Written and SHA-256 hashed BEFORE any job is submitted. Heuremen LLC, 21 Oct 2026. Chips: ibm_marrakesh, ibm_kingston, ibm_fez (one job each).

## Design (per chip)
- Emitters E: the 39 qubits in countwide_roles_<chip>.json "emitters". Spectators S: all others.
- Pairs: greedy matching over coupling-map edges with both ends in S, sorted edge order. Paired spectators: h(a); cx(a,b) -> |Phi+>.
  Unpaired spectators ("singles"): sx -> superposition on the equator.
- Barrier on all qubits, then N steps of MODE, barrier after each:
  D: emitters idle for exactly L (the published measure duration). M0: emitters measured (held in |0>).
- N in {0, 8, 16, 32}; N=0 shared. Conditions: D0, D8, D16, D32, M8, M16, M32 (7).
- Readout settings (applied to every pair (first,second) and to every single using the FIRST letter): ZZ, XX, YY, XY, YX
  (X: h; Y: sdg,h; for singles, the sx preparation is followed by the basis rotation of the first letter; Z gives populations).
  7 conditions x 5 settings = 35 pubs, 2000 shots each.
## Analysis (fixed now)
ENTANGLEMENT, per pair per condition: populations p00..p11 from ZZ; <XX>,<YY>,<XY>,<YX>;
  rho03 = (XX-YY-i(XY+YX))/4, rho12 = (XX+YY+i(XY-YX))/4, C = 2 max(0, |rho03|-sqrt(p01 p10), |rho12|-sqrt(p00 p11)).
  Local frequency shifts (local Z rotations) leave C unchanged.
SUPERPOSITION, per single per condition: coherence c = |<X> + i<Y>| from the XX/XY (X) and YY/YX (Y) settings (pooled over the two settings each).
  A frequency shift rotates the phase and leaves c unchanged.
Per qubit/pair and mode: loss L_mode(N) = value(0) - value(N). Ratio Q(N) = L_M(N) / L_D(N).
Primary readings at N = 32 (and printed at 8, 16): keep units with value(0) >= 0.5 and L_D(32) >= 0.02 (declared; also printed with no cut).
Per chip: median Q, 68 % / 95 % intervals from 2000 bootstrap resamples of units. Also median of value_M(32) - value_D(32), bootstrap 95 %.
## Readings (four states: ABOVE / BELOW / CONSISTENT / CANNOT_TELL)
1. Entanglement Q(32) vs 1.   2. Superposition Q(32) vs 1.   3. Both, on each of three chips.
- ABOVE 1 on a chip: on that chip a readout destroys more of the untouched qubits' entanglement / superposition than an equal idle.
- CONSISTENT with 1: it costs what the time costs; the phase excess is a shift.
- Headline holds if ABOVE on at least 2 of 3 chips for each quantity and BELOW on none; otherwise reported as it falls.
Exploratory (labelled): units adjacent to an emitter vs not adjacent.
