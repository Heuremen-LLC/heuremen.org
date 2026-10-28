# Excitation test: is the early extra coherence loss of untouched qubits during readouts explained by measured qubits accidentally left in |1>?
Written and SHA-256 hashed BEFORE any job is submitted. Heuremen LLC, 21 Oct 2026.

## Design (one job per chip; ibm_kingston then ibm_marrakesh, then ibm_fez if budget remains)
- Measured qubits E: the 39 in countwide_roles_<chip>.json "emitters". Spectators S: all other 117, each prepared on the equator (sx).
- Four modes, N steps each, barrier on all qubits after every step:
  D0: E idle exactly L (published measure duration), E in |0>.      M0: E measured each step, E prepared |0>.
  D1: E prepared |1> (x), then idle exactly L per step.               M1: E prepared |1>, measured each step.
- N in {0, 4, 8}. Spectators read in X or Y. Every mid-circuit readout of E is stored in its own classical bits (one register per step).
- 4 modes x 3 N x 2 quadratures = 24 pubs, 1500 shots.
## Quantities (fixed now)
- Spectator coherence c = |<X> + i<Y>|. Loss L_mode(N) = c(D0,0) - c(mode,N).
- Readout excess X = L_M0(8) - L_D0(8). Excitation cost Y = L_D1(8) - L_D0(8) (an E fully in |1> while idle).
- p = fraction of 1 outcomes in M0's mid-circuit readouts of E, averaged over steps; per spectator, p_nb = mean p over its directly coupled E (coupling edges from the roles file); spectators with no coupled E use the chip mean p.
- Mundane prediction per spectator: X_pred = p_nb * Y.
## Readings (four states: ABOVE / BELOW / CONSISTENT / CANNOT_TELL)
1. Median over spectators of (X - X_pred), 2000 bootstrap resamples, per chip:
   ABOVE 0 -> the readout excess is larger than accidental excitation can explain (mundane cause refuted on that chip);
   CONSISTENT with 0 -> accidental excitation is sufficient; BELOW 0 -> it over-predicts.
   If median X itself is CONSISTENT with 0 the test CANNOT_TELL (nothing to explain).
2. Same split: spectators directly coupled to an E vs not coupled (reported both).
3. Median X and median Y with intervals, and p per chip, printed.
