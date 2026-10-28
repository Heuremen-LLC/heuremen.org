# Deciding run: is the extra readout time real duration, a fixed per-measurement delay, or a readout-induced frequency shift?
Written and SHA-256 hashed BEFORE any job is submitted. Heuremen LLC, 21 Oct 2026.

## Design (per chip: ibm_fez L=1.66 us, ibm_marrakesh L=2.684 us, ibm_kingston L=2.18 us; order fez, marrakesh, kingston)
- Emitters E: the 39 qubits in data/reads/countwide_roles_<chip>.json "emitters". Spectators S: every other qubit (117).
- Every spectator: sx (Ramsey start) ... N steps ... quadrature (X: sx | Y: rz(pi/2), sx) ... measure.
- Step, by MODE:
  - D  (reference): emitters idle for exactly L (delay = published measure duration), no measurement.
  - M0: emitters measured (held in |0>).
  - M1: emitters prepared once in |1> (x at start), then measured (QND keeps them mostly in |1>).
  A barrier on all qubits after every step, so spectators see identical scheduling in all modes.
- N in {0, 2, 4, 8, 16}. Pubs per chip: 5 N x 3 modes x 2 quadratures = 30. Shots 2000. One job per chip.
## Analysis (fixed now)
- Per spectator, per mode: phase(N) = atan2(2*P1_Y - 1, 2*P1_X - 1), unwrapped in increasing N; slope s = least-squares phase per step.
- Ratios R0 = s_M0 / s_D, R1 = s_M1 / s_D. Spectator kept if |s_D| >= 1.0 deg/step (declared now; results also printed with NO cut).
- Per chip: median R0, median R1, median (R1 - R0); 68% and 95% intervals from 2000 bootstrap resamples of spectators.
  (One job per chip: chips are the independent units; within-chip intervals are stated as within-job.)
- Extra time per measurement dT = (R0 - 1) * L, per chip.
## Readings, each in one of four states (ABOVE / BELOW / CONSISTENT / CANNOT_TELL), never two
1. R0 vs 1 (measurement vs an equal-length idle, same job, same spectators).
2. R1 - R0 vs 0 (an emitter-state-dependent frequency shift shows here; pure extra time does not).
3. Across the three readout lengths: dT constant in ns (fixed per-measurement delay) vs R0 constant (stretch proportional to L) vs neither/both (CANNOT_TELL).
## What it can NOT separate (stated before data)
- A frequency shift of the spectator that does not depend on the emitter's state remains degenerate with extra time.
