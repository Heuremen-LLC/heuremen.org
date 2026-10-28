# DD run: does an echo on the untouched qubits remove the readout kick and the readout coherence cost?
Pre-registered. Written, SHA-256 hashed (PLAN_DD.sha256; code in CODE_DD.sha256) BEFORE any job is submitted. Heuremen LLC, Wayfinder and
Bones, 22 October 2026. Circuits and fire path: dd.py (dry run against the live chips, nothing submitted). Analysis: analyze_dd.py, which
imports the deciding runs' ruler (deciding_runs/analyze.py) unchanged. Synthetic test: synth_dd.py (synth_output.txt, null_floor.txt, nulls/).

## What is on the record
On ibm_fez, ibm_kingston and ibm_marrakesh (control run, 21 Oct), mid-circuit readout of 39 measured qubits moved the phase of the 117
untouched qubits by -1.6 / -2.9 / -2.4 deg per readout against an idle of exactly the same scheduled length, and cost them 0.004-0.021 of
coherence |<X> + i<Y>|. The kick also reached untouched qubits with no coupler to a measured qubit.

## The question
Does dynamical decoupling of the untouched qubits during each step remove (a) the phase kick and (b) the coherence cost?
Like for like: the echo is applied identically in the idle arm and in the readout arm, so every reading is still readout minus idle.

## Design (per chip; one job per chip)
Roles of the control run: E = the same 39 measured qubits (countwide_roles_<chip>.json), S = the other 117. Spectators start with SX (equator).
N in {0, 2, 4, 8, 16} steps. Each step: barrier on all 156 qubits; E idle exactly L (arm D) or measured, held in 0 (arm M); during the step
S either waits (DD off: no instruction, the scheduler pads the idle) or runs a symmetric X-X echo that fills exactly L (DD on):
  delay d1, X, delay d2, X, delay d3,   d1 = d3 = (L/4 - g/2) rounded down to the alignment, d2 = L - 2 g - 2 d1, total = L (asserted)
where L = the published measure duration (identical on all 39 measured qubits, asserted) and g = the X duration (identical on all 117
spectators, asserted). Barrier after the step. End: S rotated to X or Y (as run.py), every qubit measured into one final register.
Arms: (D, off), (M, off), (D, on), (M, on) at N = 2, 4, 8, 16; N = 0 built once under each DD label (no step, so the two N = 0 circuits are
IDENTICAL: a built-in null, R0). 36 circuits x 2 copies = 72 pubs, 4000 shots each, shuffled with a seeded RNG (20261033).
Transpile: optimization level 1, identity layout, ALAP scheduling (as every run before).

Timing from the live targets (dry run, 22 Oct; dt = 4 ns, granularity 1, pulse alignment 1, acquire alignment 1 on all three):
  ibm_fez        L 415 dt (1.66 us), X 6 dt -> echo 100 / X / 203 / X / 100 dt
  ibm_kingston   L 545 dt (2.18 us), X 8 dt -> echo 132 / X / 265 / X / 132 dt
  ibm_marrakesh  L 671 dt (2.68 us), X 9 dt -> echo 163 / X / 327 / X / 163 dt
The echo is exactly symmetric about L/2 (d1 = d3); rounding d1 to whole dt puts the X centres 0.75 dt (fez) / 0.25 dt (kingston,
marrakesh) from L/4 and 3L/4, so a constant detuning delta leaves a residual of delta x 1.5 dt (fez) at most per step.
VERIFIED in the dry run on the transpiled, scheduled circuits of all three chips (dd.py verify(), raises on failure):
  every spectator carries exactly 2N X gates with DD on and 0 with DD off, no X on a measured qubit (the delays between the two X keep them
  from being cancelled); every step, full barrier to full barrier, lasts exactly L; the step block and the total scheduled duration are
  identical between DD on and DD off for every (arm, N, setting, copy). Circuit durations (both DD settings):
  fez 1.7 / 5.0 / 8.3 / 15.0 / 28.3 us, kingston 2.2 / 6.6 / 11.0 / 19.7 / 37.1 us, marrakesh 2.8 / 8.1 / 13.5 / 24.2 / 45.7 us (N = 0..16).

## Seconds (dry run, 22 Oct 2026; per-shot model of the deciding runs, shots x (duration + 290 us))
  ibm_fez 87 s, ibm_kingston 88 s, ibm_marrakesh 89 s. Total about 264 s. Remaining at writing: 1732 s.

## The ruler (same as the deciding runs)
Every reading is a GAP between arms of the same job, per spectator, median over spectators, 95 % bootstrap (2000 resamples), one of four
outcomes: ABOVE (interval above 0), BELOW (below 0), CONSISTENT (contains 0; the interval is the bound), CANNOT TELL (fewer than 5
spectators, or the copies shifted). Copies: per spectator copy 0 - copy 1; interval contains 0 -> the pooled value is the reading;
excludes 0 -> CANNOT TELL with both copies printed. Signs carried, never folded. Chips combined only by inverse variance (sigma = width/3.92),
with chi2/dof of the chips about the combination printed as its own reading. Phase convention of the control run: phase = atan2(<Y>, <X>),
gap = wrap(phase_M - phase_D) in degrees, unwrapped along N by adding wrapped increments. Every reading is printed three times: all
spectators, WIRED (cz-coupled in the saved properties to a measured qubit) and UNWIRED (the kick reaches unwired qubits).

## Readings (per chip, then across chips)
R0  Built-in null: N = 0 "DD on" minus "DD off", phase and coherence. Identical circuits; truth 0.
    ABOVE/BELOW: drift inside the job (or the floor) at that size; it bounds every other reading. CONSISTENT: expected.
R1  Kick per readout, DD off: per spectator, least-squares slope (with intercept) of phase_M - phase_D on N = 2, 4, 8, 16; median.
    The record predicts BELOW (about -1.6 to -2.9 deg). CONSISTENT: this job does not reproduce the kick, and R3 cannot be read as removal.
R2  Kick per readout, DD on. Same slope.
    CONSISTENT: no kick resolved with the echo (the interval is the bound on what is left). BELOW/ABOVE: a kick survives with the echo, of that sign.
R3  R2 - R1 per spectator, then median (spectators alive in all four arms). R3s = sign(R1) x (R2 - R1), sign printed.
    R3s BELOW: the echo shrinks the kick (if R2 is also CONSISTENT: removed within the bound). R3s ABOVE: the echo GROWS the kick.
    R3 CONSISTENT: the echo does not change the kick within the interval. Read together with R2: R2 ~ -R1 (R3s ~ -2|R1|) is a sign FLIP.
R4  Coherence cost, DD off: coherence_D(N) - coherence_M(N) at N = 2, 4, 8, 16. The record predicts ABOVE.
R5  The same with DD on. ABOVE: the echo leaves a cost. CONSISTENT: no cost resolved with the echo.
R6  R5 - R4 per N. BELOW: the echo reduces the cost. ABOVE: the echo increases it. CONSISTENT: no change within the interval.
R7  The echo's own effect on the idle arm: coherence_D,on(N) - coherence_D,off(N). ABOVE: the echo protects the idle clocks (refocuses slow
    detuning); BELOW: the echo costs more (two X per step, gate error) than it refocuses. R7 is the context for R5/R6, not a readout effect.
R8  Dead-clock guard: at each N and arm, the number of spectators (wired / unwired) whose coherence (both copies) is below 5 x shot noise
    (5 / sqrt(8000) = 0.056). Phase readings (R1-R3) use only spectators alive in every arm and N their slope uses (R3: alive in all four
    arms). Coherence readings use all spectators.

## Physics note (a guide to the outcomes, not a reading)
Seen by a spectator, the echo multiplies any Z-type phase accrued during the step by s(t) = +1 in the first quarter, -1 in the middle half,
+1 in the last quarter (the second X restores the frame). So the echo cancels a CONSTANT frequency shift over the step exactly, cancels any
detuning whose s-weighted integral over the window is zero, and leaves anything else. A phase applied as an INSTANTANEOUS frame update (or a
shift concentrated in time) survives with sign +1 if it falls in an outer quarter and is reversed (-1) if it falls in the middle half. A
phase applied outside the step (at its fences) is not touched by the echo at all. Hence, with care:
 - R2 CONSISTENT and R3s BELOW: the kick behaves like a refocusable detuning during the window (constant, or balanced across the quarters).
   This would NOT identify the source (measurement-induced Stark shift, a ZZ shift from a measured qubit's state that holds through the
   step, or anything else refocusable); it says only how the kick sits in time.
 - R3 CONSISTENT (R2 ~ R1): the kick is not refocused: compatible with a phase in the outer quarters, at the step fences, or not of Z type.
   It would NOT show the kick is a software frame update; it shows the echo does not reach it.
 - R2 of opposite sign to R1 (R3s near -2|R1|): the kick lands mostly in the middle half of the step.
 - Partial (R3s BELOW, R2 still BELOW): a mixture of time profiles; the fraction is not a mechanism.
 - Two caveats on the comparison itself: all 117 spectators are flipped together, so the coupling BETWEEN spectators (ZZ) is not refocused
   by this echo (it is refocused against the measured qubits, which are not flipped); and the X pulses fall while the readout tones play,
   so a readout effect on the X pulses themselves would appear in R5/R6 and not in R7.
For the coherence cost: R6 BELOW says the echo removes part of the cost, which points to the cost being dephasing by a refocusable
(slow on the scale of L) detuning; R6 CONSISTENT with R5 ABOVE says the cost is not refocused by one echo per step (fast noise, an
excitation of the spectator, or a phase that varies shot to shot within the window). Neither alone names a mechanism.

## FLOOR (measured on the analysis before any job)
Synthetic data in the job's own schema (synth_dd.py, FakeFez circuits, 4000 shots per pub, two copies, per-spectator coherence 0.55-0.95
with 6 dead clocks, per-spectator start phase, per-spectator drift that DIFFERS between DD on and off but is equal in D and M):
 - PLANTED (kick -2.0 deg per readout with DD off, removed with DD on; coherence cost 0.010 at every N in both DD settings; echo with no
   effect of its own; synth_output.txt): R1 -2.004 [-2.027, -1.983] BELOW; R2 +0.015 [-0.018, +0.036] CONSISTENT; R3 +2.022 [+1.968, +2.048]
   ABOVE (R3s BELOW); R4 0.009-0.012 ABOVE at every N; R5 0.010-0.011 ABOVE at every N; R6 and R7 CONSISTENT at every N for all spectators;
   R0 CONSISTENT; R8 found the 6 dead clocks at every N and arm (111 used for phase). Recovered. On the 31 UNWIRED spectators alone the
   coherence readings are wide (R5 N=16 CONSISTENT despite the 0.010 plant; 3 copies-shifted CANNOT TELL; R6 N=2 BELOW): unwired-only
   coherence readings carry about a quarter of the spectators and should be read with that width.
 - NULL, six seeds (11-16; null_floor.txt, nulls/): 396 readings whose truth is 0: 20 resolved (5.1 %; 10 ABOVE, 10 BELOW), 364
   CONSISTENT, 12 CANNOT TELL (copies shifted, 3.0 %). Null seed 5 (synth_output.txt): 5 of 66 resolved, 2 CANNOT TELL. The resolved null
   phase readings were at most 0.03 deg (0.08 deg copy shift on unwired); the resolved null coherence readings at most 0.009.
   So about 1 in 20 of the 66 readings per chip will resolve with nothing there; a single resolved reading at the floor's size is not a finding.
Real hardware adds drift inside the job; R0 measures that drift on each chip in this very job.

## Not in this run
No other DD sequence (XY4, CPMG with more pulses), no echo on the measured qubits, no variation of where in the window the X pulses sit.
Those are the follow-ups R3 would point to; none is pre-registered here.
