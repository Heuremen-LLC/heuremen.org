# Deciding runs: six readings already on the record, each read again beside the ordinary causes that could make it
Written, SHA-256 hashed and pushed to GitHub BEFORE any job is submitted. Version 3: the across-chip table combines chips by VALUE and WIDTH (not by counting labels) and prints the chip-to-chip spread;
the floor re-measured on REAL hardware copies. Version 2 (0ea33d5) added the measured FLOOR, the trimmed mean for R4.2/R4.3 and
the run 4 wording (version 1, commit 63482af, never fired). Heuremen LLC (Wayfinder and Bones), 22 October 2026.
Circuits: runs.py (dry-run against the live chips, nothing submitted). Analysis: analyze.py, tested on synthetic data with planted effects
(synth.py, output synth/synth_output.txt) and committed with this plan.

## The ruler, for every reading
Every reading is a GAP between two arms of the same job, per qubit, median over qubits, 95 % bootstrap (2000 resamples), one of four outcomes:
ABOVE (interval above 0), BELOW (below 0), CONSISTENT (contains 0; the interval is printed as the bound), CANNOT TELL (fewer than 5 qubits, or
the copies shifted). Every circuit runs twice in its job (copy 0, copy 1, shuffled together). Copies are compared by their difference
(per qubit copy 0 - copy 1): interval contains 0 -> the pooled value is the reading; excludes 0 -> CANNOT TELL with both copies printed.
On synthetic data with identical copies this flag fired on 2 of 52 readings (3.8 %); it can only turn a reading into CANNOT TELL.
No sigma thresholds. No outcome is decided by comparing two labels. Signs are carried: a size is the kick times the sign of its own median.

## The runs (chips; seconds from the dry run)
1 DISTANCE (fez, kingston, marrakesh; 29 s each). 39 measured qubits packed in one region (breadth-first ball from qubit 0), 117 untouched at
  distance 1 to the far side; the control run's spread layout in the same job. Idle vs readout, N = 8, untouched qubits read in X and Y.
  R1.1 kick; R1.2 kick size at distance <= 2, 3-4, >= 5; R1.3 size(>= 5) - size(<= 2) (ABOVE larger far, BELOW larger near);
  R1.4 slope of size on the sum over measured qubits of 1/d^2 (the record's stated form); R1.5 size(<= 2) spread minus packed layout.
2 READOUT COUNT (fez, kingston, marrakesh; 17-19 s). Control-run layout, idle vs readout at N = 1, 3, 9, 27, 54, 81.
  R2.1 kick per readout at each N; R2.2 kick per readout over readouts 27->81 minus over 1->9; R2.3 intercept of phase gap on N;
  R2.4 coherence loss at each N; R2.5 SUM form loss(81) - 3 loss(27); R2.6/2.7 PRODUCT form -ln(readout/idle)/N at 9, 27, 81 and 81 minus 9.
  Both forms are printed side by side: coherence may combine per readout as a product, not a sum.
3 READ ORDER (fez, kingston, marrakesh; 28-29 s). All qubits in 0, three seeded groups, three readout slots. All six orders; each group read
  ALONE at each slot (nobody else read); all three together at slot 1. R3.1 P(1) after others minus alone, slots 1-3 (slot 1 is a built-in
  none: nothing is read before it); R3.2 slot 3 minus slot 2 of that; R3.3 alone slot 3 minus alone slot 1 (waiting longer, measured
  directly); R3.4 together minus alone at slot 1; R3.5 last minus first in the full orders.
4 GATES (fez, kingston, marrakesh, phoenix; 8 s each). X trains of 2, 6, 10 steps 30 us apart, each step X or X X X; the waits are the same,
  so the two trains differ in time only by two gate lengths per step (at most about 0.7 us of 300 us). A barrier after every X so X X X is
  not compiled to X (caught in the dry run). References X then wait 30, 90, 300 us.
  R4.1 P(1) X X X minus X per train; R4.2 cost per extra X gate (slope; the decay shape is the same in both trains and cancels);
  R4.3 that minus IBM's stated X error. R4.2 and R4.3 use a 20 % trimmed mean over qubits, not the median (see FLOOR).
  This run does NOT re-read the record's 0.0014-0.0020 shortfall against a decay-only prediction; it reads the quantity that decides whether
  an X gate costs more than stated once the decay shape is held equal.
5 NEIGHBOUR STATE (fez, kingston, marrakesh; 14 s each). Targets at least 3 apart run m = 0, 40, 80, 160 X gates (barrier after each), read in
  Z; their neighbours prepared in 0, 1 or + and never read. R5.1 error per gate per arm; R5.2 1 minus 0 (static neighbour state);
  R5.3 + minus 0 (the record's comparison); R5.4 + minus the mean of 0 and 1 (what a static neighbour state predicts for +); R5.5 ratio.
6 REVIVAL (phoenix and fez; 84 s each). Coupling graph two-coloured: colour A runs free clocks (SX, wait tau, read X and Y), tau = 2 to 150 us
  in steps of 4; colour B in 0, 1 or + and never read. R6.1 come-back = mean coherence 58-78 us minus 26-50 us, per arm; R6.2 + minus 0;
  R6.3 1 minus 0; R6.4 held-out comparison: fit plain decay and damped beat on odd taus, score on even taus, SSE(decay) - SSE(beat).

## FLOOR (measured on the analysis before any job; Guard 4)
Six null synthetic seeds, every run, no effect planted, same shots and analysis: 10 of 210 readings whose truth is 0 came out ABOVE or BELOW
(4.8 %, about half each side). The gate cost R4.2 with a median came out ABOVE in 2 of 6 null seeds: its per-qubit values sit on a grid of
shot counts and the median's bootstrap is too narrow on tied values; with a 20 % trimmed mean, 1 of 20 null seeds (outputs in nulls/).
REAL HARDWARE FLOOR: the same statistics on tonight's control run, copy 0 minus copy 1 of every circuit group (truth 0, including any drift
between circuits inside a job), 147 comparisons on fez, kingston, marrakesh: median 2.0 % resolved (1 above, 2 below), 20 % trimmed mean
4.8 % (4 above, 3 below). The drift inside a job does not widen the floor past the 95 % interval's own rate.
ACROSS CHIPS: for each reading, every chip's value and interval; the inverse-variance combination (sigma = interval width / 3.92) with its
interval and outcome; and the chip-to-chip spread, chi2 of the chips about the combination per dof, printed as its own reading (near 1 when
the chips agree within their widths). Chips are never combined by counting their outcomes.

## Not run
The GHZ flatness ("pinata") metric P(000) + P(111) with tau/3 offsets is exactly constant (1/16) in textbook quantum mechanics: the sum over
pairs of cos(ti)cos(tj) is -3/4 for every theta. Its textbook flatness ratio is infinite; the measured 1521x is the hardware's distance from
an identity. No run.
Total from the dry run: about 468 s. Remaining at writing: 2293 s.
