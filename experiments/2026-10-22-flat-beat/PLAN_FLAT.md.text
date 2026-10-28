# Flat beat on ibm_phoenix: on FRESH pairs, does an untouched qubit beat against one neighbour in superposition at about 12.5 kHz whatever IBM's calibrated coupling, or at the calibrated coupling?
Written, SHA-256 hashed and pushed to GitHub BEFORE any job is submitted. Heuremen LLC (Wayfinder and Bones), 22 October 2026.
Circuits: flat.py (from zz.py of the one-neighbour run). Analysis: analyze_flat.py (uses the one-neighbour run's coupling-name parser as
fixed and disclosed there). Synthetic test: synth_flat.py with the beat planted as the calibration (track), as 12.5 kHz (flat), and absent
(none), two seeds each; outputs synth_*.txt committed with this plan.

## Why
In the one-neighbour run (22 Oct) an outside referee fitted each phoenix pair's beat freely after the data and found it near 12.5 kHz
(interquartile 11-14) for pairs IBM calibrates anywhere from 10 to 57 kHz (rank correlation 0.15), while on ibm_fez the beat tracked the
calibration (log-log slope 1.08). A finding made after the data is a hypothesis. This run states it first and tests it on pairs it was not
found on.

## Design (ibm_phoenix, one job; every circuit twice, copy 0 and copy 1, shuffled together; 2000 shots)
Triples (B1, A, B2) chosen as in zz.py, seeded differently (20261024), and NO qubit that was an A in the one-neighbour run is an A here.
Dry run: 20 triples, 40 pairs. A: SX, wait tau, read X or Y; B1, B2 read in Z at the end. Arms: Z (B1, B2 in 0), P1 (B1 in +), P2 (B2 in +).
tau = 2-162 us in steps of 4 (41 values). Dry-run estimate 368 s (the last phoenix job was charged about half its estimate).

## Readings (unit = one pair; 95 % bootstrap over pairs, 2000; copies compared)
Per pair the beat f is fitted freely: C_P(tau) = C_Z(tau) |cos(pi f tau)|, f = 1-100 kHz in 0.05 kHz steps, on taus where C_Z > 5x its floor.
Q2 (THE DECIDING READING) log-log slope of fitted beat on IBM's calibrated zz across pairs. Prediction FLAT: consistent with 0 and BELOW 1.
   Alternative TRACK: consistent with 1 and ABOVE 0. Synthetic: track 1.019, 1.023 (both ABOVE 0, slope - 1 within +0.03); flat -0.001,
   -0.003 (CONSISTENT with 0, slope - 1 BELOW); none 0.
Q1 median fitted beat. Prediction FLAT: inside the synthetic flat band, 11.75-12.10 kHz (a planted 12.5 kHz reads about 0.45 kHz low through
   this fit, because the neighbour relaxes); TRACK: near the median calibrated value. 'Q1 minus 12.5' is printed but is BELOW even when the
   truth is 12.5, so it is read against the synthetic band, not against 0.
Q3 rank correlation, fitted beat vs calibrated zz. FLAT: consistent with 0. TRACK: near 1 (synthetic 0.9996-0.9998).
Q4 median fitted beat / calibrated zz (level, printed).
Both outcomes are reported as they come. If FLAT holds, it is a pre-registered, out-of-sample statement that IBM's calibrated coupling on
ibm_phoenix does not give the beat an untouched qubit shows against a neighbour in superposition in these circuits. If TRACK holds, the
22 Oct flat beat was a product of fitting after the data and is withdrawn from the data page. Anything between is reported as measured.
