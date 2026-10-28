# Delivered Hahn-echo T2 against IBM's stated T2, same hour
Heuremen LLC, Wayfinder and Bones, 22 October 2026 (night). Written, SHA-256 hashed (PLAN_HAHN.sha256; code CODE_HAHN.sha256) and pushed
to heuremen.org BEFORE any job of this run is submitted.

## Why
QIA-2026-001 holds a fez line as CANNOT_TELL: untouched-qubit coherence on ibm_fez sits below what IBM's stated T2 and stated zz predict
(with an echo, -0.095 [-0.144, -0.034] at 16 readouts). The model is built from IBM's numbers (not IBM's number). The deciding test asks
IBM's figure directly: IBM's stated T2 is a Hahn-echo T2, so measure the Hahn-echo T2 of the same qubits in the same hour and compare.

## Design (hahn.py), ibm_fez, ibm_kingston, ibm_marrakesh, open instance of .ibm_quantum_key_3
Two sparse sets per chip (no two qubits of a set coupled; every other qubit idle in |0>, so static ZZ to idle neighbours is constant and the
echo refocuses it). Per measured qubit: SX, wait tau/2, X, wait tau/2, SX-dagger, measure. tau = 0, 10, 20, 40, 60, 90, 120, 160, 220, 300,
400 us. Each circuit twice (copy 0 / copy 1), shuffled (seed 20261024), 2000 shots. IBM's calibration record saved at firing.

## Readings (hahn.py --analyze CHIP)
Per qubit, P(1) = 1/2 + A exp(-tau/T2) fitted by least squares (qubits with A > 0.2 kept); ratio delivered / stated T2; median over qubits
with a 4000-resample bootstrap; ABOVE / BELOW / CONSISTENT against 1; copies compared by the median per-qubit gap (an interval excluding 0
-> CANNOT_TELL). Count of qubits below their stated T2.

KNOWN-GOOD: the median P(1) at tau = 0 must be at least 0.90 (the echo sequence and readout work); if not met, that chip is not decided.

## Synthetic test (hahn.py --synth WORLD -> synth_hahn.text), run before this plan was hashed
match (true = stated): 0.999 [0.992, 1.007] CONSISTENT; short (true = 0.7 x stated): 0.700 [0.695, 0.703] BELOW; mixed (half the qubits at
0.7): 0.845 [0.713, 0.985] BELOW. Copies agree in every world.

## What this run cannot do
The fit is single-exponential on a 0 to 400 us window; IBM's own fit window and model are not published, so a ratio below 1 on qubits whose
decay is not exponential reads the shape as well as the time. Qubits with stated T2 well above 400 us are fitted on the early decay only.
