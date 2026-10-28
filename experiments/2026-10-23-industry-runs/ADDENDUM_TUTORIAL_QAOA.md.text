# Addendum to PLAN_TUTORIAL_QAOA.md (23 Oct 2026, /q pass; both runs in flight, NO result read)
1. UNEQUAL APPARATUS (check 9): the hardware runs with IBM's dynamical decoupling (XY4) and gate twirling; the forecast simulation has
   neither, and both normally help the hardware. So R2/R3 BELOW (hardware worse than forecast) is robust to this difference; R2/R3
   CONSISTENT or ABOVE cannot be credited to the forecast being pessimistic or optimistic, and is printed as such.
2. RECORD (check 3): forecast_final builds its noise model from IBM's record fetched at the END of the run; the record at the START is
   saved (props_tqaoa_<chip>.json). A reading counts only if every gate error, readout error, T1 and T2 on the compiled layout's qubits and
   edges is identical in the two records; otherwise the forecast is recomputed from the start record and both are printed.
3. Wording (check 4): the forecast is 'Qiskit Aer's noise model built from IBM's published figures', never 'IBM's model'.
