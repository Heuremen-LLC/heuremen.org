# Check C: does the readout cost survive in shots where every measured neighbour read 0 on every readout?
Defined after the control-run results were seen, before this check was computed. Heuremen LLC, 21 Oct 2026. Script check_c_clean_shots.py.
- Untouched qubit (arm S, N = 2, 4, 8) or Bell pair (arm P, N = 2, 4); neighbours = measured qubits wired to it (to either member of a pair),
  wiring from backend_properties_<chip>.json (cz). Copies pooled.
- Readout arm: keep only the shots in which every neighbour read 0 on every one of the N mid-circuit readouts ("clean shots").
  Idle arm: all shots (it has no mid-circuit readouts). Loss = idle value minus clean-shot readout value.
- Reading per processor and N: median over qubits (pairs) of the clean-shot loss, 95 % bootstrap; four outcomes; fraction of shots kept printed.
  ABOVE 0  -> the cost survives with every neighbour reading 0: neighbours stuck in 1 are excluded as its cause. Leakage to a level
             that reads as 0 is NOT excluded by this check.
  CONSISTENT with 0 where the all-shot loss was ABOVE 0 -> the cost is consistent with stuck neighbours.
  Only qubits with at least one wired measured qubit and pairs with idle concurrence > 0 enter.
