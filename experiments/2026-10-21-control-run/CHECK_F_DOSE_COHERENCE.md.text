# Check F: does the coherence loss depend on how many qubits are read out, at the same schedule?
Defined after the data, before computing. Heuremen LLC, 21 Oct 2026. Script check_f_dose_coherence.py.
- Arms S (39 read out) and C (1 or 10 read out), N = 8, same 117 untouched qubits on the equator, same step schedule (each step is one
  readout slot of the stated length followed by a barrier; qubits not read out idle exactly that length). Any extra elapsed time per
  readout step, or any other cause tied to the schedule, is the same for K = 1, 10 and 39.
- Per untouched qubit: loss(K) = coherence after N = 8 idle steps minus after N = 8 readout steps, each against its own arm's idle (copies pooled).
- Reading per processor: median over qubits of loss(39) - loss(1), 95 % bootstrap. ABOVE 0 -> the loss grows with the number read out at a
  fixed schedule, so extra elapsed time per step (and anything else tied to the schedule) is EXCLUDED as its explanation.
  CONSISTENT -> cannot tell. Also printed: median loss(1), loss(10), loss(39) each vs 0.
