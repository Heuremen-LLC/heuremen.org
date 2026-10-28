# Check E: is the superposition loss larger on untouched qubits wired to a measured qubit than on those with none, in the same circuits?
Defined after the data, before computing. Heuremen LLC, 21 Oct 2026. Script check_e_wired_vs_unwired.py.
- Extra elapsed time in the readout steps, or any chip-wide cause, acts on every untouched qubit; a cause acting through measured
  neighbours (a readout effect, or stuck or leaked neighbours) acts on wired qubits.
- Per untouched qubit (arm S, N = 2, 4, 8, copies pooled): implied extra time per readout tau = -T2/N * ln(c_M / c_D)
  (the extra idle time that would turn the idle coherence into the readout coherence), T2 from IBM's calibration saved at firing.
  Qubits with c_M or c_D <= 0.05 or no T2 are dropped (count printed).
- Reading: median tau over wired qubits minus median tau over unwired qubits, 95 % bootstrap resampling each group.
  ABOVE 0 -> wired qubits lose more than extra elapsed time (or any chip-wide cause) can explain for both groups: EXCLUDED as the
  explanation of the wired excess. CONSISTENT -> cannot tell. BELOW -> reported as such.
