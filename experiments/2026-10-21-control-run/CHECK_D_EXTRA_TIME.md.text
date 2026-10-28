# Check D: could readout steps that secretly last longer than the idle steps account for the superposition cost?
Defined after the data, before computing. Heuremen LLC, 21 Oct 2026. Script check_d_extra_time.py.
- Extra time dt per readout step would add phase rate*dt to every untouched qubit, so the kick would rise with precession rate with
  slope a = dt / L. The upper end of the 95 % interval of a from the committed analysis (R3, both cross-copy fits; the larger is used)
  gives dt_max = max(a_hi, 0) * L.
- Extra time can cost coherence at most c_D(N) * (1 - exp(-N * dt_max / T2)) for each untouched qubit, with T2 from IBM's calibration saved
  at firing (backend_properties_<chip>.json). Upper bound uses each qubit's own T2 and the idle-arm coherence c_D(N).
- Reading per processor and N (2, 4, 8): median over untouched qubits of (measured loss - extra-time bound).
  ABOVE 0 -> extra elapsed time cannot account for the loss (EXCLUDED); otherwise NOT EXCLUDED.
