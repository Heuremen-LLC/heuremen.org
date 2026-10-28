# Check H: does the kick fall off with distance from the measured qubits, is it the same per readout at 2 and at 8 readouts, and does the superposition cost grow with the number of readouts?
Defined from data already taken (control run, raw-readout run), before this was computed. Heuremen LLC, 22 Oct 2026. Script check_h_distance_and_count.py.
Data: arm S (39 measured, 117 untouched spectators on the equator, idle D vs readout M at N = 2, 4, 8; settings XX, YY; both copies pooled).
Run 1 = control run (standard readout bits). Run 2 = raw-readout run (bits from the IQ projection t > 0.5, as PREREG_IQ.md).
Phase = atan2(<Y>, <X>); kick per readout = wrap(phase_M - phase_D) / N in degrees. Coherence = |<X> + i<Y>|. Distance = shortest path in IBM's
cz graph (calibration saved at firing) from the spectator to the nearest measured qubit. 95 % bootstrap, 2000 resamples; four outcomes
(ABOVE / BELOW / CONSISTENT / CANNOT TELL, the last when a group has fewer than 5 qubits).

H1 distance, N = 8: median kick for spectators at distance 1, 2 and 3 or more, each against 0; and median(distance 1) - median(distance >= 2),
   bootstrap resampling within each group. ABOVE -> the kick falls off away from the measured qubits. CONSISTENT -> no fall-off resolved; the
   interval is the bound. BELOW -> larger far away.
H2 per-readout constancy: per spectator kick(N = 8) - kick(N = 2); median. CONSISTENT -> the same kick per readout at 2 and 8 readouts (the
   interval bounds any change). Also per spectator least-squares line of total phase excess against N (2, 4, 8): median intercept against 0.
   Intercept CONSISTENT -> no fixed offset per circuit beside the per-readout kick.
H3 superposition cost against readout count: per spectator loss(N) = coherence(idle) - coherence(readout); median of loss(8) - loss(2).
   ABOVE -> the cost grows with the number of readouts. Run 1 all spectators, all shots. Run 2 wired spectators, clean shots (every
   measured neighbour verified 0 on every readout, as R1 of PREREG_IQ.md).
A reading is reported for a processor when run 1 and run 2 give the same outcome; otherwise both are printed and it is CANNOT TELL.
