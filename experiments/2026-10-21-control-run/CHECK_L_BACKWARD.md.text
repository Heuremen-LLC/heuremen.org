# Check L: do untouched qubits spinning backward get the kick too? With intervals, on the two timestamped runs.
Defined 22 Oct 2026 before computing; the readout-vs-idle run printed this as a median without an interval. Heuremen LLC. Script check_l_backward.py.
Data: arm S of the control run (standard bits) and of the raw-readout run (IQ projection), same 117 spectators, idle D and readout M at
N = 2, 4, 8, both copies pooled. Per spectator: spin rate s_D = slope of the idle phase over N (deg per step; phase = atan2(<Y>, <X>), unwrapped);
kick = slope of (readout phase - idle phase) over N (deg per step). Groups as declared in the readout-vs-idle plan: forward s_D >= +1,
backward s_D <= -1 deg per step. Median kick per group, 95 % bootstrap (2000), each against 0, and backward minus forward.
Extra elapsed time per step predicts a kick proportional to s_D: opposite signs in the two groups. A real kick that ignores the spin predicts the
same sign in both and backward minus forward CONSISTENT. Sign carried: the forward group's median sign is printed and the backward group is read
against it.
