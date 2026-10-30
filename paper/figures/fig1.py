"""Fig. 1. (a) One step of the idle arm and of the readout arm on a shared time axis of exactly the readout length L.
(b) ibm_fez coupling map (156 qubits, heavy-hex) with the 39 measured qubits and the 117 untouched qubits.
Inputs (read only): control_run/backend_properties_ibm_fez.json (cz couplers),
~/bones-forge/heuremen.org/experiments/2026-10-04-three-positions/data/reads/countwide_roles_ibm_fez.json (the ROLES file of control_run/run.py),
control_run/job_ibm_fez.json (E, S and L_dt as written into the job; checked against the roles file)."""
import json, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
import figstyle as fs

fs.setup()
CR = os.path.join(fs.ROOT, 'control_run')
props = json.load(open(os.path.join(CR, 'backend_properties_ibm_fez.json')))
roles = json.load(open(os.path.join(fs.FORGE, 'heuremen.org/experiments/2026-10-04-three-positions/data/reads/countwide_roles_ibm_fez.json')))
job = json.load(open(os.path.join(CR, 'job_ibm_fez.json')))
E = sorted(roles['emitters']); n = roles['n']
assert sorted(job['E']) == E and len(E) == 39 and len(job['S']) == 117 and n == 156
edges = sorted({tuple(sorted(g['qubits'])) for g in props['gates'] if g['gate'] == 'cz'})
adj = {}
for a, b in edges: adj.setdefault(a, set()).add(b); adj.setdefault(b, set()).add(a)

# heavy-hex coordinates from the coupling list: 20-qubit blocks = a row of 16 plus 4 bridge qubits below it
pos = {}
for q in range(n):
    blk, r = divmod(q, 20)
    if r < 16: pos[q] = (r, -2.0 * blk)
for q in range(n):
    if q not in pos:
        up = min(adj[q]); pos[q] = (pos[up][0], pos[up][1] - 1.0)
for a, b in edges:   # every coupler must join neighbouring sites in this layout
    (xa, ya), (xb, yb) = pos[a], pos[b]; assert abs(xa - xb) + abs(ya - yb) == 1, (a, b)

fig = plt.figure(figsize=(fs.COL1, 3.75))
axa = fig.add_axes([0.0, 0.53, 1.0, 0.47]); axb = fig.add_axes([0.0, 0.0, 1.0, 0.52])

# ---------------- (a) schematic ----------------
ax = axa; ax.set_xlim(-1.35, 9.0); ax.set_ylim(0.55, 5.75); ax.axis('off')
x_prep, x_b1, x_b2, x_rot, x_meas, x_end = 0.0, 1.0, 5.0, 6.75, 7.75, 8.75   # step spans [x_b1, x_b2] = L
L_dt = job['L_dt']

def barrier(x, y0, y1):
    ax.add_patch(Rectangle((x - 0.07, y0), 0.14, y1 - y0, facecolor='#d9d9d9', edgecolor='none', zorder=1))
    ax.plot([x, x], [y0, y1], color='k', lw=0.6, ls=(0, (2, 1.5)), zorder=2)

def gate(x, y, txt, w=0.62, h=0.46, fc='white', ec='k', fs_=6.5, color='k'):
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, facecolor=fc, edgecolor=ec, lw=0.7, zorder=3))
    ax.text(x, y, txt, ha='center', va='center', fontsize=fs_, color=color, zorder=4)

def meter(x, y, w=0.62, h=0.46, ec='k'):
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, facecolor='white', edgecolor=ec, lw=0.7, zorder=3))
    t = np.linspace(np.pi * 0.15, np.pi * 0.85, 30); r = 0.2
    ax.plot(x + r * np.cos(t), y - 0.12 + r * np.sin(t), color=ec, lw=0.6, zorder=4)
    ax.plot([x, x + 0.16], [y - 0.12, y + 0.13], color=ec, lw=0.6, zorder=4)

arms = [('idle arm', 4.25, 'D'), ('readout arm', 1.85, 'M')]
for label, y0, mode in arms:
    yE, yS = y0 + 0.5, y0 - 0.25
    ax.text(-0.75, y0 + 0.92, label, ha='left', va='bottom', fontsize=7.5, fontweight='bold', color=fs.IDLE if mode == 'D' else fs.READ)
    ax.text(-0.45, yE, 'measured\n(39)', ha='right', va='center', fontsize=6.3)
    ax.text(-0.45, yS, 'untouched\n(117)', ha='right', va='center', fontsize=6.3)
    for y in (yE, yS): ax.plot([-0.35, x_end], [y, y], color='k', lw=0.6, zorder=0)
    gate(x_prep + 0.25, yS, r'$\sqrt{X}$')
    barrier(x_b1, yS - 0.35, yE + 0.35); barrier(x_b2, yS - 0.35, yE + 0.35)
    if mode == 'D':
        ax.add_patch(Rectangle((x_b1 + 0.12, yE - 0.23), x_b2 - x_b1 - 0.24, 0.46, facecolor='#d1e3f0', edgecolor=fs.IDLE, lw=0.7, zorder=3))
        ax.text((x_b1 + x_b2) / 2, yE, 'delay $L$', ha='center', va='center', fontsize=6.5, zorder=4)
    else:
        ax.add_patch(Rectangle((x_b1 + 0.12, yE - 0.23), x_b2 - x_b1 - 0.24, 0.46, facecolor='#f7dccb', edgecolor=fs.READ, lw=0.7, zorder=3))
        ax.text((x_b1 + x_b2) / 2, yE, 'measure (duration $L$)', ha='center', va='center', fontsize=6.5, zorder=4)
    ax.text((x_b1 + x_b2) / 2, yS + 0.06, 'no gates', ha='center', va='bottom', fontsize=6, color=fs.GREY, style='italic')
    ax.text(x_b2 + 0.42, y0 + 0.12, r'$\times N$', ha='center', va='center', fontsize=7)
    gate(x_rot, yS, r'$X/Y$', w=0.72)
    meter(x_meas, yS)
# shared time axis of exactly L under the steps
ya = 1.05
ax.annotate('', xy=(x_b1, ya), xytext=(x_b2, ya), arrowprops=dict(arrowstyle='<->', lw=0.6, shrinkA=0, shrinkB=0))
ax.text((x_b1 + x_b2) / 2, ya - 0.08, r'one step $=L$ (%d dt $=1.660\,\mu$s on fez)' % L_dt, ha='center', va='top', fontsize=6.3)
for x in (x_b1, x_b2): ax.plot([x, x], [ya, 1.25], color=fs.GREY, lw=0.4, ls=':')
ax.text(x_b2 + 0.15, 5.5, 'barriers (all qubits)', ha='left', va='top', fontsize=6, color=fs.GREY)
ax.text(-1.35, 5.75, '(a)', ha='left', va='top', fontsize=9, fontweight='bold')

# ---------------- (b) coupling map ----------------
ax = axb; Eset = set(E)
for a, b in edges:
    (xa, ya_), (xb, yb) = pos[a], pos[b]; ax.plot([xa, xb], [ya_, yb], color='#bbbbbb', lw=0.7, zorder=1)
xs = np.array([pos[q][0] for q in range(n)]); ys = np.array([pos[q][1] for q in range(n)]); isE = np.array([q in Eset for q in range(n)])
ax.scatter(xs[~isE], ys[~isE], s=11, facecolor='white', edgecolor=fs.BLUE, linewidth=0.8, zorder=2, label='untouched (117)')
ax.scatter(xs[isE], ys[isE], s=14, facecolor=fs.VERM, edgecolor=fs.VERM, linewidth=0.6, marker='s', zorder=3, label='measured (39)')
ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-2.2, 15.6); ax.set_ylim(-14.6, 0.6)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.0), ncol=2, handletextpad=0.3, columnspacing=1.5, fontsize=7)
ax.text(-0.04, 0.99, '(b)', transform=ax.transAxes, ha='left', va='top', fontsize=9, fontweight='bold')
fs.save(fig, 1)
