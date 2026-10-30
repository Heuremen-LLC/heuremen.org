"""Shared style for the figures of 'One Readout, Three Clocks' (revtex4-2, two-column).
Okabe-Ito colourblind-safe palette; 8-9 pt fonts; vector PDF + PNG output."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                      # .../arxiv_one_readout
FORGE = os.path.dirname(ROOT)                     # .../bones-forge
COL1, COL2 = 3.4, 7.0                             # single / double column width, inches

# Okabe-Ito
BLUE, VERM, GREEN, ORANGE, SKY, PURPLE, YELLOW, BLACK = ('#0072B2', '#D55E00', '#009E73', '#E69F00',
                                                          '#56B4E9', '#CC79A7', '#F0E442', '#000000')
IDLE, READ = BLUE, VERM                           # the two arms, everywhere
GREY = '#888888'
CHIPS = ('ibm_fez', 'ibm_kingston', 'ibm_marrakesh')
SHORT = {c: c.replace('ibm_', '') for c in CHIPS}

def setup():
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
        'mathtext.fontset': 'dejavuserif',
        'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8, 'legend.fontsize': 7,
        'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
        'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'in', 'ytick.direction': 'in',
        'lines.linewidth': 1.0, 'lines.markersize': 3.5, 'errorbar.capsize': 1.5,
        'legend.frameon': False, 'legend.handlelength': 1.6, 'legend.borderaxespad': 0.3,
        'savefig.dpi': 600, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    })

def save(fig, n):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(HERE, 'fig%s.%s' % (n, ext)), bbox_inches='tight', pad_inches=0.02)
    print('wrote', os.path.join(HERE, 'fig%s.pdf' % n), 'and .png')

def panel(ax, s, x=-0.02, y=1.0):
    ax.text(x, y, '(%s)' % s, transform=ax.transAxes, ha='right', va='bottom', fontsize=9, fontweight='bold')
