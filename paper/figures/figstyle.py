"""Shared style for the figures of 'One Readout, Three Clocks' (revtex4-2, two-column).
Okabe-Ito colourblind-safe palette; 8-9 pt fonts; vector PDF + PNG output."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# Public layout: this folder is paper/figures in the heuremen.org repository. The scripts were written against short run names (control_run, deciding_runs, ...).
# A folder of links maps each short name to the dated folder under experiments/ (2026-10-21-control-run -> control_run). Files are looked up in this repository first
# and then in READOUT_EXPERIMENTS, the experiments/ folder of the unpacked raw bundle (Zenodo concept DOI 10.5281/zenodo.23169755), which holds the large .npz arrays.
SITE = os.path.dirname(os.path.dirname(HERE))
def _farm():
    farm = os.path.join(HERE, '_runs'); os.makedirs(farm, exist_ok=True)
    roots = [os.environ.get('READOUT_EXPERIMENTS'), os.path.join(SITE, 'experiments')]
    for root in [r for r in roots if r and os.path.isdir(r)]:
        for d in sorted(os.listdir(root)):
            m = re.match(r'\d{4}-\d{2}-\d{2}-(.+)', d)
            if not m or not os.path.isdir(os.path.join(root, d)): continue
            short = m.group(1).replace('-', '_'); dst = os.path.join(farm, short); os.makedirs(dst, exist_ok=True)
            for f in os.listdir(os.path.join(root, d)):
                for name in ([f, f[:-5] + '.txt'] if f.endswith('.text') else [f]):      # the site stores text results with a .text suffix, the scripts ask for .txt
                    link = os.path.join(dst, name)
                    if not os.path.lexists(link): os.symlink(os.path.join(root, d, f), link)
    # in the original working folder some results were kept together; here they sit in their own dated folders
    for short, extras in (('control_run', ['iq_run', 'pairs_run']), ('excitation_kingston', ['kingston_claim4_repeat'])):
        dst = os.path.join(farm, short); os.makedirs(dst, exist_ok=True)
        for e in extras:
            src = os.path.join(farm, e)
            if not os.path.isdir(src): continue
            for f in os.listdir(src):
                link = os.path.join(dst, f)
                if not os.path.lexists(link): os.symlink(os.path.realpath(os.path.join(src, f)), link)
    return farm
import re
ROOT = _farm(); FORGE = SITE
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
