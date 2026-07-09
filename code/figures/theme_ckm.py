"""theme_ckm.py — matplotlib mirror of the locked CKM ggplot theme
(paper 3.eas/analysis/tools/theme_publication_ckm.R). R is not installed on this box, so Paper B
figures are rendered in matplotlib with the SAME tokens as Paper A's R figures for cross-paper
consistency: harmonized 2-tier +1pt fonts (body 9pt, panel labels 11pt bold), Arial, gridlines OFF,
white background (theme_classic equivalent), Okabe-Ito layer colors. Dims 89mm single / 183mm double
column, max 240mm height; save 300dpi PNG + vector PDF.
"""
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager

# --- fonts: Arial if available, else closest sans fallback (harmonized 2-tier +1pt) ---
_HAVE_ARIAL = any('Arial' == f.name for f in font_manager.fontManager.ttflist)
BODY, PANEL = 9, 11
mpl.rcParams.update({
    'font.family': 'Arial' if _HAVE_ARIAL else 'DejaVu Sans',
    'font.size': BODY, 'axes.titlesize': BODY, 'axes.labelsize': BODY,
    'xtick.labelsize': BODY, 'ytick.labelsize': BODY, 'legend.fontsize': BODY,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42,  # editable text in vector output
})
MM = 1/25.4
def figsize(w_mm, h_mm): return (w_mm*MM, h_mm*MM)

# --- LOCKED Okabe-Ito color tokens (layer-keyed; identical to Paper A) ---
LAYER = {'proteome':'#0072B2','transcriptome':'#009E73','metabolome':'#CC79A7','methylome':'#D55E00',
         'pooled':'#000000','neutral':'#333333'}
SEM = {'up':'#D55E00','down':'#56B4E9','ns':'#999999','yes':'#D55E00','no':'#56B4E9',
       'sig_band':'#fde0dd','equiv_band':'#ECECEC'}

def panel_label(ax, letter, dx=-0.18, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=PANEL, fontweight='bold',
            va='bottom', ha='left')

def save_ckm(fig, path_noext, dpi=300):
    fig.savefig(path_noext+'.png', dpi=dpi, bbox_inches='tight')
    fig.savefig(path_noext+'.pdf', bbox_inches='tight')
    return path_noext+'.png', path_noext+'.pdf'
