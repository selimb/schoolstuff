# MATPLOTLIB CONFIGURATION FOR ALL PLOTS
# Execute this with the following command:
# >>> execfile('<path-to-this-file>')

import matplotlib

params={
    'backend': 'ps',
    'axes.labelsize' : 15,
    'axes.titlesize' : 18,
    'legend.fontsize' : 14,
    'legend.fancybox' : True,
    'legend.framealpha' : 0.5,
    'legend.loc' : 'best',
    'text.latex.preamble': [r'\usepackage{gensymb}'],
    'text.usetex': True,
    'font.family': 'serif'
}
matplotlib.rcParams.update(params)
