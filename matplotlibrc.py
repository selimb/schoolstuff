# MATPLOTLIB CONFIGURATION FOR ALL PLOTS
# Execute this with the following command:
# >>> execfile('<path-to-this-file>')

import matplotlib

params={
    'axes.labelsize' : 18,
    'axes.titlesize' : 18,
    'legend.fontsize' : 14,
    'legend.fancybox' : True,
    'legend.framealpha' : 0.5,
    'legend.loc' : 'best',
}
matplotlib.rcParams.update(params)
