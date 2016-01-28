# MATPLOTLIB CONFIGURATION FOR ALL PLOTS
# Execute this with the following command:
# >>> execfile('<path-to-this-file>')

import matplotlib

params={
    'axes.labelsize' : 14,
    'axes.titlesize' : 16,
    'legend.fontsize' : 14,
    'legend.fancybox' : True,
    'legend.framealpha' : 0.5,
}
matplotlib.rcParams.update(params)
