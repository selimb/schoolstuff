import matplotlib
def custom():
    params = {
        'axes.labelsize'  : 18,
        'axes.titlesize'  : 16,
        'legend.fontsize' : 15,
        'xtick.labelsize' : 12,
        'ytick.labelsize' : 12,
        'font.family'     : 'serif',
    }
    matplotlib.rcParams.update(params)
