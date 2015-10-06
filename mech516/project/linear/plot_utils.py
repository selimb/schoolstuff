import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
from math import sqrt
from functools import partial


SPINE_COLOR = 'gray'

def latexify(fig_width=None, fig_height=None, columns=1):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float,  optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'text.latex.preamble' : '\\usepackage{gensymb}',
              'backend': 'ps',
              'axes.labelsize': 12, # fontsize for x and y labels (was 10)
              'axes.titlesize': 8,
              'font.size': 10, # was 10
              'legend.fontsize': 12, # was 10
              'xtick.labelsize': 11,
              'ytick.labelsize': 11,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'font.family': 'serif'
    }

    matplotlib.rcParams.update(params)

def format_axes(ax):

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    for spine in ['left', 'bottom']:
        ax.spines[spine].set_color(SPINE_COLOR)
        ax.spines[spine].set_linewidth(0.5)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_tick_params(direction='out', color=SPINE_COLOR)

    return ax

def make_five(solver, file_name):
    triangle = lambda x: 0.5 + 0.075*x
    fig_width = 20
    ratio = 1.5
    fig, axes = plt.subplots(ncols=5, sharey=True,
                             figsize=(20, fig_width/ratio))
    plot = partial(solver.plot, left=-2, right=35)
    for i, ax in enumerate(axes):
        plot(time=i, ax=ax)
        plot_exact(i, triangle, ax)
        format_axes(ax)
    plt.tight_layout()
    # plt.show()
    fig.savefig(file_name, dpi=350)

def calc_exact(time, signal, dx=0.1, left=None, right=None, x=None):
 # Exact solution
# x0 = [0..20]. At a = 2m/s, 200m have been traveled.
# xf = [200 .. 220] with same triangular shape.
  signal = np.vectorize(signal)
  start = 2*time
  end = start + 20
  if x is None:
    if left is None:
      left = start - 20
    if right is None:
      right = end + 20
    x = np.arange(left, right + dx, dx)
  x_signal_range = np.logical_and(x >= start, x <= end)
  x_signal_offset = x[x_signal_range][0]
  u_exact_signal = np.ones_like(x)*0.5
  u_exact_signal[x_signal_range] = signal(x[x_signal_range] - x_signal_offset)
  return x,  u_exact_signal

def plot_exact(time, signal, ax, dx=0.1, left=None, right=None):
  x, u_exact_signal = calc_exact(time, signal, dx=dx, left=left,
                                 right=right)
  ax.plot(x, u_exact_signal, 'r-')

