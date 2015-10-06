import numpy as np
from matplotlib.ticker import AutoMinorLocator

def pi_plot(ax, pi1, pi2):
    a = np.pi*pi1
    b = np.pi*pi2
    ax.set_xlim(a, b)
    ticks = np.arange(a, b + np.pi/2, np.pi/2)
    make_label = lambda p: '$%.1f\pi$' % (p/np.pi)
    labels = list(map(make_label, ticks))
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    minor_locator = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.grid(which='minor')

def plot_init(ax, x, f, x0, legend=True):
    y = f(x)
    ax.plot(x, y)
    ax.plot(x0, f(x0), 'r.', label='$x_0$')
    if legend:
        ax.legend(loc='best', fancybox=True, numpoints=1)

def showoff(ax, all_x, f, fp):
    ltxlabel = lambda n: '$x_%d$' % n
    a, b = ax.get_xlim()
    ticks = []
    tick_labels = []
    x0 = all_x[0]
    y0 = f(x0)
    ticks.append(x0)
    tick_labels.append(ltxlabel(0))
    for i in range(len(all_x) - 1):
        xi = all_x[i]
        d = fp(xi)
        c = f(xi) - d*xi
        y = lambda x: d*x + c
        ya = y(a)
        yb = y(b)
        r = -c/d
        yr = f(all_x[i + 1])
        ax.plot([a, b], [ya, yb], '--', linewidth=1)
        ax.plot((r, r), (0.0, yr), 'k-', linewidth=1)
        ax.plot(r, yr, 'r.')
        ticks.append(r)
        tick_labels.append(ltxlabel(i+1))
    ax.grid(which='both')
    ax.grid(which='major', axis='y')
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=16, verticalalignment='top')
    ax.set_yticks([0])

