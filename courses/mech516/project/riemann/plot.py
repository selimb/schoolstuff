"""Plotting utilities."""

import matplotlib.pyplot as plt
import numpy as np
from . import waves

STYLES = dict(shock={'linestyle': 'solid',
                     'color': 'red'},
              cs={'linestyle': 'dashed',
                  'color': 'black',
                  'linewidth': 2,
                  'dashes': (15, 10)},
              fan={'linestyle': 'solid',
                   'color': 'blue'})


def plot_line(x0, S, t, ax=None, **kwargs):
    """Plot generic straight line.

    Parameters
    ----------
    x0 : float
        Origin of wave.
    S : float or column ndarray
        Wave speed.
    t : ndarray
        Array of time onto which wave position is mapped.
    ax : matplotlib.pyplot.Axes, optional
        Axes to plot on.
    **kwargs
        Line2D properties.
    """
    # Transpose needed if S is ndarray. Does not affect S = float.
    x = (t * S + x0).T
    args = (x, t)
    if ax:
        ax.plot(*args, **kwargs)
    else:
        plt.plot(*args, **kwargs)


def plot_shock(x0, S,  t, ax=None):
    """Plot shockwave.

    Parameters
    ----------
    x0 : float
        Origin of wave.
    S : float
        Shock speed.
    t : ndarray
        Array of time onto which wave position is mapped.
    ax : matplotlib.pyplot.Axes, optional
        Axes to plot on.
    """
    style = 'shock'
    plot_line(x0, S, t, ax, **STYLES[style])


def plot_fan(x0, S_head, S_tail, t, ax=None, N=5):
    """Plot expansion fan.

    Parameters
    ----------
    x0 : float
        Origin of wave.
    S_head : float
        Head expansion wave speed.
    S_tail : float
        Tail expansion wave speed.
    t : ndarray
        Array of time onto which wave position is mapped.
    ax : matplotlib.pyplot.Axes, optional
        Axes to plot on.
    N : int, optional
        Number of waves to plot, including head and tail.
    """
    style = 'fan'
    S = np.linspace(S_head, S_tail, N).reshape(N, 1)
    plot_line(x0, S, t, ax, **STYLES[style])


def plot_cs(x0, S, t, ax=None):
    """Plot contact surface.

    Parameters
    ----------
    x0 : float
        Origin of contact surface.
    S : float
        Contact surface speed.
    t : ndarray
        Array of time onto which wave position is mapped.
    ax : matplotlib.pyplot.Axes, optional
        Axes to plot on.
    """
    style = 'cs'
    plot_line(x0, S, t, ax, **STYLES[style])


def plot_riemann(riemann, t, ax=None):
    """Plot waves of particular Riemann problem.

    Parameters
    ----------
    riemann : Riemann
        Riemann problem to plot.
    t : ndarray
        Array of time onto which wave position is mapped.
    ax : matplotlib.pyplot.Axes, optional
        Axes to plot on.
    """
    riemann.update()
    plot_cs(riemann.x0, riemann.u_star, t, ax)
    for wave in [riemann.left_wave, riemann.right_wave]:
        if isinstance(wave, waves.Shockwave):
            plot_shock(riemann.x0, wave.speed, t, ax)
        else:
            plot_fan(riemann.x0, wave.speed[0], wave.speed[1], t, ax)

def plot_particles(riemann, x0, t_max, N=50, ax=None):
    """Plot particle trajectories of particular Riemann problem.

    Parameters
    ----------
    riemann : Riemann
        Riemann problem.
    x0 : ndarray
        Initial particle positions.
    t_max : float
        Time to plot until.
    N : int, optional
        Number of samples to generate. Default is 50.
    ax : matplotlib.pyplot.Axes, optional
        Axes to plot on.
    """
    riemann.update()
    time = np.linspace(0, t_max, num=N)
    t_step = time[1] - time[0]
    x = np.empty((N, len(x0)))
    x[0, :] = x0
    for i, t in enumerate(time[:-1]):
        u = riemann.evaluate_states(x[i], t)[:, 1]
        x[i+1] = x[i] + u*t_step
    if ax:
        ax.plot(x, time)
    else:
        plt.plot(x, time)


