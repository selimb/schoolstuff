"""Common utilities."""

import numpy as np
import scipy as sp


def compute_c(rho, p, gamma):
    """Compute sound speed.

    Parameters
    ----------
    rho : float or ndarray
        Density.
    p : float or ndarray
        Pressure.
    gamma : float or ndarray
        Heat capacity ratio.

    Returns
    -------
    float or ndarray
        Speed of sound.
    """
    rho, p, gamma = np.array(rho), np.array(p),  np.array(gamma)
    return np.sqrt(gamma*p/rho)

def compute_c_from_primitives(W_l, W_r, gamma):
    """Compute sound speed.

    Parameters
    ----------
    W_l : [rho, u, p] array_like
        Left primitive values.
    W_r : [rho, u, p] array_like
        Right primitive values.
    gamma : float or ndarray
        Heat capacity ratio.

    Returns
    -------
    c : 1x2 ndarray
        Speed of sound.
    """
    return compute_c(np.array([W_l[0], W_r[0]]), np.array([W_l[2], W_r[2]]),
                     gamma)

