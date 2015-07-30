"""Implement utility to solve general Riemann problem exactly.
"""
# TODO: Finish docstring above.
# QUESTION: Should input x/t implicitly take into account x0?
from __future__ import division, print_function
import numpy as np
from .utils import compute_c_from_primitives
from . import waves


class Riemann(object):

    """Represents a Riemann Problem.

    Attributes
    ----------
    W_l : [rho, u, p] array_like
        Left primitive values.
    c_l : float
        Sound speed on left.
    W_r : [rho, u, p] array_like
        Right primitive values.
    c_r : float
        Sound speed on right.
    gamma : float
        Heat capacity ratio.
    x0 : float
        Initial discontinuity position.
    p_star : float
        Pressure in middle.
    u_star : float
        Velocity in middle.
    left : gasdynamics.Wave
        Left Wave.
    right : gasdynamics.Wave
        Right Wave.
    """

    def __init__(self, W_l, W_r, gamma, x0=0):
        """Define Riemann problem using left and right values and gamma.

        Parameters
        ----------
        W_l : [rho, u, p] array_like
            Left primitive values.
        W_r : [rho, u, p] array_like
            Right primitive values.
        gamma : float
            Heat capacity ratio.
        x0 : float, optional
            Initial discontinuity position.
        """
        # TODO: Test store c_l, c_r versus compute them in p_star.
        # TODO: Raise error for negative pressure/density (not velocity!)
        self.W_l = W_l
        self.W_r = W_r
        self.gamma = gamma
        self.x0 = x0
        self.c_l, self.c_r = compute_c_from_primitives(self.W_l, self.W_r,
                                                       self.gamma)
        self.p_star = -1
        self.u_star = -1
        self.left_wave = None
        self.right_wave = None

    def update(self):
        """Update p_star, u_star and determine wave configuration."""
        self.p_star, self.u_star = waves.compute_star(
            W_l=self.W_l, W_r=self.W_r, c_l=self.c_l, c_r=self.c_r,
            gamma=self.gamma)
        # Determine wave configuration.
        self.left_wave = self._choose_wave(self.W_l, self.c_l, -1)
        self.right_wave = self._choose_wave(self.W_r, self.c_r, 1)
        self.left_wave.update_speed()
        self.right_wave.update_speed()

    def _choose_wave(self, W0, c0, sgn):
        kwargs = dict(W0=W0, c0=c0, g=self.gamma, sgn=sgn, p_star=self.p_star,
                      u_star=self.u_star, x0=self.x0)
        if W0[2] < self.p_star:
            return waves.Shockwave(**kwargs)
        else:
            return waves.Fan(**kwargs)

    def evaluate_single_state(self, x, t):
        """Evaluate primitive values at single unknown state.

        Parameters
        ----------
        x : float
            Position of unknown state.
        t : float
            Time at unknown state.

        Returns
        -------
        W : [rho,  u, p] row ndarray
            Primitive values at unknown state.
        """
        if self.p_star == -1:
            self.update()
        S = (x - self.x0)/t
        if S < self.u_star:
            wave = self.left_wave
        else:
            wave = self.right_wave
        return wave.evaluate_single_state(S)

    def evaluate_states(self, x, t):
        """Evaluate primitive values at unknown states.

        Parameters
        ----------
        x : float or row ndarray
            Position of unknown states.
        t : float
            Time at unknown states.

        Returns
        -------
        W : [rho, u, p] ndarray
            Primitive values at unknown states. Shape is Nx3, where N = len(x).
        """
        if self.p_star == -1:
            self.update()
        S = (x - self.x0)/t
        try:
            len(S)
        except TypeError:
            return self.evaluate_single_state(x, t)
        S_left = S < self.u_star
        S_right = np.logical_not(S_left)
        return np.vstack((self.left_wave.evaluate_states(S[S_left]),
                          self.right_wave.evaluate_states(S[S_right])))

if __name__ == '__main__':
    # From bin/RStest3.txt
    # x0 = 0
    gamma = 1.4
    x0 = 0
    W_l = np.array([1.0, 0.0, 0.01])
    W_r = np.array([1.0, 0.0, 100.0])
    riemann = Riemann(W_l, W_r, gamma, x0)
    x = 2
    t = 0.035
    riemann.evaluate_state(x, t)
