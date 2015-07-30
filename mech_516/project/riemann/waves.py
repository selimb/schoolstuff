"""Equations related to waves (rarefaction and shock)"""

from __future__ import division, print_function
import abc
from math import sqrt
import numpy as np


#########################
#  Evaluation of p_star #
#########################
class VacuumCondition(Exception):

    def __init__(self):
        self.msg = "Vacuum State not supported."


def _compute_p_star0(W_l, c_l, W_r, c_r):
    """Compute first guess for p_star.

    Use acoustic wave approximation.

    Parameters
    ----------
    W_l : [rho, u, p] array_like
        Left primitive values.
    c_l : float
        Sound speed on left.
    W_r : [rho, u, p] array_like
        Right primitive values.
    c_r : float
        Sound speed on right.

    Returns
    -------
    p_star0 : float
        First approximation for p.
    """
    rho_l, u_l, p_l = W_l
    rho_r, u_r, p_r = W_r
    num = p_l*rho_r*c_r + p_r*rho_l*c_l + (u_l - u_r)*rho_l*c_l*rho_r*c_r
    den = rho_l * c_l + rho_r * c_r
    return num / den


def compute_star(W_l, c_l, W_r, c_r, gamma, tol=10**-6):
    """Compute p_star, u_star iteratively using newton-rapshon method.

    Parameters
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
    tol : float, optional
        Minimum error/tolerance between two iterations.

    Returns
    -------
    p_star : float
    u_star : float
    """
    rho_l, u_l, p_l = W_l
    rho_r, u_r, p_r = W_r

    def shock(p_star, rho, u, p, c):
        """Compute f_k(p_star, W) and f'_k(p_star, W) for shockwave."""
        # f  = (b_num/b_den)**0.5 * a = b_5*a
        # f' = (b_num/b_den)**0.5 * (1 - a/(2*b_den)) = b_5*(1 - a/(2*b_den))
        a = p_star - p
        b_num = 2/((gamma + 1) * rho)
        b_den = p_star + p*(gamma - 1)/(gamma + 1)
        b_5 = (b_num/b_den)**0.5
        return b_5*a, b_5*(1 - a/(2*b_den))

    def fan(p_star, rho, u, p, c):
        """Compute f_k(p_star, W) and f'_k(p_star, W) for fan wave."""
        # f  = a*(b**pwr - 1) = a*(b_pwr - 1)
        # f' = c*b_pwr/b
        a = 2*c/(gamma - 1)
        b = (p_star/p)
        pwr = (gamma - 1)/(2*gamma)
        b_pwr = b**pwr
        c = 1/(rho*c)
        return a*(b_pwr - 1), c*b_pwr/b

    def wave_cases(p_star, rho, u, p, c):
        """Compute f_k(p_star, W) and f'_k(p_star, W) depending on p_star
        and p."""
        if p_star > p:
            f_k, fder_k = shock(p_star, rho, u, p, c)
        else:
            f_k, fder_k = fan(p_star, rho, u, p, c)
        return f_k, fder_k

    def compute_f_and_fder(p_star):
        """Compute f(p_star, W_l, W_r) and f'(p_star, W_L, W_r)."""
        f_l, fder_l = wave_cases(p_star, rho_l, u_l, p_l, c_l)
        f_r, fder_r = wave_cases(p_star, rho_r, u_r, p_r, c_r)
        return f_l + f_r + u_r - u_l, fder_l + fder_r

    def newton(p0):
        """Newton root-finding method.

        Parameters
        ----------
        p0 : float
            First approximation for p_star.
        """
        for i in range(100):  # 100 maximum iterations.
            fval, fder = compute_f_and_fder(p0)
            p = p0 - fval/fder
            # if abs(p - p0)/p0 < tol:
            if abs(p - p0) < tol:
                return p
            p0 = p

    # Compute p_star
    p0 = _compute_p_star0(W_l, c_l, W_r, c_r)
    if p0 < 0.0:
        # Check vacuum condition
        if (W_r[1] - W_l[1] >= 2/(gamma - 1)*(c_l + c_r)):
            raise VacuumCondition
        else:
            p0 = 0.001
    p_star = newton(p0)

    # Compute u_star from coupled equation
    f_l, fder_l = wave_cases(p_star, rho_l, u_l, p_l, c_l)
    f_r, fder_r = wave_cases(p_star, rho_r, u_r, p_r, c_r)
    u_star = 0.5*(u_l + u_r + f_r - f_l)
    return p_star, u_star

#################
# Riemann waves #
#################


class BaseWave(object):
    # TODO: Reword docstring

    """Base wave representation for Riemann problem where two states are
       separated by a wave.

    Attributes
    ----------
    W0 : [rho, u, p] array_like
        Primitive values ahead of the wave.
    c0 : [rho, u, p] array_like
        Speed of sound ahead of the wave.
    g : float
        Heat capacity ratio.
    p_star : float
        Pressure behind wave.
    u_star : float
        Velocity behind wave.
    sgn : int
        Sign for distinguishing "left" (-1) from "right" (1) wave.
    x0 : float
        Initial wave position.
    speed : float
        Wave speed. Only computed when update_speed() method is called.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, W0,  c0, g, p_star, u_star, sgn, x0=0):
        """Constructor for base wave class.

        Parameters
        ----------
        W0 : [rho, u, p] array_like
            Primitive values ahead of the wave.
        c0 : [rho, u, p] array_like
            Speed of sound ahead of the wave.
        g : float
            Heat capacity ratio.
        p_star : float
            Pressure behind wave.
        u_star : float
            Velocity behind wave.
        sgn : int
            Sign for distinguishing "left" (-1) from "right" (1) wave.
        x0 : float, optional
            Initial wave position.
        """
        self.W0 = W0
        self.c0 = c0
        self.g = g
        self.p_star = p_star
        self.u_star = u_star
        self.sgn = sgn
        self.x0 = x0
        self.speed = None

    @abc.abstractmethod
    def update_speed(self):
        """Compute wave speed. """
        pass

    @abc.abstractmethod
    def evaluate_single_state(self, S):
        """Evaluate primitive values at single unknown state.

        Parameters
        ----------
        S : float or ndarray
            Speed of unknown state.

        Returns
        -------
        W : [rho, u, p] ndarray
            Primitive values at unknown state.
        """
        if not self.speed:
            self.update_speed()

    @abc.abstractmethod
    def evaluate_states(self, S):
        """Evaluate primitive values at unknown states.

        Parameters
        ----------
        S : float or ndarray
            Speed of unknown states.

        Returns
        -------
        W : [rho, u, p] ndarray
            Primitive values at unknown states. Shape is Nx3, where N = len(S).
        """
        if not self.speed:
            self.update_speed()

    @abc.abstractmethod
    def __repr__(self):
        pass


class Shockwave(BaseWave):

    """Shockwave."""

    def update_speed(self):
        Ms_a = (self.g + 1)*self.p_star/(2*self.g*self.W0[2])
        Ms_b = (self.g - 1)/(2*self.g)
        Ms = sqrt(Ms_a + Ms_b)
        self.speed = self.W0[1] + self.sgn*self.c0*Ms

    def evaluate_single_state(self, S):
        super(Shockwave, self).evaluate_single_state(S=S)
        # Outside
        if self.sgn*S > self.sgn*self.speed:
            return self.W0

        # Inside
        return self._compute_inside()

    def evaluate_states(self, S):
        super(Shockwave, self).evaluate_states(S=S)
        try:
            rows = len(S)
        except TypeError:  # Float
            rows = 1
        cols = 3
        W = np.empty((rows, cols))
        S_outside = self.sgn*S > self.sgn*self.speed
        S_inside = self.sgn*S < self.sgn*self.speed
        # Outside
        W_outside = self.W0
        # Inside
        W_inside = self._compute_inside()
        # Assign
        try:
            W[S_outside] = W_outside
            W[S_inside] = W_inside
        except IndexError:  # Float
            if S_outside:
                W = W_outside
            else:
                W = W_inside
        return W

    def _compute_inside(self):
        g_m1_p1 = (self.g - 1)/(self.g + 1)
        p_ratio = self.p_star/self.W0[2]
        rho_star_num = p_ratio + g_m1_p1
        rho_star_den = g_m1_p1*p_ratio + 1
        rho_star = self.W0[0]*rho_star_num/rho_star_den
        return np.array([rho_star, self.u_star, self.p_star])

    def __repr__(self):
        side = "Left" if self.sgn == -1 else "Right"
        return ("{0:s} Shock Wave.".format(side))


class Fan(BaseWave):

    """Expansion fan.

    Attribute 'speed' is a tuple (S_head, S_tail).
    """

    def update_speed(self):
        S_head = self.W0[1] + self.sgn*self.c0
        c_star = self.c0*(self.p_star/self.W0[2])**((self.g - 1)/(2*self.g))
        S_tail = self.u_star + self.sgn*c_star
        self.speed = (S_head, S_tail)

    def evaluate_single_state(self, S):
        super(Fan, self).evaluate_single_state(S)
        # Outside
        if self.sgn*S > self.sgn*self.speed[0]:
            return self.W0

        # Inside
        if self.sgn*S < self.sgn*self.speed[1]:
            return self._compute_inside()

        # In fan
        return self._compute_in_fan(S)

    def evaluate_states(self, S):
        super(Fan, self).evaluate_states(S=S)
        try:
            rows = len(S)
        except TypeError:  # Float
            rows = 1
        cols = 3
        W = np.empty((rows, cols))
        S_outside = self.sgn*S > self.sgn*self.speed[0]
        S_inside = self.sgn*S < self.sgn*self.speed[1]
        S_in_fan = np.logical_not(np.logical_or(S_outside, S_inside))
        # Outside
        W_outside = self.W0
        # Inside
        W_inside = self._compute_inside()
        # In fan
        try:
            W_in_fan = self._compute_in_fan(S[S_in_fan]).T
        except TypeError:  # Float
            W_in_fan = self._compute_in_fan(S)
        # Assign
        try:
            W[S_outside] = W_outside
            W[S_in_fan] = W_in_fan
            W[S_inside] = W_inside
        except IndexError:  # Float
            if S_outside:
                W = W_outside
            elif S_in_fan:
                W = W_in_fan
            else:
                W = W_inside
        return W

    def _compute_inside(self):
        rho_star = self.W0[0]*(self.p_star/self.W0[2])**(1/self.g)
        return np.array([rho_star, self.u_star, self.p_star])

    def _compute_in_fan(self, S):
        gp1 = self.g + 1
        gm1 = self.g - 1
        rho_temp = 2/gp1 - self.sgn*gm1/(gp1*self.c0)*(self.W0[1] - S)
        rho = (self.W0[0]*rho_temp**(2/gm1))
        u = 2/gp1*(gm1/2*self.W0[1] + S - self.sgn*self.c0)
        p_temp = 2/gp1 - self.sgn*gm1/(gp1*self.c0)*(self.W0[1] - S)
        p = self.W0[2]*p_temp**(2*self.g/gm1)
        return np.array([rho, u, p])

    def __repr__(self):
        side = "Left" if self.sgn == -1 else "Right"
        return ("{0:s} Rarefaction Wave.".format(side))
