import numpy as np
from numpy.testing import assert_allclose, assert_equal

from ..riemann import waves
from ..riemann.utils import compute_c, compute_c_from_primitives


def compare(actual, target, rtol):
    msg = "Relative Error: {0:f}".format(abs(actual - target) / abs(target))
    assert_allclose(actual, target, rtol=rtol, err_msg=msg)


def test_p_star0():
    gamma = 1.4
    W_l = np.array([1, 6, 4])
    rho_l, u_l, p_l = W_l
    W_r = np.array([2, 6, 2])
    rho_r, u_r, p_r = W_r
    c_l, c_r = compute_c([rho_l, rho_r], [p_l, p_r], gamma)
    assert_equal(rho_l * c_l, rho_r * c_r)
    pstar0 = waves._compute_p_star0(W_l, c_l, W_r, c_r)
    assert_equal(pstar0, 1 / 2 * (p_l + p_r))


def test_p_star_uniform():
    places = 6
    tol = 10**(-places)
    gamma = 1.4
    W_l = np.array([1, 2, 3])
    W_r = np.array([1, 2, 3])
    c_l, c_r = compute_c_from_primitives(W_l, W_r, gamma)
    pstar, ustar = waves.compute_star(W_l, c_l, W_r, c_r, gamma, tol=tol)
    compare(W_l[2], pstar, rtol=10**-4)
    compare(W_l[1], ustar, rtol=10**-4)


def test_p_star_1():
    # From bin/RStest1.txt
    gamma = 1.4
    W_l = np.array([5.99924, 19.5975, 460.894])
    W_r = np.array([5.99242, -6.19633, 46.0950])
    target_p_star, target_u_star = 0.16916 * 10**4, 0.86899 * 10**1
    c_l, c_r = compute_c_from_primitives(W_l, W_r, gamma)
    p_star, u_star = waves.compute_star(W_l, c_l, W_r, c_r, gamma, tol=10**-7)
    compare(p_star, target_p_star, rtol=0.002)
    compare(u_star, target_u_star, rtol=0.002)


def test_p_star_2():
    # From bin/RStest2.txt
    gamma = 1.4
    W_l = np.array([1.0, 0.75, 1.0])
    W_r = np.array([0.125, 0.0, 0.1])
    target_p_star, target_u_star = 0.46667 * 10**0, 0.13618 * 10**1
    c_l, c_r = compute_c_from_primitives(W_l, W_r, gamma)
    p_star, u_star = waves.compute_star(W_l, c_l, W_r, c_r, gamma, tol=10**-7)
    compare(p_star, target_p_star, rtol=0.002)
    compare(u_star, target_u_star, rtol=0.002)


def test_p_star_3():
    # From bin/RStest3.txt
    gamma = 1.4
    W_l = np.array([1.0, 0.0, 0.01])
    W_r = np.array([1.0, 0.0, 100.0])
    target_p_star, target_u_star = 0.46160 * 10**2, -0.62007 * 10**1
    c_l, c_r = compute_c_from_primitives(W_l, W_r, gamma)
    p_star, u_star = waves.compute_star(W_l, c_l, W_r, c_r, gamma, tol=10**-7)
    compare(p_star, target_p_star, rtol=0.002)
    compare(u_star, target_u_star, rtol=0.002)


def test_star_vacuum():
    gamma = 1.4
    W_l = np.array([1.0, -20, 1.0])
    W_r = np.array([1.0, 20, 1.0])
    c_l, c_r = compute_c_from_primitives(W_l, W_r, gamma)
    try:
        waves.compute_star(W_l, c_l, W_r, c_r, gamma)
    except waves.VacuumCondition:
        return
    raise AssertionError


def test_star_almost_vacuum():
    # This is actually from Problem 2 of Mini-Project 1.
    W_l = np.array([1.0, -2.0, 0.4])
    W_r = np.array([1.0, 2.0, 0.4])
    gamma = 1.4
    c_l, c_r = compute_c_from_primitives(W_l, W_r, gamma)
    p_star, u_star = waves.compute_star(W_l, c_l, W_r, c_r, gamma)
    assert(p_star > 0.0)
    assert_equal(u_star, 0.0)
