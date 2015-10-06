import numpy as np
from numpy.testing import assert_equal, assert_allclose

try:
    from ..riemann.riemann import Riemann
    from ..riemann.utils import compute_c_from_primitives
    from ..riemann.waves import compute_star
    from ..riemann import waves
except ValueError:
    from riemann.riemann import Riemann
    from riemann.utils import compute_c_from_primitives
    from riemann.waves import compute_star
    from riemann import waves


def compare(actual, target, rtol):
    try:
        rerr = abs(actual - target) / abs(target + 1 * 10**-16)
    except RuntimeWarning:
        msg = ''
    else:
        msg = "Relative Error: {0}".format(rerr)
    assert_allclose(actual, target, rtol=rtol, err_msg=msg)


def test_constructor():
    W_L = np.array([1, 1, 1])
    W_R = np.array([1, 1, 1])
    riemann = Riemann(W_L, W_R, 1.4)
    assert riemann


def test_p_star():
    gamma = 1.4
    W_l = np.array([1, 6, 4])
    W_r = np.array([2, 6, 2])
    c_l, c_r = compute_c_from_primitives(W_l, W_r, gamma)
    pstar_waves, ustar_waves = compute_star(W_l, c_l, W_r, c_r, gamma)
    riemann = Riemann(W_l, W_r, gamma)
    riemann.update()
    assert_equal(pstar_waves, riemann.p_star)
    assert_equal(ustar_waves, riemann.u_star)


def _test_sample_singleeval(test_file, W_l, W_r, gamma, t, x0=0):
    target_results = np.genfromtxt(test_file, skip_header=1)
    riemann = Riemann(W_l, W_r, gamma,  x0)
    for x, rho, u, p in target_results:
        W_actual = riemann.evaluate_single_state(x, t)
        W_target = np.array([rho, u, p])
        try:
            compare(W_actual, W_target, 0.002)
        except AssertionError:
            print("x:", x)
            print("t:", t)
            print("S:", x / t)
            raise


def _test_sample_multipleeval(test_file, W_l, W_r, gamma, t, x0=0):
    target_results = np.genfromtxt(test_file, skip_header=1)
    riemann = Riemann(W_l, W_r, gamma,  x0)
    x = target_results[:, 0]
    W_target = target_results[:, 1:]
    W = riemann.evaluate_states(x, t)
    assert_allclose(W, W_target, rtol=0.002)


def test_sample3_singleeval():
    test_file = 'bin/RStest3.txt'
    gamma = 1.4
    W_l = np.array([1.0, 0.0, 0.01])
    W_r = np.array([1.0, 0.0, 100.0])
    t = 0.035
    _test_sample_singleeval(test_file, W_l, W_r, gamma, t)


def test_sample3_multipleeval():
    test_file = 'bin/RStest3.txt'
    gamma = 1.4
    W_l = np.array([1.0, 0.0, 0.01])
    W_r = np.array([1.0, 0.0, 100.0])
    t = 0.035
    _test_sample_multipleeval(test_file, W_l, W_r, gamma, t)


def test_sample1_singleeval():
    test_file = 'bin/RStest1.txt'
    x0 = -0.1
    gamma = 1.4
    W_l = np.array([5.99924, 19.5975, 460.894])
    W_r = np.array([5.99242, -6.19633, 46.0950])
    t = 0.035
    _test_sample_singleeval(test_file, W_l, W_r, gamma, t, x0)


def test_sample1_multipleeval():
    test_file = 'bin/RStest1.txt'
    x0 = -0.1
    gamma = 1.4
    W_l = np.array([5.99924, 19.5975, 460.894])
    W_r = np.array([5.99242, -6.19633, 46.0950])
    t = 0.035
    _test_sample_multipleeval(test_file, W_l, W_r, gamma, t, x0)


def test_sample2_singleeval():
    test_file = 'bin/RStest2.txt'
    x0 = -0.2
    gamma = 1.4
    W_l = np.array([1.0, 0.75, 1.0])
    W_r = np.array([0.125, 0.0, 0.1])
    t = 0.2
    _test_sample_singleeval(test_file, W_l, W_r, gamma, t, x0)


def test_sample2_multipleeval():
    test_file = 'bin/RStest2.txt'
    x0 = -0.2
    gamma = 1.4
    W_l = np.array([1.0, 0.75, 1.0])
    W_r = np.array([0.125, 0.0, 0.1])
    t = 0.2
    _test_sample_multipleeval(test_file, W_l, W_r, gamma, t, x0)


def test_vacuum():
    x0 = 0.5
    gamma = 1.4
    W_l = np.array([1.0, -20, 1.0])
    W_r = np.array([1.0, 20, 1.0])
    riemann = Riemann(W_l, W_r, gamma, x0)
    try:
        riemann.update()
    except waves.VacuumCondition:
        return
    raise AssertionError


def test_almost_vacuum():
    W_l = np.array([1.0, -2.0, 0.4])
    W_r = np.array([1.0, 2.0, 0.4])
    gamma = 1.4
    riemann = Riemann(W_l, W_r, gamma)
    riemann.update()
    assert(riemann.p_star > 0.0)
    assert_equal(riemann.u_star, 0.0)
    assert_equal(np.array(riemann.left_wave.speed),
                 -np.array(riemann.right_wave.speed))
