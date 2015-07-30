from nose.tools import assert_equals

from ..riemann.utils import compute_c, compute_c_from_primitives


def setup():
    print("SETUP!")


def teardown():
    print("TEAR DOWN!")


def test_basic():
    print("I RAN!")
