import sympy as sp
import numpy as np
from utils import printsp

x, y = sp.symbols('x y')
_F = sp.Matrix([
    3*x**2 - y**2,
    3*x*y**2 - x**3 - 1
])
_v = sp.Matrix([x, y])
_J = _F.jacobian(_v)

def sub(A, X):
    return np.array(A.subs(x, X[0]).subs(y, X[1]).tolist(), dtype=np.float32)

def F(X=None):
    if X is None:
        return printsp(_F, 'F')
    return sub(_F, X)

def J(X=None):
    if X is None:
        return printsp(_J, 'J')
    return sub(_J, X)
