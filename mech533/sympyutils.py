import sympy as sp
import IPython.display as disp

def show(expr, name):
    return disp.Math(name + '=' + sp.latex(expr))
