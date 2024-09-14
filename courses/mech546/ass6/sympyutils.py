import sympy as sp
import IPython.display as disp

def mk_expr(expr, name, do_align=False):
    eq = '='
    if do_align:
        eq = '&' + eq

    return name + eq + sp.latex(expr)

def mk_align(body):
    return '\n'.join([r'\begin{aligned}', body, r'\end{aligned}'])

def show(expr, name):
    return disp.Math(mk_expr(expr, name))

def show_add(exprs, name):
    body = name + '=&' + r'\\ &+ '.join([sp.latex(expr) for expr in exprs])
    return disp.Math(mk_align(body))

def show_list(exprs, names):
    assert len(exprs) == len(names)
    body = r'\\'.join(
        [mk_expr(expr, name, True) for expr, name in zip(exprs, names)]
    )
    return disp.Math(mk_align(body))
