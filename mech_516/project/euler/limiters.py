"""Flux limiters"""

import numpy as np
import pdb


def limiter(W, name):
    name = name.lower()
    # Initialize return
    ret = np.zeros_like(W)

    if 'zero' == name:
        return ret

    if 'average' == name:
        delta_W_p = W[2:, :] - W[1:-1, :]
        delta_W_m = W[1:-1, :] - W[:-2, :]
        ret[1:-1, :] = 0.5*(delta_W_p + delta_W_m)
        # pdb.set_trace()
        return ret

    ########################
    ## Compute r, a, b, s ##
    ########################
    # r = U_{i+1} - U_{i}
    r = np.zeros_like(W)
    r[1:-1, :] = W[2:, :] - W[1:-1, :]
    r[0, :] = W[1, :] - W[0, :]
    r[-1, :] = W[-1, :] - W[-2, :]

    # a = U_{i} - U_{i-1}
    a = r[0:-1, :]
    abs_a = np.abs(a)

    # b = U_{i+1} - U_{i}
    b = r[1:, :]
    abs_b = np.abs(b)

    # s = (sgn(a) + sgn(b))/2
    s = (np.sign(a) + np.sign(b))/2

    if name in ('minmod', 'superbee', 'sweby'):
        if 'minmod' == name:
            beta = 1
        elif 'superbee' == name:
            beta = 2
        else:
            beta = 1.5

        # min1 = min(|a|, beta*|b|)
        min1 = np.minimum(abs_a, beta*abs_b)
        # min2 = min(beta*|a|, |b|)
        min2 = np.minimum(beta*abs_a, abs_b)

        ret[1:, :] = s*np.maximum(min1, min2)
        # pdb.set_trace()
        return ret

    if 'vanleer' == name:
        # raise NotImplementedError
        ret[1:, :] = 2*s*abs_a*abs_b/(abs_a + abs_b +
                                      np.ones_like(abs_a)*10**-16)
        return ret

    raise ValueError("Name %s not valid limiter." % name)
