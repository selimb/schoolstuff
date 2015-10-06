"""Calculate values of u at next time step.
u: values of u at current time step
We don't update the values of u at i=0 and i=-1
"""

import numpy as np

def lax_friedrichs(u, a, dt, dx):
    u_new = np.empty_like(u)
    for i in range(1,len(u_new)-1):
        first = dt/dx*a*(u[i+1] - u[i-1])
        second = u[i+1] - 2*u[i] + u[i-1]
        u_new[i] = u[i] - first/2.0 + second/2.0
    u_new[0] = u[0]
    u_new[-1] = u[-1]
    return u_new

def lax_wendroff(u, a, dt, dx):
    u_new = np.empty_like(u)
    for i in range(1, len(u_new)-1):
        first = dt/dx*a*(u[i+1] - u[i-1])/2
        second = (dt/dx*a)**2*(u[i+1] - 2*u[i] + u[i-1])/2
        u_new[i] = u[i] - first + second
    u_new[0] = u[0]
    u_new[-1] = u[-1]
    return u_new
