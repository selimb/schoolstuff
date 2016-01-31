import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

TMAX = 10
A = 0.5
GRID_LENGTH = 40

def signal(x):
    # tanh = np.tanh(250.0*(x - 20.0))
    # return 0.5*(1.0 + tanh)
    w = 4
    xrange = np.logical_and(x >= 0, x <= w*2*np.pi)
    u = np.zeros_like(x)
    u[xrange] = np.sin(x[xrange]/w)
    return u

def u_exact(x, t):
    return signal(x - A*t)

# Functions to solve flow
def upwind(um, ui, up, dt, dx, u):
    u_new = u.copy()
    u_new[1:-1] = ui - A*(dt/dx)*(ui - um)
    return u_new

def lax_friedrichs(um, ui, up, dt, dx, u):
    u_new = u.copy()
    u_new[1:-1] = 0.5*(up + um) - A*0.5*(dt/dx)*(up - um)
    return u_new

def lax_wendroff(um, ui, up, dt, dx, u):
    u_new = u.copy()
    u_new[1:-1] = (ui
            - 0.5*(dt/dx)*A*(up - um)
            + 0.5*(dt/dx*A)**2*(up - 2*ui + um))
    return u_new

def time_step(u, dt, dx, scheme):
    u = u.copy()
    ui = u[1:-1]
    um = u[:-2]
    up = u[2:]
    u_new = scheme(um, ui, up, dt, dx, u)
    return u_new

def solve(scheme, dt, dx):
    x = np.arange(0, GRID_LENGTH + dx, dx, dtype=float)
    tsteps = int(TMAX/dt)
    u = signal(x)
    t = dt
    for i in range(tsteps - 1):
        t += dt
        u = time_step(u, dt, dx, scheme)
    return x, u

def mccormack(dt, dx):
    x = np.arange(0, GRID_LENGTH + dx, dx, dtype=float)
    tsteps = int(TMAX/dt)
    u = signal(x)
    for n in range(tsteps - 1):
        u_bar = u.copy()
        for i in range(1, len(x) - 1):
            u_bar[i] = u[i] - A*(dt/dx)*(u[i+1] - u[i])
        u_new = u.copy()
        for i in range(1, len(x) - 1):
            u_new[i] = 0.5*(u[i] + u_bar[i] - A*(dt/dx)*(u_bar[i] - u_bar[i-1]))
        u = u_new
    return x, u

def leap_frog(dt, dx):
    def timestep(u_jm1, u_j):
        u_new = u_j.copy()
        u_new[1:-1] = (
                u_jm1[1:-1]
                - A*(dt/dx)*(u_j[2:] - u_j[:-2]))
        return u_new

    x = np.arange(0, GRID_LENGTH + dx, dx, dtype=float)
    tsteps = int(TMAX/dt)
    u = np.zeros((2, len(x)))
    u[0] = signal(x)
    u[1] = timestep(u[0], u[0])
    for i in range(1, tsteps - 1):
        u_new = timestep(u[0], u[1])
        u[0] = u[1].copy()
        u[1] = u_new.copy()
    return x, u[1]

def plot_exact(ax):
    dx = 0.001
    x = np.arange(0, GRID_LENGTH + dx, dx)
    ux = u_exact(x, TMAX)
    ax.plot(x, ux, 'k-', label='Exact')

def plot_options(ax, dt, dx):
    title = r'$\Delta x$ = %.1f, $\Delta t$ = %.1f, CFL = %.2f, t = %i' % (dx, dt, A*dt/dx, TMAX)
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.set_title('%s' % title)
    ax.set_ylim(-1.0, 1.5)
    ax.legend()

def q1(dt, dx, xlim=None):
    fig, ax = plt.subplots(figsize=(20, 15))
    x, u_upwind = solve(upwind, dt, dx) ; ax.plot(x, u_upwind, 'r--', label='Upwind')
    _, u_lax = solve(lax_friedrichs, dt, dx) ; ax.plot(x, u_lax, 'b-.', label='Lax-Friedrichs')
    _, u_lw = solve(lax_wendroff, dt, dx) ; ax.plot(x, u_lw, 'g:', label='Lax-Wendroff')
    _, u_mc = mccormack(dt, dx) ; ax.plot(x, u_mc, '-', label='McCormack')
    x, u_lf = leap_frog(dt, dx) ; ax.plot(x, u_lf, '-', label='Leap Frog')
    plot_exact(ax)
    plot_options(ax, dx=dx, dt=dt)
    if xlim is not None:
        ax.set_xlim(xlim)

    plt.show()

MAX_PTS = 4001
DXq2 = [1.0, 0.5, 0.1, 0.01, 0.001]
def calcs_q2():
    CFL = 0.75
    for i, dx in enumerate(DXq2):
        dt = CFL*dx/A
        x, u_up = solve(upwind, dt, dx)
        _, u_lax = solve(lax_wendroff, dt, dx)
        suffix = '_%s' % dx
        np.save(''.join(['x', suffix]), x)
        np.save(''.join(['u_up', suffix]), u_up)
        np.save(''.join(['u_lax', suffix]), u_lax)

def q2_1():
    def compute_pts(pts):
        if pts <= MAX_PTS:
            return 1
        else:
            return pts/MAX_PTS
    CFL = 0.75
    fig, ax = plt.subplots(figsize=(20, 15))
    colors = ['b', 'g', 'r', 'c', 'm', 'darkorange']
    line_handles = []
    for i, dx in enumerate(DXq2):
        suffix = '_%s.npy' % dx
        x = np.load(''.join(['x', suffix]))
        u_up = np.load(''.join(['u_up', suffix]))
        u_lax = np.load(''.join(['u_lax', suffix]))
        pts = compute_pts(len(x))
        color = colors[i]
        ax.plot(x, u_up, color=color, marker='o', markevery=pts)
        ax.plot(x, u_lax, color=color, marker='s', markevery=pts)
        line = mlines.Line2D([], [], color=color, label='$\Delta x$ = %s' % dx)
        line_handles.append(line)
    plot_exact(ax)
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.set_ylim(-1.0, 1.5)
    # ax.set_xlim(24, 26)
    title = '''Grid/Timestep Refinement with constant CFL = %s
    Upwind: Circles -- Lax: Squares
    ''' % (CFL)
    ax.set_title(title)
    ax.legend(handles=line_handles)

def q2_2():
    def calc_err(u, ux):
        n = len(u)
        err = u - ux
        return np.linalg.norm(err)/(np.sqrt(n))

    fig, ax = plt.subplots(figsize=(20, 15))
    err_up = []
    err_lax = []
    for i, dx in enumerate(DXq2):
        suffix = '_%s.npy' % dx
        x = np.load(''.join(['x', suffix]))
        u_up = np.load(''.join(['u_up', suffix]))
        u_lax = np.load(''.join(['u_lax', suffix]))
        ux = u_exact(x, TMAX)
        err_up.append(calc_err(u_up, ux))
        err_lax.append(calc_err(u_lax, ux))

    ax.semilogy(DXq2, err_up, label='Upwind')
    ax.semilogy(DXq2, err_lax, label='Lax')
    ax.set_ylabel(r'$\log\Vert e\Vert$', fontsize='x-large')
    ax.set_xlabel(r'$\Delta x$', fontsize='x-large')
    ax.legend()
