import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

TMAX = 10
A = 0.5
GRID_LENGTH = 40

def signal(x):
    tanh = np.tanh(250.0*(x - 20.0))
    return 0.5*(1.0 + tanh)
    # w = 4
    # xrange = np.logical_and(x >= 0, x <= w*2*np.pi)
    # u = np.zeros_like(x)
    # u[xrange] = np.sin(x[xrange]/w)
    # return u

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
    t = 0
    for i in range(tsteps):
        t += dt
        u = time_step(u, dt, dx, scheme)
    return x, u

def mccormack(dt, dx):
    x = np.arange(0, GRID_LENGTH + dx, dx, dtype=float)
    tsteps = int(TMAX/dt)
    u = signal(x)
    for n in range(tsteps):
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
    t = dt
    for i in range(tsteps - 1):
        t += dt
        u_new = timestep(u[0], u[1])
        u[0] = u[1].copy()
        u[1] = u_new.copy()
    return x, u[1]

def plot_exact(ax, **kwargs):
    dx = 0.001
    x = np.arange(0, GRID_LENGTH + dx, dx)
    ux = u_exact(x, TMAX)
    ax.plot(x, ux, 'k-', label='Exact', **kwargs)

def plot_options(ax, dt, dx):
    title = r'$\Delta x$ = %.1f, $\Delta t$ = %.1f, CFL = %.2f, t = %i' % (dx, dt, A*dt/dx, TMAX)
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.set_title('%s' % title)
    ax.set_ylim(-0.5, 1.1)
    ax.legend()

def q1(dt, dx, xlim=None):
    fig, ax = plt.subplots(figsize=(20, 15))
    x, u_upwind = solve(upwind, dt, dx) ; ax.plot(x, u_upwind, 'm', label='Upwind')
    _, u_lax = solve(lax_friedrichs, dt, dx) ; ax.plot(x, u_lax, c='darkorange',
                                                       label='Lax-Friedrichs')
    _, u_lw = solve(lax_wendroff, dt, dx) ; ax.plot(x, u_lw, ':', c='g', lw=10,
                                                    label='Lax-Wendroff')
    _, u_mc = mccormack(dt, dx) ; ax.plot(x, u_mc, 'r', label='McCormack')
    x, u_lf = leap_frog(dt, dx) ; ax.plot(x, u_lf, label='Leap Frog')
    plot_exact(ax)
    plot_options(ax, dx=dx, dt=dt)
    if xlim is not None:
        ax.set_xlim(xlim)

    plt.show()

MAX_PTS = 4001
DXq2 = [1.0, 0.75, 0.5, 0.25,
        0.1, 0.075, 0.05, 0.025,
        0.01, 0.0075, 0.005, 0.0025,
        0.001, 0.00075, 0.0005, 0.00025,
        0.0001]
# DXq2 = [0.0025,]
# DXq2 = [1.0, 0.5, 0.1, 0.01, 0.001]
def calcs_q2():
    CFL = 0.75
    for dx in DXq2:
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
    for i, dx in enumerate(DXq2[::4]):
        suffix = '_%s.npy' % dx
        x = np.load(''.join(['x', suffix]))
        u_up = np.load(''.join(['u_up', suffix]))
        u_lax = np.load(''.join(['u_lax', suffix]))
        pts = compute_pts(len(x))
        color = colors[i]
        ax.plot(x, u_up, '--', color=color)
        ax.plot(x, u_lax, '-', color=color)
        line = mlines.Line2D([], [], color=color, label='$\Delta x$ = %s' % dx)
        line_handles.append(line)
    plot_exact(ax, linestyle='dotted')
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.set_ylim(-0.5, 1.1)
    # ax.set_xlim(24, 26)
    ax.set_xlim(24.5, 25.5)
    title = '''Grid/Timestep Refinement with constant CFL = %s
    Upwind: Dashed -- Lax-Wendroff: Solid
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

    nx = [GRID_LENGTH/dx for dx in DXq2]
    ax.loglog(nx[6:], err_up[6:], '-o', label='Upwind')
    ax.loglog(nx[6:], err_lax[6:], '-o', label='Lax-Wendroff')
    ax.set_ylabel(r'Error $e$', fontsize='x-large')
    ax.set_xlabel(r'Number of grid points', fontsize='x-large')
    ax.legend()

def solve_with_history(scheme, dt, dx):
    x = np.arange(0, GRID_LENGTH + dx, dx, dtype=float)
    tsteps = int(TMAX/dt)
    u = np.zeros((tsteps + 1, len(x)))
    u[0] = signal(x)
    t = 0
    for i in range(tsteps):
        t += dt
        u[i+1] = time_step(u[i], dt, dx, scheme)
    return x, u

def q4(dx):
    def cmp_amp(u):
        offset = 0.5
        return np.abs(u - offset).max(axis=1)

    all_CFL = [0.75, 1.0, 1.001, 1.005, 1.01, 1.02]
    fig, ax = plt.subplots(figsize=(20, 15))
    colors = ['b', 'g', 'r', 'c', 'm', 'darkorange']
    line_handles = []
    for i, CFL in enumerate(all_CFL):
        dt = CFL*dx/A
        x, u_up = solve_with_history(upwind, dt, dx)
        _, u_lax = solve_with_history(lax_wendroff, dt, dx)
        amp_up = cmp_amp(u_up)
        amp_lax = cmp_amp(u_lax)
        color = colors[i]
        ax.plot(amp_up, '--', color=color)
        ax.plot(amp_lax, '-', color=color)
        line = mlines.Line2D([], [], color=color, label='CFL = %s' % CFL)
        line_handles.append(line)
    ax.set_xlabel('Time step')
    ax.set_ylabel('Wave amplitude')
    title = '''CFL increment with constant $\Delta x$ = %s
    Upwind: Solid -- Lax: Dashed
    ''' % (dx)
    ax.set_ylim(0.4, 1.6)
    ticks = list(ax.get_yticks())
    ticks.append(0.5)
    ax.set_yticks(ticks)
    ax.set_title(title)
    ax.legend(handles=line_handles)
