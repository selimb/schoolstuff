import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import scipy.integrate
import scipy.optimize

class Collector(list):
    def add(self, f):
        self.append(f)
        return f

solvers = Collector()

@solvers.add
def forward_euler(func, y0, t):
    steps = len(t)
    y = np.empty((len(t), len(y0)))
    y[0] = y0
    for i in range(1, steps):
        dt = t[i] - t[i-1]
        y[i] = y[i-1] + dt*func(y[i-1], t[i-1])
    return y

@solvers.add
def backward_euler(func, y0, t):
    steps = len(t)
    y = np.empty((len(t), len(y0)))
    y[0] = y0
    for i in range(1, steps):
        dt = t[i] - t[i-1]
        y[i] = scipy.optimize.fsolve(
            func=lambda x: x - y[i-1] - dt*func(x, t[i]),
            x0=y[i-1]
        )
    return y

@solvers.add
def crank(func, y0, t):
    steps = len(t)
    y = np.empty((len(t), len(y0)))
    y[0] = y0
    for i in range(1, steps):
        dt = t[i] - t[i-1]
        y[i] = scipy.optimize.fsolve(
            func=lambda x: x - y[i-1] - 0.5*dt*(func(x, t[i]) + func(y[i-1], t[i-1])),
            x0=y[i-1]
        )
    return y

def pi_plot(ax, pi1, pi2):
    a = np.pi*pi1
    b = np.pi*pi2
    ax.set_xlim(a, b)
    ticks = np.arange(a, b + np.pi/2, np.pi/2)
    make_label = lambda p: '$%.1f\pi$' % (p/np.pi)
    labels = list(map(make_label, ticks))
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    minor_locator = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.grid(which='major')

def f(x, t):
    return np.array([
        -x[1],
         x[0],
    ])


times = np.linspace(0, np.pi*5, 100)
x0 = np.array([1, 0])
y = scipy.integrate.odeint(f, x0, times)
fig, ax = plt.subplots(figsize=(12, 8))
for solver in solvers:
    name = solver.func_name
    y = solver(f, x0, times)
    ax.plot(times, y[:, 1], label=name)
ax.plot(times, np.sin(times), 'k--', label='analytical')
pi_plot(ax, 0, 5)
ax.legend(fancybox=True, loc='best')
ax.set_xlabel('$t$')
ax.set_ylabel('x(t)')
ax.set_title('Comparison of various ODE solvers')
plt.show()
