import numpy as np
import matplotlib.pyplot as plt

class Solver(object):

    a = 2
    t_max = 100
    u0 = 0.5

    def __init__(self, scheme, disturbance):
        self.scheme = scheme
        self.disturbance = disturbance # Disturbance function

    def solve(self, courant, dt=1, x_lims=[-20.0, 400.0]):
        """scheme = func(u, a, dt, dx)"""
        self.courant = courant
        self.dt = dt
        self.x_lims = x_lims
        self._initialize_grid()
        for n in range(len(self.t)-1):
            self.u[n+1] = self.scheme(self.u[n], self.a,
                                      self.dt, self.dx)

    def _initialize_grid(self):
        self.dx = self.a*self.dt/self.courant
        x = np.arange(self.x_lims[0], self.x_lims[1],
                      self.dx)
        t = np.arange(0, self.t_max+self.dt, self.dt)

        x_begin = np.logical_and(x<=20, x>=0)
        u = np.zeros((len(t), len(x)))
        u[0] = self.u0
        u[0, x_begin] = self.disturbance(x[x_begin])
        self.u, self.x, self.t = u, x, t

    def plot(self, time, left=None, right=None, ax=None, *args, **kwargs):
        if ax is None:
            ax = plt.gca()
        if left is None:
            left = 0
        else:
            left = (self.x > left).argmax()
        if right is None:
            right = -1
        else:
            right = (self.x > right).argmax()
        ax.plot(self.x[left:right], self.u[time, left:right], 'b.',
                markersize=10, *args, **kwargs)
