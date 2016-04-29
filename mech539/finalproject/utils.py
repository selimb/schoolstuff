import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import subprocess
import re

def hide_ax(ax):
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

def mk_colors(l):
    num = len(l)
    if num == 1:
        raise ValueError('Really?')
    if num == 2:
        return ['b', 'g']
    if num == 3:
        return ['b', 'g', 'r']
    if num == 4:
        return ['m', 'b', 'g', 'r']
    if num == 5:
        return ['m', 'b', 'g', 'darkorange', 'r']
    if num > 5:
        cm_subsection = np.linspace(0, 1, num)
        return [ cm.rainbow(x) for x in cm_subsection ]

FLOWPRM_FILE = 'flow.prm'
PTOT_IN = 2117

class SolverError(Exception):
    pass

def read_params():
    return open(FLOWPRM_FILE).read()

def modify_param(param, new):
    flowprm = read_params()
    open('%s.bak' % FLOWPRM_FILE, 'w').write(flowprm)

    pattern = r'( params%{name}=).*(,|\/)\n'.format(name=param)
    r = re.search(pattern, flowprm)
    if r is None:
        print('Could not find a match for \n%s' % pattern)
        return
    replace = r'\g<1>{new}\g<2>\n'.format(new=str(new))
    newflowprm = re.sub(
        pattern=pattern,
        repl=replace,
        string=flowprm)
    open(FLOWPRM_FILE, 'w').write(newflowprm)

def run():
    cmd = './bin/solver.out'
    proc = subprocess.Popen(cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = proc.stdout.read().decode()
    err = proc.stderr.read().decode()
    if err:
        raise SolverError(err)
    return out

def read_data():
    state_file = 'state.csv'
    with open(state_file) as f:
        cols = f.readline().split()
    state = np.genfromtxt(state_file, skip_header=1, names=cols)
    residuals = np.genfromtxt('residuals.csv')
    time = np.genfromtxt('time.csv')
    return state, residuals, time

def calc_mach(state):
    return state['u']/state['c']

def plot_mach(state, ax, title_args=None, **kwargs):
    d = dict(linestyle='solid', marker='.', color='b')
    d.update(**kwargs)
    mach = calc_mach(state)
    ax.plot(state['x'], mach, **d)
    ax.set_yticks([1.0], minor=True)
    ax.grid(which='minor', axis='y')
    ax.set_ylabel('$M$')
    ax.set_xlabel('$x$')
    title = 'Mach Number distribution'
    if title_args:
        title += title_args
    ax.set_title(title)

def plot_residual(residuals, ax, do_labels=True, **kwargs):
    d = dict(linestyle='solid', color='b')
    d.update(**kwargs)
    ax.semilogy(residuals, **d)
    if do_labels:
        ax.set_ylabel('Residuals')
        ax.set_xlabel('Iterations')
        title = 'Convergence of density residual'
        ax.set_title(title)

def plot_pressure(state, ax, title_args=None, **kwargs):
    d = dict(linestyle='solid', marker='.', color='b')
    d.update(**kwargs)
    ax.plot(state['x'], state['p']/PTOT_IN, **d)
    ax.set_ylabel('$p/p_t$')
    ax.set_xlabel('$x$')
    title = 'Pressure distribution'
    if title_args:
        title += title_args
    ax.set_title(title)

def plot_residual_time(residuals, time, ax, do_labels=True, **kwargs):
    d = dict(linestyle='solid', color='b')
    d.update(**kwargs)
    ax.semilogy(time, residuals, **d)
    if do_labels:
        ax.set_ylabel('Residuals')
        ax.set_xlabel('CPU time [s]')
#     ax.set_ylim(0, 1)
        title = 'Convergence of density residual'
        ax.set_title(title)

