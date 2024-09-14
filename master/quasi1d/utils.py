import numpy as np
import matplotlib.pyplot as plt
import subprocess
import re

FLOWPRM_FILE = 'flow.prm'

def read_params():
    return open(FLOWPRM_FILE).read()

def modify_param(param, new):
    flowprm = read_params()
    open('%s.bak' % FLOWPRM_FILE, 'w').write(flowprm)

    param_line = r' params%{0}={{0}},\n'.format(param)
    pattern = param_line.format('.*')
    if re.search(pattern, flowprm) is None:
        print('Could not find a match for \n%s' % pattern)
        return

    newflowprm = re.sub(
        pattern=pattern.format('.*'),
        repl=param_line.format(str(new)),
        string=flowprm)
    open(FLOWPRM_FILE, 'w').write(newflowprm)

def run():
    cmd = './q1d.out'
    proc = subprocess.Popen(cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = proc.stdout.read().decode()
    err = proc.stderr.read().decode()
    if err:
        return err
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

PLOT_KWARG = dict(
    linestyle='solid',
    marker='.'
)

def plot_mach(state, ax, title_args=None, **kwargs):
    PLOT_KWARG.update(**kwargs)
    mach = calc_mach(state)
    ax.plot(state['x'], mach, **PLOT_KWARG)
    ax.set_ylabel('$M$')
    ax.set_xlabel('$x$')
    title = 'Mach Number distribution'
    if title_args:
        title += title_args
    ax.set_title(title)

def plot_residual(residuals, ax, title_args=None, **kwargs):
    PLOT_KWARG.update(**kwargs)
    ax.semilogy(residuals, **PLOT_KWARG)
    ax.set_ylabel('Residuals')
    ax.set_xlabel('Iterations')
#     ax.set_ylim(0, 1)
    title = 'Convergence of density residual'
    if title_args:
        title += title_args
    ax.set_title(title)

def plot_pressure(state, ax, title_args=None, **kwargs):
    PLOT_KWARG.update(**kwargs)
    ax.plot(state['x'], state['p'], **PLOT_KWARG)
    ax.set_ylabel('$P$')
    ax.set_xlabel('$x$')
    title = 'Pressure distribution'
    if title_args:
        title += title_args
    ax.set_title(title)

def plot_residual_time(residuals, time, ax, title_args=None, **kwargs):
    PLOT_KWARG.update(**kwargs)
    ax.semilogy(time, residuals, **PLOT_KWARG)
    ax.set_ylabel('Residuals')
    ax.set_xlabel('CPU time [s]')
#     ax.set_ylim(0, 1)
    title = 'Convergence of density residual'
    if title_args:
        title += title_args
    ax.set_title(title)

