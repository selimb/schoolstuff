from collections import namedtuple
import numpy as np
from glob import glob
import os
from pprint import pformat
import shutil
import subprocess


params = dict(
   mach=0,
   solverID=1,
   tol=2,
)
Solver = namedtuple('Solver', ['name', 'ID'])
SOLVERS = [
    Solver('GS', 1),
    Solver('LI', 2),
]
GRIDS = [
    'grid_coarse',
    'grid_medium',
    'grid_fine',
]
dm = 0.02
MACHS = np.arange(0.8, 0.9+dm, dm)

MACHINE = dict(
    folder='dat_machine',
    tol=1e-15,
    grids=[GRIDS[0]],
    solvers=SOLVERS,
    machs=MACHS
)
FINER = dict(
    folder='dat_finer',
    grids=GRIDS[1:],
    solvers=[SOLVERS[0]],
    machs=MACHS,
)

__all__ = []
for var in locals().keys():
    if var == var.upper():
        __all__.append(var)


def modify_param(name, val):
    filename = 'input.prm'
    vals = open(filename).read().split('\n')
    idx = params[name]
    vals[idx] = str(val)
    open(filename, 'w').write('\n'.join(vals))

def modify_grid(name):
    shutil.copy(name, 'grid')

def run():
    out = subprocess.check_output('./bin/solver.out')
    return out

def mkdir(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)

def move_dat(folder, prefix):
    datafiles = glob('*.dat')
    for oldpath in datafiles:
        filename = os.path.basename(oldpath)
        newfilename = '_'.join((prefix, filename))
        newpath = os.path.join(folder, newfilename)
        shutil.move(oldpath, newpath)

def get_time_taken():
    timefile = 't.dat'
    times = open(timefile).readlines()
    times = list(times)
    iters = len(times)
    abs_time_taken = float(times[-1].strip())
    time_per_iter = abs_time_taken/iters
    return iters, time_per_iter, abs_time_taken

def mk_prefix(grid, solver, mach):
    grid_suffix = grid.split('_')[1]
    prefix = '%s_%s_%.2f' % (grid_suffix, solver, mach)
    return prefix

def write_tol(folder, tol):
    fname = os.path.join(folder, 'tol')
    s = '%.2e' % tol
    with open(fname, 'w') as f:
        f.write(s)

def run_cases(folder, tol, grids, solvers, machs):
    mkdir(folder)
    modify_param('tol', tol)
    write_tol(folder, tol)
    print('Tolerance: %.2e' % tol)
    for grid in grids:
        modify_grid(grid)
        for solver in solvers:
            modify_param('solverID', solver.ID)
            for mach in machs:
                modify_param('mach', mach)
                prefix = mk_prefix(grid, solver.name, mach)
                print('Running %s' % prefix)
                print(run())
                move_dat(folder, prefix)

def run_machine():
    print('Running Machine Precision stuff')
    run_cases(**MACHINE)

def run_finer():
    print('Running finer grids (GS only)')
    run_cases(**FINER)

def run_all():
    print('Running all cases')
    run_machine()
    run_finer()

if __name__ == '__main__':
    run_all()
