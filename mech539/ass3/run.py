import numpy as np
from glob import glob
import os
from pprint import pformat
import shutil
import subprocess
params = dict(
   mach=0,
   solverID=1,
)
solvers = ['Regular', 'Line-Implicit']
grids = ['grid_coarse', 'grid_medium', 'grid_fine']
dm = 0.02
MACHS = np.arange(0.8, 0.9+dm, dm)

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

def move_dat(folder, prefix):
    if not os.path.exists(folder):
        os.mkdir(folder)
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

def run_all():
    print('Running all cases')
    for grid in grids:
        modify_grid(grid)
        for solverID, solver in enumerate(solvers, start=1):
            modify_param('solverID', solverID)
            for mach in MACHS:
                modify_param('mach', mach)
                print('Running %s' % prefix)
                print(run())
                move_dat('dat', prefix)

if __name__ == '__main__':
    run_all()
