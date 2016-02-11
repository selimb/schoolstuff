import numpy as np
from glob import glob
import os
import shutil
import subprocess

params = dict(
   nx=0,
   solverID=1,
   relax=2,
)
solvers = ['Jacobi', 'Gauss', 'SOR']
NXs = [100, 200, 400]

def modify_param(name, val):
    filename = 'input.prm'
    vals = open(filename).read().split('\n')
    idx = params[name]
    vals[idx] = str(val)
    open(filename, 'w').write('\n'.join(vals))

def run():
    out = subprocess.check_output('./bin/solver.out')
    return out

def move_dat(prefix):
    datafiles = glob('*.dat')
    for oldpath in datafiles:
        filename = os.path.basename(oldpath)
        newfilename = '_'.join((prefix, filename))
        newpath = os.path.join('dat', newfilename)
        shutil.move(oldpath, newpath)


def run_all():
    print('Running all cases')
    modify_param('relax', 1.5)
    for solverID, solver in enumerate(solvers, start=1):
        for nx in NXs:
            modify_param('nx', nx)
            modify_param('solverID', solverID)
            prefix = '%s_%i' % (solver, nx)
            print('Running %s' % prefix)
            print(run())
            move_dat(prefix)
    print('Done')
if __name__ == '__main__':
    run_all()
