import numpy as np
from glob import glob
import os
from pprint import pformat
import shutil
import subprocess

params = dict(
   nx=0,
   solverID=1,
   relax=2,
   do_cond=3,
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

def compile(olevel):
    subprocess.call("cd src; make compile olevel=%s" % olevel, shell=True)

def time_analysis():
    modify_param('do_cond', 0)
    olevels = ['-g', '-O1', '-O2', '-O3']
    solvers = [1, 2] # Jacobi and Gauss
    NXs = [100,200,400]
    times = {}
    for olevel in olevels:
        times[olevel] = {}
        print('Compiling with %s' % olevel)
        compile(olevel)
        for solver in solvers:
            times[olevel][solver] = {}
            for nx in NXs:
                modify_param('solverID', solver)
                modify_param('nx', nx)
                prefix = '%s_%i' % (solver, nx)
                print('Running %s' % prefix)
                print(run())
                times[olevel][solver][nx] = get_time_taken()
    open('time_data.py', 'w').write(pformat(times))

def run_all():
    modify_param('do_cond', 0)
    print('Running all cases')
    modify_param('relax', 1.5)
    for solverID, solver in enumerate(solvers, start=1):
        for nx in NXs:
            modify_param('nx', nx)
            modify_param('solverID', solverID)
            prefix = '%s_%i' % (solver, nx)
            print('Running %s' % prefix)
            print(run())
            move_dat('dat', prefix)
    print('Done')

def run_sor():
    modify_param('do_cond', 0)
    modify_param('solverID', 3)
    all_relax = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 1.95]
    for relax in all_relax:
        modify_param('relax', relax)
        for nx in NXs:
            modify_param('nx', nx)
            prefix = '%.2f_%i' % (relax, nx)
            print('Running %s' % prefix)
            print(run())
            move_dat('sordat', prefix)

if __name__ == '__main__':
    # time_analysis()
    # run_all()
    run_sor()
