import os
import shutil

FOLDER=''
files = ['airfoilpressure.dat', 'cl-cd.dat', 'transition.dat', 'upper_cf.dat', 'lower_cf.dat']

def cp(prefix):
    folder = os.path.join('dat', FOLDER)
    print(folder)
    if not os.path.exists(folder):
        os.mkdir(folder)
    for f in files:
        filename = '_'.join([prefix, f])
        filepath = os.path.join(folder, filename)
        og = os.path.join('pablo', f)
        shutil.move(og, filepath)
