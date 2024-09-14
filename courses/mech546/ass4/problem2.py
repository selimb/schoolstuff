from IPython.display import clear_output
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin, cos
from pymatbridge import Matlab
get_ipython().magic(u'matplotlib inline')
execfile('../../matplotlibrc.py')

# Run Matlab code and fetch relevant data
mlab = Matlab()
mlab.start()
results = mlab.run_code(open('fem1d.m').read())
K = mlab.get_variable('K')
U = mlab.get_variable('U')
nodeLocs = mlab.get_variable('nodeLocs')
mlab.stop()
clear_output()
print('K')
print(K)
print('U')
print(U)

x = nodeLocs[:,0]

xex = np.linspace(0, 1);
w = np.sqrt(2);
uex = (cos(w) - 1)*sin(w*xex)/(2.0*sin(w)) - 0.5*cos(w*xex) + 0.5

fig, ax = plt.subplots(figsize=(14, 10))
ax.plot(xex, uex, 'k-', label='Exact')
ax.plot(x, U, 'b-o', label='FEM')
ax.set_xlabel('Position $x$')
ax.set_ylabel('Displacement $u(x)$')
ax.set_title('Comparison of analytical solution and 5-element FEM solution.')
ax.legend()
plt.show()

