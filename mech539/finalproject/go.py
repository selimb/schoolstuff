import matplotlib.pyplot as plt
from utils import *

out = run()
print out
state, resi, time = read_data()
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(10, 15))
plot_mach(state, ax1)
plot_pressure(state, ax2)
plot_residual(resi, ax3)
plt.show()
