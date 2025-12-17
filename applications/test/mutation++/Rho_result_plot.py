import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("twoT_energy_based.dat")

t = data[:,0]

rho_N2 = data[:,9]
rho_O2 = data[:,10]

plt.plot(t, rho_N2, label="N2")
plt.plot(t, rho_O2, label="O2")

plt.xscale("log")
plt.xlabel("Time [s]")
plt.ylabel(r"Species density [kg/m$^3$]")
plt.legend()
plt.grid(True)
plt.show()
