import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------
# Load data from the C++ solver
# -------------------------------------------------
filename = "twoT_energy_based.dat"

data = np.loadtxt(filename)

time = data[:, 0]
Ttr  = data[:, 1]
Tv   = data[:, 2]

# -------------------------------------------------
# Plot: Two temperatures (paper-style)
# -------------------------------------------------
plt.figure(figsize=(7, 5))

plt.plot(time, Ttr, 'k-', linewidth=2, label=r'$T_{tr}$')
plt.plot(time, Tv,  'r--', linewidth=2, label=r'$T_v$')

plt.xscale('log')
plt.xlabel('Time, $t$ [s]')
plt.ylabel('Temperature [K]')
plt.title('Two-Temperature Thermal Relaxation')
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend()

plt.tight_layout()
plt.savefig("figure_2T_relaxation.png", dpi=300)
plt.show()

print("Saved: figure_2T_relaxation.png")
