import pandas as pd
import matplotlib.pyplot as plt

# 1. Load Data (Updated filename for the 2-Temperature simulation)
filename = 'air_2T_data.csv'

try:
    df = pd.read_csv(filename)
except FileNotFoundError:
    print(f"Error: '{filename}' not found. Run the 2-Temperature C++ code first.")
    exit()

# 2. Constants for Normalization (Molar Masses kg/kmol)
Mw = {
    'rho_N2': 28.0134,
    'rho_O2': 31.9988,
    'rho_NO': 30.0061,
    'rho_N':  14.0067,
    'rho_O':  15.9994
}

# Calculate Molar Concentrations [kmol/m3] and Normalize
molar_conc = pd.DataFrame()
for col, mass in Mw.items():
    if col in df.columns:
        molar_conc[col] = df[col] / mass

# Total Initial Molar Concentration (n0)
n0 = molar_conc.iloc[0].sum()
norm_density = molar_conc / n0

# 3. Create Plots (Replicating Figure 9 style but with 2 Temperatures)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# --- Subplot A: Temperatures vs Time ---
# Now plotting BOTH T and Tv
ax1.plot(df['time'], df['T'], 'k-', linewidth=2, label='$T$ (Trans-Rot)')
ax1.plot(df['time'], df['Tv'], 'r--', linewidth=2, label='$T_v$ (Vib-Elec)')

ax1.set_xscale('log')
ax1.set_xlabel('Time, t [s]', fontsize=12)
ax1.set_ylabel('Temperature [K]', fontsize=12)
ax1.set_title('Thermal Relaxation ($T$ vs $T_v$)', fontsize=14)
ax1.set_xlim(1e-9, 1e-4) # Adjusted to 1e-4 as per 2T sim end time
ax1.set_ylim(0, 12000)
ax1.grid(True, which="both", linestyle='--', alpha=0.6)
ax1.legend(loc='best')

# --- Subplot B: Normalized Number Density vs Time ---
# Markers/Lines for Species
ax1_lines = []
if 'rho_N2' in norm_density:
    ax2.plot(df['time'], norm_density['rho_N2'], 'k-', label='$N_2$')
if 'rho_O2' in norm_density:
    ax2.plot(df['time'], norm_density['rho_O2'], 'k--', label='$O_2$')
if 'rho_NO' in norm_density:
    ax2.plot(df['time'], norm_density['rho_NO'], 'r-.', label='NO')
if 'rho_N' in norm_density:
    ax2.plot(df['time'], norm_density['rho_N'], 'g:', linewidth=2, label='N')
if 'rho_O' in norm_density:
    ax2.plot(df['time'], norm_density['rho_O'], 'b:', linewidth=2, label='O')

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Time, t [s]', fontsize=12)
ax2.set_ylabel('Normalised number density, $n/n^0$', fontsize=12)
ax2.set_title('Species Evolution', fontsize=14)
ax2.set_xlim(1e-9, 1e-4)
ax2.set_ylim(1e-5, 1.2)
ax2.grid(True, which="both", linestyle='--', alpha=0.6)
ax2.legend(loc='lower left')

# 4. Save and Show
plt.tight_layout()
plt.savefig('figure_2T_relaxation.png', dpi=300)
print("Plot saved as 'figure_2T_relaxation.png'")
plt.show()