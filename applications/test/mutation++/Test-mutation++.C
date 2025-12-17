#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mutation++.h>
#include <algorithm> // For std::max, std::min

// Helper: Ensure density never goes negative
static inline double pos(double x) { return (x < 1.0e-30) ? 1.0e-30 : x; }

struct CellState
{
    double rho;
    double T;
    double Tv;
    std::vector<double> rhoY;
};

void chemistry_relaxation_step(Mutation::Mixture &mix, CellState &cell, double dt)
{
    // 1. Calculate current total internal energy (e_target)
    // We must conserve this energy (Adiabatic process)
    double Tvec[] = {cell.T, cell.Tv};
    mix.setState(cell.rhoY.data(), Tvec, 2);
    double e_target = mix.mixtureEnergyMass(); // J/kg

    // 2. Get Production Rates (wdot)
    std::vector<double> wdot(mix.nSpecies());
    mix.netProductionRates(wdot.data());

    // 3. Update Densities (Explicit Euler)
    double rho_sum = 0.0;
    for (int i = 0; i < mix.nSpecies(); ++i)
    {
        // Update density
        double change = wdot[i] * dt;

        // Safety: Don't consume more than we have (prevents negative mass instabilities)
        if (change < 0.0 && std::abs(change) > cell.rhoY[i])
        {
            change = -0.99 * cell.rhoY[i];
        }

        cell.rhoY[i] += change;
        cell.rhoY[i] = pos(cell.rhoY[i]);
        rho_sum += cell.rhoY[i];
    }

    // Correct total density drift
    double correction = cell.rho / rho_sum;
    for (int i = 0; i < mix.nSpecies(); ++i)
        cell.rhoY[i] *= correction;

    // 4. Solve for Temperature (Newton-Raphson with Clamping)
    double T_new = cell.T;

    for (int iter = 0; iter < 20; ++iter)
    {
        // Assume 1-Temperature model for this validation
        double T_dual[] = {T_new, T_new};
        mix.setState(cell.rhoY.data(), T_dual, 2);

        double e_current = mix.mixtureEnergyMass();
        double Cv = mix.mixtureFrozenCvMass(); // Specific heat

        double error = e_target - e_current;

        if (std::abs(error) < 1.0)
            break; // Converged

        // Calculate step
        double dT = error / Cv;

        // --- SAFETY CLAMPING ---
        // 1. Limit step size to avoid wild oscillations
        double max_step = 1000.0;
        if (dT > max_step)
            dT = max_step;
        if (dT < -max_step)
            dT = -max_step;

        // 2. Apply step
        T_new += dT;

        // 3. Hard limits on Temperature
        T_new = std::max(300.0, std::min(T_new, 15000.0));
    }

    cell.T = T_new;
    cell.Tv = T_new;
}

int main()
{
    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    // Indices
    int iN2 = mix.speciesIndex("N2");
    int iO2 = mix.speciesIndex("O2");
    int iNO = mix.speciesIndex("NO");
    int iN = mix.speciesIndex("N");
    int iO = mix.speciesIndex("O");

    // Initial Conditions
    double P_init = 0.063 * 101325.0; // Pa
    double T_init = 10000.0;          // K

    std::vector<double> X(mix.nSpecies(), 0.0);
    X[iN2] = 0.79;
    X[iO2] = 0.21;

    // Initial Density
    double Mw_mix = 0.0;
    for (int i = 0; i < mix.nSpecies(); ++i)
        Mw_mix += X[i] * mix.speciesMw(i);
    double R_univ = 8.31446;
    double rho_init = P_init / ((R_univ / Mw_mix) * T_init);

    CellState cell;
    cell.rho = rho_init;
    cell.T = T_init;
    cell.Tv = T_init;
    cell.rhoY.resize(mix.nSpecies());
    for (int i = 0; i < mix.nSpecies(); ++i)
    {
        cell.rhoY[i] = cell.rho * (X[i] * mix.speciesMw(i) / Mw_mix);
    }

    // ---- UPDATED TIME SETTINGS ----
    double t = 0.0;
    double t_end = 1.0e-3;
    double dt = 1.0e-9; // Reduced to 1 nanosecond for stability

    std::cout << "time,T,rho_N2,rho_O2,rho_NO,rho_N,rho_O" << std::endl;

    long long step = 0;
    while (t < t_end)
    {
        // Output every 10 microseconds (10,000 steps)
        if (step % 10000 == 0)
        {
            std::cout << t << "," << cell.T
                      << "," << cell.rhoY[iN2] << "," << cell.rhoY[iO2]
                      << "," << cell.rhoY[iNO] << "," << cell.rhoY[iN]
                      << "," << cell.rhoY[iO] << std::endl;
        }

        chemistry_relaxation_step(mix, cell, dt);

        t += dt;
        step++;
    }

    return 0;
}