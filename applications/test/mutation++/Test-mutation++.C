#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <mutation++.h>

// Helper: Ensure density never goes negative
static inline double pos(double x) { return (x < 1.0e-30) ? 1.0e-30 : x; }

struct CellState
{
    double rho;
    double T;
    double Tv;
    std::vector<double> rhoY;
};

// --- ADAPTIVE TIME STEPPING ---
double get_adaptive_dt(const std::vector<double> &rhoY, const std::vector<double> &wdot, double current_dt)
{
    double min_timescale = 1.0e-6;

    for (size_t i = 0; i < rhoY.size(); ++i)
    {
        if (std::abs(wdot[i]) > 1.0e-5)
        {
            // 0.5% density change limit (Very conservative)
            double timescale = 0.005 * (rhoY[i] + 1.0e-10) / (std::abs(wdot[i]) + 1.0e-20);
            if (timescale < min_timescale)
                min_timescale = timescale;
        }
    }

    double target_dt = std::max(1.0e-15, std::min(min_timescale, 1.0e-6));
    // Prevent dt from growing too fast (max 2% growth)
    if (target_dt > current_dt * 1.02)
        return current_dt * 1.02;
    return target_dt;
}

void chemistry_relaxation_step(Mutation::Mixture &mix, CellState &cell, double &dt)
{
    // --- STEP 0: INITIAL ENERGIES ---
    double Tvec[] = {cell.T, cell.Tv};
    mix.setState(cell.rhoY.data(), Tvec, 2);
    double e_total_target = mix.mixtureEnergyMass(); // Conserved Quantity

    // Current Vib Energy
    std::vector<double> energies(mix.nEnergyEqns());
    mix.mixtureEnergies(energies.data());
    double Ev_density_old = energies[1] * cell.rho;

    // --- STEP 1: RATES & SOURCES ---
    std::vector<double> wdot(mix.nSpecies());
    mix.netProductionRates(wdot.data());

    std::vector<double> source_E(mix.nEnergyEqns());
    mix.energyTransferSource(source_E.data());
    double Q_vib = source_E[1];

    dt = get_adaptive_dt(cell.rhoY, wdot, dt);

    // --- STEP 2: UPDATE SPECIES ---
    std::vector<double> rhoY_new = cell.rhoY;
    double rho_sum = 0.0;
    for (int i = 0; i < mix.nSpecies(); ++i)
    {
        double change = wdot[i] * dt;
        // Hard Clamp: Max 5% change per step
        double limit = 0.05 * cell.rhoY[i];
        if (change > limit)
            change = limit;
        if (change < -limit)
            change = -limit;

        rhoY_new[i] += change;
        rhoY_new[i] = pos(rhoY_new[i]);
        rho_sum += rhoY_new[i];
    }
    // Mass correction
    double factor = cell.rho / rho_sum;
    for (int i = 0; i < mix.nSpecies(); ++i)
        rhoY_new[i] *= factor;
    cell.rhoY = rhoY_new;

    // --- STEP 3: UPDATE VIB ENERGY ---
    // Explicit update
    double Ev_density_new = Ev_density_old + Q_vib * dt;
    if (Ev_density_new < 1e-5)
        Ev_density_new = 1e-5;

    // Find Tv matching new Energy
    double Tv_new = cell.Tv;
    for (int k = 0; k < 10; ++k)
    {
        double T_dummy[] = {cell.T, Tv_new};
        mix.setState(cell.rhoY.data(), T_dummy, 2);
        mix.mixtureEnergies(energies.data());
        double Ev_curr = energies[1] * cell.rho;

        double diff = Ev_density_new - Ev_curr;
        if (std::abs(diff) < 1.0)
            break;

        double Cv_approx = Ev_curr / Tv_new;
        if (Cv_approx < 1e-3)
            Cv_approx = 1e-3;

        double dTv = diff / Cv_approx;
        // Clamp dTv
        if (dTv > 100.0)
            dTv = 100.0;
        if (dTv < -100.0)
            dTv = -100.0;

        Tv_new += dTv;
    }
    cell.Tv = Tv_new;

    // --- STEP 4: UPDATE TRANSLATIONAL T ---
    // Solve Energy Equation: E_total(T) = e_total_target
    double T_new = cell.T;
    for (int iter = 0; iter < 20; ++iter)
    {
        double T_dual[] = {T_new, cell.Tv};
        mix.setState(cell.rhoY.data(), T_dual, 2);

        double e_curr = mix.mixtureEnergyMass();
        double Cv = mix.mixtureFrozenCvMass();

        double error = e_total_target - e_curr;
        if (std::abs(error) < 1.0)
            break;

        double dT = error / Cv;

        // --- STIFF LIMITER ---
        // If the physics wants to change T by > 50K, we stop it.
        // This forces the loop to run stable.
        if (dT > 50.0)
            dT = 50.0;
        if (dT < -50.0)
            dT = -50.0;

        T_new += dT;
    }
    cell.T = T_new;
}

int main()
{
    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    // Initial Conditions
    double P_init = 0.063 * 101325.0;
    double T_init = 10000.0;
    double Tv_init = 1000.0; // Non-Equilibrium Start

    int iN2 = mix.speciesIndex("N2");
    int iO2 = mix.speciesIndex("O2");
    int iNO = mix.speciesIndex("NO");
    int iN = mix.speciesIndex("N");
    int iO = mix.speciesIndex("O");

    std::vector<double> X(mix.nSpecies(), 0.0);
    X[iN2] = 0.79;
    X[iO2] = 0.21;

    double Mw_mix = 0.0;
    for (int i = 0; i < mix.nSpecies(); ++i)
        Mw_mix += X[i] * mix.speciesMw(i);
    double rho_init = P_init / ((8.31446 / Mw_mix) * T_init);

    CellState cell;
    cell.rho = rho_init;
    cell.T = T_init;
    cell.Tv = Tv_init;
    cell.rhoY.resize(mix.nSpecies());
    for (int i = 0; i < mix.nSpecies(); ++i)
        cell.rhoY[i] = cell.rho * (X[i] * mix.speciesMw(i) / Mw_mix);

    std::ofstream outFile("air_2T_data.csv");
    outFile << "time,T,Tv,rho_N2,rho_O2,rho_NO,rho_N,rho_O" << std::endl;

    double t = 0.0;
    double t_end = 1.0e-3;
    double dt = 1.0e-12; // Start extremely small

    std::cout << "Starting 2T Sim. T=" << cell.T << " Tv=" << cell.Tv << std::endl;
    int step = 0;
    int print_freq = 100;

    while (t < t_end)
    {
        if (step % print_freq == 0)
        {
            outFile << t << "," << cell.T << "," << cell.Tv
                    << "," << cell.rhoY[iN2] << "," << cell.rhoY[iO2]
                    << "," << cell.rhoY[iNO] << "," << cell.rhoY[iN]
                    << "," << cell.rhoY[iO] << std::endl;

            // Console Debug
            if (step % 1000 == 0)
                std::cout << "t=" << t << " T=" << cell.T << " Tv=" << cell.Tv << std::endl;

            if (dt > 1e-10)
                print_freq = 1000;
        }

        chemistry_relaxation_step(mix, cell, dt);
        t += dt;
        step++;
    }

    outFile.close();
    std::cout << "Done." << std::endl;
    return 0;
}