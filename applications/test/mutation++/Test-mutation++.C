#include "mutationMixture.H"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// ---------- Constants ----------
static constexpr double Ru = 8.31446261815324; // J/mol/K
static constexpr double ATM = 101325.0;        // Pa

int main()
{
    try
    {
        // ------------------------------------------------------------
        // Create wrapper (no direct Mutation::Mixture in this file)
        // ------------------------------------------------------------
        mutationMixture mix("air_5");

        const int ns = mix.nSpecies();
        if (ns <= 0)
            throw std::runtime_error("Mixture has no species.");

        // ------------------------------------------------------------
        // Build species name list (for file header)
        // ------------------------------------------------------------
        std::vector<std::string> spNames(ns);
        for (int s = 0; s < ns; ++s)
            spNames[s] = mix.speciesName(s);

        // ------------------------------------------------------------
        // Composition (air)
        // ------------------------------------------------------------
        std::vector<double> Y(ns, 0.0);

        const int iN2 = mix.speciesIndex("N2");
        const int iO2 = mix.speciesIndex("O2");

        if (iN2 < 0 || iO2 < 0)
            throw std::runtime_error("N2 and/or O2 not found in air_5.");

        Y[iN2] = 0.79;
        Y[iO2] = 0.21;

        // ------------------------------------------------------------
        // Compute mixture gas constant Rmix = sum(Ys * Ru/Ms)
        // ------------------------------------------------------------
        double Rmix = 0.0;
        for (int s = 0; s < ns; ++s)
        {
            if (Y[s] <= 0.0)
                continue;
            const double Mw = mix.speciesMw(s); // kg/mol
            if (Mw <= 0.0)
                continue;
            Rmix += Y[s] * Ru / Mw; // J/kg/K
        }

        if (Rmix <= 0.0)
            throw std::runtime_error("Invalid Rmix (check Y and molecular weights).");

        // ------------------------------------------------------------
        // Initial conditions (paper heat-bath style)
        // ------------------------------------------------------------
        double Ttr = 12000.0;
        double Tv = 2000.0;

        const double p0 = ATM;
        const double rho = p0 / (Rmix * Ttr); // constant-volume test

        // Initial energies (translation-only model used by your wrapper)
        double Et = rho * (1.5 * Rmix * Ttr);
        double Ev = 0.0; // OK: wrapper will move energy into Ev via Qve

        // ------------------------------------------------------------
        // Time loop
        // ------------------------------------------------------------
        const double dt = 1.0e-8;
        const int nSteps = 4000;
        double t = 0.0;

        std::ofstream out("twoT_energy_based.dat");
        if (!out)
            throw std::runtime_error("Cannot open output file twoT_energy_based.dat");

        // Header
        out << "# t[s] Ttr[K] Tv[K] Et[J/m3] Ev[J/m3] p[Pa]";
        for (int s = 0; s < ns; ++s)
            out << " rho_" << spNames[s];
        out << "\n";

        for (int n = 0; n < nSteps; ++n)
        {
            // One relaxation step (updates Et/Ev and returns updated Ttr/Tv)
            mix.step(dt, rho, Y, Et, Ev, Ttr, Tv);

            // Pressure (ideal gas, based on translational temperature)
            const double p = rho * Rmix * Ttr;

            // Write line
            out << t << " " << Ttr << " " << Tv << " " << Et << " " << Ev << " " << p;

            // Species densities (constant here unless you later evolve Y)
            for (int s = 0; s < ns; ++s)
                out << " " << (rho * Y[s]);

            out << "\n";

            t += dt;
        }

        out.close();
        std::cout << "Wrote: twoT_energy_based.dat\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
