#include <mutation++.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

using namespace Mutation;

// ------------------------------------------------------------
// Constants
// ------------------------------------------------------------
static constexpr double Ru = 8.31446261815324; // J/mol/K
static constexpr double ATM = 101325.0;        // Pa

using Mutation::KB;
using Mutation::NA;
using Mutation::PI;

// ------------------------------------------------------------
// Vibrational species data (paper Appendix A)
// ------------------------------------------------------------
struct VibSpecies
{
    const char *name;
    double M;       // kg/mol
    double theta_v; // K
};

static const VibSpecies vibList[] =
    {
        {"N2", 28.0134e-3, 3371.0},
        {"O2", 31.9988e-3, 2256.0},
        {"NO", 30.0061e-3, 2719.0}};

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------
inline double Rs(double M) { return Ru / M; }

// Eq.(5): vibrational energy per unit mass
inline double ev(double T, double theta_v, double M)
{
    const double x = theta_v / T;
    return Rs(M) * theta_v / (std::exp(x) - 1.0);
}

// d(ev)/dT for Newton inversion
inline double devdT(double T, double theta_v, double M)
{
    const double x = theta_v / T;
    const double ex = std::exp(x);
    return Rs(M) * theta_v * (ex * theta_v / (T * T)) /
           ((ex - 1.0) * (ex - 1.0));
}

// ------------------------------------------------------------
// Invert Ev -> Tv using Newton
// ------------------------------------------------------------
double invertTv(
    double Ev_target,
    double rho,
    const std::vector<double> &Y,
    const std::vector<int> &vibIdx,
    const std::vector<double> &vibM,
    const std::vector<double> &vibTheta)
{
    double Tv = 3000.0;

    for (int it = 0; it < 30; ++it)
    {
        double Ev = 0.0;
        double dEv = 0.0;

        for (size_t i = 0; i < vibIdx.size(); ++i)
        {
            const int s = vibIdx[i];
            if (Y[s] <= 0.0)
                continue;

            Ev += rho * Y[s] * ev(Tv, vibTheta[i], vibM[i]);
            dEv += rho * Y[s] * devdT(Tv, vibTheta[i], vibM[i]);
        }

        const double f = Ev - Ev_target;
        if (std::abs(f) < 1e-6 * std::max(1.0, Ev_target))
            break;

        Tv -= f / dEv;
        Tv = std::max(300.0, std::min(25000.0, Tv));
    }
    return Tv;
}

// ------------------------------------------------------------
// Get vibrational energy source from Mutation++
// ------------------------------------------------------------
double getQve(Mixture &mix)
{
    std::vector<double> src(mix.nEnergyEqns(), 0.0);
    mix.energyTransferSource(src.data());
    return src[0]; // Vib–electronic source [J/m^3/s]
}

// ------------------------------------------------------------
// MAIN
// ------------------------------------------------------------
int main()
{
    // ------------------------------------------------------------
    // Mutation++ setup (REQUIRED)
    // ------------------------------------------------------------
    MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    Mixture mix(opts);

    const int ns = mix.nSpecies();

    // ------------------------------------------------------------
    // Composition (air)
    // ------------------------------------------------------------
    std::vector<double> Y(ns, 0.0);
    const int iN2 = mix.speciesIndex("N2");
    const int iO2 = mix.speciesIndex("O2");

    if (iN2 < 0 || iO2 < 0)
        throw std::runtime_error("N2 or O2 not found");

    Y[iN2] = 0.79;
    Y[iO2] = 0.21;

    // ------------------------------------------------------------
    // Molar masses
    // ------------------------------------------------------------
    std::vector<double> Mmol(ns, 0.0);
    auto setM = [&](const char *name, double M)
    {
        int idx = mix.speciesIndex(name);
        if (idx >= 0)
            Mmol[idx] = M;
    };

    setM("N2", 28.0134e-3);
    setM("O2", 31.9988e-3);
    setM("NO", 30.0061e-3);
    setM("N", 14.0067e-3);
    setM("O", 15.9994e-3);

    // ------------------------------------------------------------
    // Mixture gas constant
    // ------------------------------------------------------------
    double Rmix = 0.0;
    for (int s = 0; s < ns; ++s)
        if (Y[s] > 0.0)
            Rmix += Y[s] * Rs(Mmol[s]);

    // ------------------------------------------------------------
    // Initial conditions
    // ------------------------------------------------------------
    double Ttr = 12000.0;
    double Tv = 2000.0;
    double p0 = ATM;

    double rho = p0 / (Rmix * Ttr); // constant-volume heat bath

    // ------------------------------------------------------------
    // Vibrational species present
    // ------------------------------------------------------------
    std::vector<int> vibIdx;
    std::vector<double> vibM, vibTheta;

    for (const auto &v : vibList)
    {
        int idx = mix.speciesIndex(v.name);
        if (idx >= 0)
        {
            vibIdx.push_back(idx);
            vibM.push_back(v.M);
            vibTheta.push_back(v.theta_v);
        }
    }

    // ------------------------------------------------------------
    // Initial energies (per unit volume)
    // ------------------------------------------------------------
    double Et = rho * (1.5 * Rmix * Ttr);
    double Ev = 0.0;

    for (size_t i = 0; i < vibIdx.size(); ++i)
        Ev += rho * Y[vibIdx[i]] * ev(Tv, vibTheta[i], vibM[i]);

    // ------------------------------------------------------------
    // Output
    // ------------------------------------------------------------
    std::ofstream out("twoT_energy_based.dat");
    out << "# t  Ttr  Tv  Et  Ev  p\n";

    // ------------------------------------------------------------
    // Time loop
    // ------------------------------------------------------------
    const double dt = 1.0e-8;
    const int nSteps = 4000;

    double Tstate[2];
    double t = 0.0;

    for (int n = 0; n < nSteps; ++n)
    {
        // Recover temperatures
        Ttr = Et / (rho * 1.5 * Rmix);
        double p = rho * Rmix * Ttr;
        Tv = invertTv(Ev, rho, Y, vibIdx, vibM, vibTheta);

        // Update Mutation++ state
        Tstate[0] = Ttr;
        Tstate[1] = Tv;
        std::vector<double> rho_i(ns);
        for (int s = 0; s < ns; ++s)
            rho_i[s] = rho * Y[s];

        mix.setState(rho_i.data(), Tstate, 1);

        // Write
        out << t << " " << Ttr << " " << Tv
            << " " << Et << " " << Ev << " " << p << "\n";

        // Source term from Mutation++
        double Qve = getQve(mix);

        // Conservative update
        Ev += Qve * dt;
        Et -= Qve * dt;

        t += dt;
    }

    out.close();
    std::cout << "Wrote: twoT_energy_based.dat\n";
    return 0;
}
