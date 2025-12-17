#include "mutationMixture.H"
#include <cmath>
#include <iostream>

using namespace Mutation;

// ------------------------------------------------------------
// Constants
// ------------------------------------------------------------
static constexpr double Ru = 8.31446261815324; // J/mol/K

// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
mutationMixture::mutationMixture(const std::string &mechanism)
    : mix_(
          [&]()
          {
              MixtureOptions opts(mechanism);
              opts.setStateModel("ChemNonEqTTv");
              opts.setThermodynamicDatabase("RRHO");
              return Mixture(opts);
          }())
{
    const int ns = mix_.nSpecies();

    rho_i_.resize(ns);
    Tstate_.resize(2);

    // --------------------------------------------------------
    // Paper vibrational data (Appendix A)
    // --------------------------------------------------------
    struct VibInit
    {
        const char *name;
        double M;
        double theta;
    };

    static const VibInit vibInit[] =
        {
            {"N2", 28.0134e-3, 3371.0},
            {"O2", 31.9988e-3, 2256.0},
            {"NO", 30.0061e-3, 2719.0}};

    for (const auto &v : vibInit)
    {
        int idx = mix_.speciesIndex(v.name);
        if (idx >= 0)
        {
            vibData_.push_back({idx, v.M, v.theta});
        }
    }

    std::cout << "mutationMixture initialized with "
              << vibData_.size()
              << " vibrational species\n";
}

// ------------------------------------------------------------
// ONE TIME STEP (paper heat-bath model)
// ------------------------------------------------------------
void mutationMixture::step(
    double dt,
    double rho,
    const std::vector<double> &Y,
    double &Et,
    double &Ev,
    double &Ttr,
    double &Tv)
{
    const int ns = mix_.nSpecies();

    // --------------------------------------------------------
    // Recover translational temperature (PAPER MODEL)
    // e_tr = rho * (1.5 * Rmix * Ttr)
    // --------------------------------------------------------
    double Rmix = 0.0;
    for (int s = 0; s < ns; ++s)
    {
        if (Y[s] > 0.0)
            Rmix += Y[s] * Ru / mix_.speciesMw(s);
    }

    Ttr = Et / (rho * 1.5 * Rmix);

    // --------------------------------------------------------
    // Recover vibrational temperature
    // --------------------------------------------------------
    Tv = invertTv(Ev, rho, Y);

    // --------------------------------------------------------
    // Set Mutation++ state
    // --------------------------------------------------------
    for (int s = 0; s < ns; ++s)
        rho_i_[s] = rho * Y[s];

    Tstate_[0] = Ttr;
    Tstate_[1] = Tv;

    // State model = 1 → density + temperatures
    mix_.setState(rho_i_.data(), Tstate_.data(), 1);

    // --------------------------------------------------------
    // Vibrational energy source term Qve
    // --------------------------------------------------------
    std::vector<double> src(mix_.nEnergyEqns(), 0.0);
    mix_.energyTransferSource(src.data());

    const double Qve = src[0]; // J/m^3/s

    // --------------------------------------------------------
    // Conservative update
    // --------------------------------------------------------
    Ev += Qve * dt;
    Et -= Qve * dt;
}

// ------------------------------------------------------------
// Invert Ev -> Tv (Newton method, paper Eq. 5)
// ------------------------------------------------------------
double mutationMixture::invertTv(
    double Ev_target,
    double rho,
    const std::vector<double> &Y) const
{
    double Tv = 3000.0;

    for (int it = 0; it < 30; ++it)
    {
        double Ev = 0.0;
        double dEv = 0.0;

        for (const auto &v : vibData_)
        {
            const int s = v.speciesIndex;
            if (Y[s] <= 0.0)
                continue;

            const double Rs = Ru / v.M;
            const double x = v.theta_v / Tv;
            const double ex = std::exp(x);

            const double ev =
                Rs * v.theta_v / (ex - 1.0);

            const double devdT =
                Rs * v.theta_v *
                (ex * v.theta_v / (Tv * Tv)) /
                ((ex - 1.0) * (ex - 1.0));

            Ev += rho * Y[s] * ev;
            dEv += rho * Y[s] * devdT;
        }

        const double f = Ev - Ev_target;
        if (std::abs(f) < 1e-6 * std::max(1.0, Ev_target))
            break;

        Tv -= f / dEv;

        // C++11-safe clamp
        if (Tv < 300.0)
            Tv = 300.0;
        if (Tv > 25000.0)
            Tv = 25000.0;
    }

    return Tv;
}
