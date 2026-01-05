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

    // Compute Rmix (still needed for any consistency checks if you want)
    double Rmix = 0.0;
    for (int s = 0; s < ns; ++s)
    {
        if (Y[s] > 0.0)
            Rmix += Y[s] * Ru / mix_.speciesMw(s);
    }

    // Use the temperatures provided by the caller (OpenFOAM)
    if (!(std::isfinite(Ttr) && Ttr > 0.0))
        Ttr = 300.0;
    if (!(std::isfinite(Tv) && Tv > 0.0))
        Tv = Ttr;

    // --------------------------------------------------------
    // Set Mutation++ state
    // --------------------------------------------------------
    for (int s = 0; s < ns; ++s)
    {
        rho_i_[s] = std::max(rho * Y[s], 1e-12);
    }

    Tstate_[0] = Ttr;
    Tstate_[1] = Tv;

    // State model = 1 → density + temperatures
    mix_.setState(rho_i_.data(), Tstate_.data(), 1);
    const double pIdeal = rho * Rmix * Ttr;
    std::cerr << "ideal p=" << pIdeal << "\n";

    std::cerr << "ideal p=" << pIdeal
              << " Ttr=" << Ttr << " Tv=" << Tv << "\n";

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

double mutationMixture::EvFromTv(double Tv, double rho, const std::vector<double> &Y) const
{
    double Ev = 0.0;
    for (const auto &v : vibData_)
    {
        const int s = v.speciesIndex;
        if (Y[s] <= 0.0)
            continue;

        const double Rs = Ru / v.M;
        const double x = v.theta_v / Tv;
        const double ex = std::exp(x);

        const double ev = Rs * v.theta_v / (ex - 1.0); // J/kg
        Ev += rho * Y[s] * ev;                         // J/m3
    }
    return Ev;
}

// ------------------------------------------------------------
// Invert Ev -> Tv (Newton method, paper Eq. 5)
// ------------------------------------------------------------
double mutationMixture::invertTv(
    double Ev_target,
    double rho,
    const std::vector<double> &Y,
    double Tv_init) const
{
    // Quick exits
    if (!(std::isfinite(Ev_target)) || Ev_target <= 0.0)
        return 300.0;

    double Tv = Tv_init;
    if (!(std::isfinite(Tv) && Tv > 0.0))
        Tv = 3000.0;

    // clamp initial
    if (Tv < 300.0)
        Tv = 300.0;
    if (Tv > 25000.0)
        Tv = 25000.0;

    // Newton with damping
    for (int it = 0; it < 30; ++it)
    {
        double Ev = 0.0;
        double dEv = 0.0;

        for (const auto &v : vibData_)
        {
            const int s = v.speciesIndex;
            const double Ys = Y[s];
            if (Ys <= 0.0)
                continue;

            const double Rs = Ru / v.M;

            // x = theta/T
            double x = v.theta_v / Tv;
            // avoid exp overflow in double
            if (x > 700.0)
                x = 700.0;

            // expm1(x) = exp(x)-1, stable for small x
            const double denom = std::expm1(x);
            if (!(std::isfinite(denom)) || denom <= 0.0)
                continue;

            const double ev = Rs * v.theta_v / denom; // J/kg

            const double ex = std::exp(x);
            const double devdT =
                Rs * v.theta_v *
                (ex * v.theta_v / (Tv * Tv)) /
                (denom * denom);

            Ev += rho * Ys * ev;     // J/m^3
            dEv += rho * Ys * devdT; // J/m^3/K
        }

        if (!(std::isfinite(Ev) && std::isfinite(dEv)) || dEv <= 0.0)
            break;

        const double f = Ev - Ev_target;

        if (std::abs(f) < 1e-8 * std::max(1.0, Ev_target))
            break;

        double step = f / dEv;

        // damping: don’t jump too far
        const double maxStep = 0.5 * Tv;
        if (step > maxStep)
            step = maxStep;
        if (step < -maxStep)
            step = -maxStep;

        Tv -= step;

        if (Tv < 300.0)
            Tv = 300.0;
        if (Tv > 25000.0)
            Tv = 25000.0;
    }

    return Tv;
}

// -----------------------------------------------------------------------------
// YOU MUST IMPLEMENT THIS using Mutation++ calls you have available.
// It must return translational/heavy-mode energy density Et(Ttr,Tv,rho,Y) in J/m^3.
// -----------------------------------------------------------------------------
double mutationMixture::EtFromState_(
    double Ttr,
    double Tv,
    double rho,
    const std::vector<double> &Y)
{
    // Example skeleton: setState then query a property.
    // You MUST replace the "TODO" part with your actual mix_ getter.

    for (int s = 0; s < mix_.nSpecies(); ++s)
        rho_i_[s] = std::max(rho * Y[s], 1e-12);

    Tstate_[0] = Ttr;
    Tstate_[1] = Tv;

    mix_.setState(rho_i_.data(), Tstate_.data(), 1);

    // Total mixture energy per mass from Mutation++ (J/kg)
    const double eTot_mass = mix_.mixtureEnergyMass();

    // Convert to per volume (J/m^3)
    const double ETot_vol = rho * eTot_mass;

    // Subtract your vib energy model (J/m^3) to get Et reservoir consistent with your split
    const double Ev_vol = EvFromTv(Tv, rho, Y);

    return std::max(ETot_vol - Ev_vol, 0.0);
}

double mutationMixture::invertTtr(
    double Et_target,
    double rho,
    const std::vector<double> &Y,
    double Tv,
    double Ttr_init)
{
    if (!(std::isfinite(Et_target)) || Et_target <= 0.0)
        return 300.0;

    double Tlo = 50.0;
    double Thi = 25000.0;

    double T = Ttr_init;
    if (!(std::isfinite(T) && T > 0.0))
        T = 3000.0;

    if (T < Tlo)
        T = Tlo;
    if (T > Thi)
        T = Thi;

    auto F = [&](double Ttr) -> double
    {
        const double Et_model = EtFromState_(Ttr, Tv, rho, Y);
        return Et_model - Et_target;
    };

    // Build a bracket [a,b] if possible
    double a = std::max(Tlo, 0.5 * T);
    double b = std::min(Thi, 2.0 * T);

    double fa = 0.0, fb = 0.0;
    bool bracketed = false;

    for (int k = 0; k < 20; ++k)
    {
        fa = F(a);
        fb = F(b);
        if (std::isfinite(fa) && std::isfinite(fb) && (fa * fb <= 0.0))
        {
            bracketed = true;
            break;
        }
        a = std::max(Tlo, a * 0.7);
        b = std::min(Thi, b * 1.3);
        if (a <= Tlo && b >= Thi)
            break;
    }

    // Safeguarded Newton
    double x = T;
    for (int it = 0; it < 30; ++it)
    {
        double fx = F(x);
        if (!std::isfinite(fx))
            break;

        if (std::abs(fx) < 1e-8 * std::max(1.0, Et_target))
            return x;

        // numerical derivative (finite difference)
        const double dx = std::max(1e-3, 1e-4 * x);
        double fpx = (F(std::min(Thi, x + dx)) - F(std::max(Tlo, x - dx))) / (2.0 * dx);

        double x_new = x;
        bool usedNewton = false;

        if (std::isfinite(fpx) && std::abs(fpx) > 0.0)
        {
            x_new = x - fx / fpx;
            usedNewton = true;
        }

        // If Newton step is bad, or no bracket: use bisection if bracketed
        if (!usedNewton || !std::isfinite(x_new) || x_new < Tlo || x_new > Thi)
        {
            if (bracketed)
                x_new = 0.5 * (a + b);
            else
                x_new = std::min(Thi, std::max(Tlo, x * (fx > 0.0 ? 0.9 : 1.1)));
        }

        // If we have a bracket, maintain it
        if (bracketed)
        {
            double fnew = F(x_new);
            if (std::isfinite(fnew))
            {
                if (fa * fnew <= 0.0)
                {
                    b = x_new;
                    fb = fnew;
                }
                else
                {
                    a = x_new;
                    fa = fnew;
                }
            }
        }

        x = x_new;
    }

    // fallback
    return x;
}
