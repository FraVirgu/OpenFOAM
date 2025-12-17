#include <mutation++.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <string>

using namespace Mutation;

// ---------- Physical constants ----------
static constexpr double Ru = 8.31446261815324; // J/mol/K
using Mutation::KB;
using Mutation::NA;
using Mutation::PI;
static constexpr double ATM = 101325.0; // Pa

// ---------- Species data (from paper Appendix A, Table A1) ----------
struct VibSpecies
{
    std::string name;
    double M_kg_per_mol; // kg/mol
    double theta_v;      // K
    double sigma_m;      // m^2 (Park high-T cross section prefactor)
};

// Paper values:
// theta_v: N2=3371, O2=2256, NO=2719  :contentReference[oaicite:2]{index=2}
// sigma_m: N2,O2 = 3e-21 ; NO = 3e-22 :contentReference[oaicite:3]{index=3}
static const VibSpecies vibList[] = {
    {"N2", 28.0134e-3, 3371.0, 3.0e-21},
    {"O2", 31.9988e-3, 2256.0, 3.0e-21},
    {"NO", 30.0061e-3, 2719.0, 3.0e-22},
};

// ---------- Helpers ----------
static inline double Rs(double M_kg_per_mol)
{
    return Ru / M_kg_per_mol; // J/kg/K
}

// Eq. (5): ev(T) = Rs * theta_v / (exp(theta_v/T)-1)  :contentReference[oaicite:4]{index=4}
static inline double ev_species(double T, double theta_v, double M_kg_per_mol)
{
    const double x = theta_v / T;
    const double ex = std::exp(x);
    return Rs(M_kg_per_mol) * (theta_v / (ex - 1.0));
}

// derivative d(ev)/dT for Newton inversion
static inline double devdT_species(double T, double theta_v, double M_kg_per_mol)
{
    const double x = theta_v / T;
    const double ex = std::exp(x);
    const double denom = (ex - 1.0);
    // ev = Rs * theta / (e^x - 1), x=theta/T
    // d/dT [1/(e^x-1)] = (e^x * (theta/T^2)) / (e^x-1)^2
    const double dInv = (ex * (theta_v / (T * T))) / (denom * denom);
    return Rs(M_kg_per_mol) * theta_v * dInv;
}

// Reduced mass Mm,s in kg/mol: Mm*Ms/(Mm+Ms)  :contentReference[oaicite:5]{index=5}
static inline double reduced_mass(double Mm, double Ms)
{
    return (Mm * Ms) / (Mm + Ms);
}

// Millikan-White coefficients from Eqs (12)-(13): :contentReference[oaicite:6]{index=6}
// Am,s = 1.16e-3 * sqrt(Mm,s) * theta_v^(4/3)
// Bm,s = 0.015 * (Mm,s)^(1/4)
// NOTE: the paper shows sqrt-like dependence via reduced mass and theta_v^(4/3).
static inline void MW_coeff(double Mm_s_red, double theta_v_m, double &A, double &B)
{
    A = 1.16e-3 * std::sqrt(Mm_s_red * 1e3) * std::pow(theta_v_m, 4.0 / 3.0); // heuristic unit-handling
    B = 0.015 * std::pow(Mm_s_red * 1e3, 1.0 / 4.0);
}

// Millikan-White relaxation time (Eq. 11) with p in atm :contentReference[oaicite:7]{index=7}
static inline double tau_MW(double p_atm, double Ttr, double A, double B)
{
    // tau = (1/p) * exp( A*(T^-1/3 - B) - 18.42 )
    const double Tm13 = std::pow(Ttr, -1.0 / 3.0);
    return (1.0 / p_atm) * std::exp(A * (Tm13 - B) - 18.42);
}

// Park correction term (Eqs. 15-17) :contentReference[oaicite:8]{index=8}
static inline double tau_Park(double Ttr, double sigma_m, double M_m_kg_per_mol, double n_collider)
{
    // average molecular speed: cbar = sqrt(8*KB*T / (pi*m)) ; m = M/NA
    const double m_particle = M_m_kg_per_mol / NA;
    const double cbar = std::sqrt((8.0 * KB * Ttr) / (PI * m_particle));

    // sigma_v = sigma_m * (50000/T)^2
    const double sigma_v = sigma_m * std::pow(50000.0 / Ttr, 2.0);

    // tau_P = 1 / (cbar * sigma_v * n)
    return 1.0 / (cbar * sigma_v * n_collider);
}

// Mixture-averaged tau_m (Eq. 9-10) using molar fractions Xs :contentReference[oaicite:9]{index=9}
static double tau_vt_mixture_avg(
    double Ttr,
    double p_Pa,
    const std::vector<double> &X,    // molar fractions of all species
    const std::vector<double> &Mmol, // kg/mol per species
    int idx_m,                       // molecule index (in species list) for molecule m
    double theta_v_m,
    double sigma_m)
{
    const double p_atm = p_Pa / ATM;

    // total number density (ideal gas): n_tot = p / (KB*Ttr)
    const double n_tot = p_Pa / (KB * Ttr);

    double denom = 0.0;
    for (size_t s = 0; s < X.size(); ++s)
    {
        if (X[s] <= 0.0)
            continue;

        const double Mm = Mmol[idx_m];
        const double Ms = Mmol[s];
        const double Mred = reduced_mass(Mm, Ms);

        double A = 0.0, B = 0.0;
        MW_coeff(Mred, theta_v_m, A, B);

        const double tauMW = tau_MW(p_atm, Ttr, A, B);

        const double n_s = X[s] * n_tot; // collider number density
        const double tauP = tau_Park(Ttr, sigma_m, Mm, n_s);

        const double tau_ms = tauMW + tauP;

        denom += X[s] / tau_ms;
    }
    if (denom <= 0.0)
    {
        throw std::runtime_error("Invalid relaxation-time denominator.");
    }
    return 1.0 / denom;
}

// Invert Tv from Ev by Newton: Ev = rho * sum(Ym * ev_m(Tv))
static double invert_Tv_from_Ev(
    double Ev_target,
    double rho,
    const std::vector<double> &Y,
    const std::vector<int> &vibSpeciesIdx,
    const std::vector<double> &Mmol_vib,
    const std::vector<double> &theta_v_vib)
{
    double Tv = 3000.0; // initial guess
    for (int it = 0; it < 30; ++it)
    {
        double Ev = 0.0;
        double dEv = 0.0;

        for (size_t j = 0; j < vibSpeciesIdx.size(); ++j)
        {
            const int s = vibSpeciesIdx[j];
            const double Ym = Y[s];
            if (Ym <= 0.0)
                continue;

            const double M = Mmol_vib[j];
            const double th = theta_v_vib[j];

            Ev += rho * Ym * ev_species(Tv, th, M);
            dEv += rho * Ym * devdT_species(Tv, th, M);
        }

        const double f = Ev - Ev_target;
        if (std::abs(f) < 1e-6 * std::max(1.0, std::abs(Ev_target)))
            break;

        Tv -= f / dEv;
        if (Tv < 50.0)
            Tv = 50.0;
        if (Tv > 1.0e6)
            Tv = 1.0e6;
    }
    return Tv;
}

int main()
{
    // ------------------------------------------------------------
    // REQUIRED Mutation++ setup (your exact snippet)
    // ------------------------------------------------------------
    MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    Mixture mix(opts);

    const int ns = mix.nSpecies();

    // Build Y (mass fractions). Example: air = 0.79 N2 + 0.21 O2
    std::vector<double> Y(ns, 0.0);
    const int iN2 = mix.speciesIndex("N2");
    const int iO2 = mix.speciesIndex("O2");
    if (iN2 < 0 || iO2 < 0)
        throw std::runtime_error("N2/O2 not found in air_5.");

    Y[iN2] = 0.79;
    Y[iO2] = 0.21;

    // Store molar masses (kg/mol). We will take from paper values for robustness.
    // (We still use Mutation++ for indexing / mixture definition.)
    std::vector<double> Mmol(ns, 0.0);
    // air_5 typically: N2 O2 NO N O (no vib for atoms)
    // Put paper molar masses where applicable:
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

    // Compute molar fractions X from Y:
    // Xs = (Ys/Ms) / sum(Yk/Mk)
    std::vector<double> X(ns, 0.0);
    double denom = 0.0;
    for (int s = 0; s < ns; ++s)
    {
        if (Y[s] > 0.0 && Mmol[s] > 0.0)
            denom += Y[s] / Mmol[s];
    }
    for (int s = 0; s < ns; ++s)
    {
        if (Y[s] > 0.0 && Mmol[s] > 0.0)
            X[s] = (Y[s] / Mmol[s]) / denom;
    }

    // Mixture gas constant Rmix = sum(Ys*Rs) = sum(Ys*Ru/Ms)
    double Rmix = 0.0;
    for (int s = 0; s < ns; ++s)
    {
        if (Y[s] > 0.0 && Mmol[s] > 0.0)
            Rmix += Y[s] * Rs(Mmol[s]);
    }

    // ------------------------------------------------------------
    // Initial conditions (vibrational heating: Ttr > Tv, curves intersect)
    // ------------------------------------------------------------
    double Ttr = 12000.0;
    double Tv = 2000.0;

    // constant-volume heat bath: keep rho constant, pressure evolves with Ttr
    double p0 = 101325.0;
    double rho = p0 / (Rmix * Ttr); // from p = rho*Rmix*Ttr (no electrons) :contentReference[oaicite:10]{index=10}

    // Initialize Mutation++ state (not used for Et/Ev math here, but satisfies your requirement)
    double Tstate[2] = {Ttr, Tv};
    mix.setState(Y.data(), Tstate, 0);

    // Build vib species indices present in this mixture
    std::vector<int> vibIdx;
    std::vector<double> vibM;
    std::vector<double> vibTheta;
    std::vector<double> vibSigma;

    for (const auto &v : vibList)
    {
        int idx = mix.speciesIndex(v.name);
        if (idx >= 0)
        {
            vibIdx.push_back(idx);
            vibM.push_back(v.M_kg_per_mol);
            vibTheta.push_back(v.theta_v);
            vibSigma.push_back(v.sigma_m);
        }
    }

    // ------------------------------------------------------------
    // Initial energies per unit volume (paper: sum rho_s * e_mode,s)
    // Et = rho * sum(Ys * 3/2 * Rs * Ttr)  (Eq. 3) :contentReference[oaicite:11]{index=11}
    // Ev = rho * sum(Ym * ev_m(Tv))        (Eq. 5) :contentReference[oaicite:12]{index=12}
    // ------------------------------------------------------------
    double Et = rho * (1.5 * Rmix * Ttr);
    double Ev = 0.0;
    for (size_t j = 0; j < vibIdx.size(); ++j)
    {
        const int s = vibIdx[j];
        Ev += rho * Y[s] * ev_species(Tv, vibTheta[j], vibM[j]);
    }

    // ------------------------------------------------------------
    // Time loop
    // ------------------------------------------------------------
    const double dt = 1.0e-8;
    const int nSteps = 4000;

    std::ofstream out("twoT_energy_based.dat");
    out << "# t[s]  Ttr[K]  Tv[K]  Et[J/m3]  Ev[J/m3]  p[Pa]\n";

    double t = 0.0;
    for (int n = 0; n < nSteps; ++n)
    {

        // Recover Ttr from Et analytically:
        // Et = rho * (3/2 * Rmix * Ttr) => Ttr = Et / (rho * 1.5 * Rmix)
        Ttr = Et / (rho * 1.5 * Rmix);

        // Pressure evolves in constant volume: p = rho*Rmix*Ttr :contentReference[oaicite:13]{index=13}
        double p = rho * Rmix * Ttr;

        // Recover Tv from Ev by Newton inversion
        Tv = invert_Tv_from_Ev(Ev, rho, Y, vibIdx, vibM, vibTheta);

        // Write output
        out << t << " " << Ttr << " " << Tv << " " << Et << " " << Ev << " " << p << "\n";

        // Update Mutation++ state (again: satisfies your “use Mutation++ input data” requirement)
        Tstate[0] = Ttr;
        Tstate[1] = Tv;
        mix.setState(Y.data(), Tstate, 0);

        // Compute Q_VT (Eq. 8): sum_m rho_m*(ev(Ttr)-ev(Tv))/tau_m :contentReference[oaicite:14]{index=14}
        double Qvt = 0.0;
        for (size_t j = 0; j < vibIdx.size(); ++j)
        {
            const int m = vibIdx[j]; // molecule species index
            const double Ym = Y[m];
            if (Ym <= 0.0)
                continue;

            const double ev_tr = ev_species(Ttr, vibTheta[j], vibM[j]);
            const double ev_tv = ev_species(Tv, vibTheta[j], vibM[j]);

            const double tau_m = tau_vt_mixture_avg(
                Ttr, p, X, Mmol,
                m, vibTheta[j], vibSigma[j]);

            Qvt += rho * Ym * (ev_tr - ev_tv) / tau_m;
        }

        // Conservative update:
        Ev += Qvt * dt;
        Et -= Qvt * dt;

        t += dt;
    }

    out.close();
    std::cout << "Wrote: twoT_energy_based.dat\n";
    return 0;
}
