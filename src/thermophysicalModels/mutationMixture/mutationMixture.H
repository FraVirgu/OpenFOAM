#ifndef mutationMixture_H
#define mutationMixture_H

#include <mutation++.h>
#include <vector>
#include <string>

// ------------------------------------------------------------
// Paper vibrational data container
// ------------------------------------------------------------
struct VibData
{
    int speciesIndex; // index in Mutation++ mixture
    double M;         // kg/mol
    double theta_v;   // K
};

// ------------------------------------------------------------
// mutationMixture: single-step 2T energy relaxation operator
// ------------------------------------------------------------
class mutationMixture
{
public:
    // Constructor
    explicit mutationMixture(const std::string &mechanism);

    // Perform ONE time step
    void step(
        double dt,
        double rho,
        const std::vector<double> &Y,
        double &Et,
        double &Ev,
        double &Ttr,
        double &Tv);

    // ---- pass-through accessors ----
    int nSpecies() const
    {
        return mix_.nSpecies();
    }

    int speciesIndex(const std::string &name) const
    {
        return mix_.speciesIndex(name);
    }

    std::string speciesName(int i) const
    {
        return mix_.speciesName(i);
    }

    double speciesMw(int i) const
    {
        return mix_.speciesMw(i);
    }

private:
    // Mutation++ mixture
    Mutation::Mixture mix_;

    // Working buffers
    std::vector<double> rho_i_;
    std::vector<double> Tstate_;

    // Paper vibrational species data
    std::vector<VibData> vibData_;

    // Invert Ev -> Tv (Newton method, paper definition)
    double invertTv(
        double Ev_target,
        double rho,
        const std::vector<double> &Y) const;
};

#endif
