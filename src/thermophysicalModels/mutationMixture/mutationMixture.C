#include "mutationMixture.H"

// ------------------------------------------------------------------------
// GUARD: Undefine OpenFOAM macros that conflict with Mutation++
// ------------------------------------------------------------------------
#ifdef Log
#undef Log
#endif
#ifdef log
#undef log
#endif
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#include "mutation++.h"
#include "IOobject.H"

namespace Foam
{
    defineTypeNameAndDebug(mutationMixture, 0);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    mutationMixture::mutationMixture(
        const dictionary &thermoDict,
        const fvMesh &mesh,
        const word &phaseName)
        : mesh_(mesh),
          species_(thermoDict.subDict("specie").lookup("species")),
          Y_(species_.size()),
          mixture_(nullptr),
          mechanismName_(thermoDict.lookup("mechanism"))
    {
        Info << "Initializing Mutation++ mechanism: " << mechanismName_ << endl;

        // 1. Initialize Mutation++
        Mutation::MixtureOptions opts(mechanismName_);
        opts.setStateModel("ChemNonEqTTv");
        mixture_ = new Mutation::Mixture(opts);

        // 2. Validate Species List
        if (species_.size() != static_cast<label>(mixture_->nSpecies()))
        {
            FatalErrorInFunction
                << "Species mismatch! OpenFOAM: " << species_.size()
                << ", Mutation++: " << mixture_->nSpecies()
                << exit(FatalError);
        }

        // 3. Allocate Working Vectors
        const int nSp = static_cast<int>(mixture_->nSpecies());
        Y_work_.resize(nSp);
        rho_work_.resize(nSp);
        Mw_work_.resize(nSp);
        E_work_.resize(nSp);
        H_RT_work_.resize(nSp);
        T_work_.resize(2); // [0]=T_tr, [1]=T_ve

        // 4. Pre-fill Molecular Weights
        for (int i = 0; i < nSp; ++i)
        {
            Mw_work_[i] = mixture_->speciesMw(i);
        }

        // 5. Create/Load OpenFOAM Fields
        // Uses .name() for OpenFOAM-13 compatibility
        const word timeName = mesh_.time().name();

        forAll(species_, i)
        {
            const word name = species_[i];
            Y_.set(
                i,
                new volScalarField(
                    IOobject(
                        name,
                        timeName,
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE),
                    mesh_));
        }
    }

    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    mutationMixture::~mutationMixture()
    {
        if (mixture_)
            delete mixture_;
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    void mutationMixture::setState(
        const scalar p,
        const scalar T,
        const scalar Tve,
        const scalarField &Y) const
    {
        const int nSp = static_cast<int>(species_.size());

        // 1. Fill Mass Fractions and Calc Mean MW
        double invMeanMw = 0.0;
        for (int i = 0; i < nSp; ++i)
        {
            Y_work_[i] = static_cast<double>(Y[i]);
            invMeanMw += Y_work_[i] / Mw_work_[i];
        }

        // 2. Ideal Gas Law for Mixture Density
        // rho = p / (R_specific_mix * T)
        const double Ru = Mutation::RU;
        const double rho = static_cast<double>(p) / (Ru * invMeanMw * static_cast<double>(T));

        // 3. Calculate Species Partial Densities
        for (int i = 0; i < nSp; ++i)
        {
            rho_work_[i] = rho * Y_work_[i];
        }

        // 4. Set Temperatures
        T_work_[0] = static_cast<double>(T);   // T_tr
        T_work_[1] = static_cast<double>(Tve); // T_ve

        // 5. Update Mutation++ State
        // stateModel 1 corresponds to Density + Temperature input
        mixture_->setState(rho_work_.data(), T_work_.data(), 1);
    }

    scalar mutationMixture::H() const
    {
        return static_cast<scalar>(mixture_->mixtureHMass());
    }

    scalar mutationMixture::Cp() const
    {
        return static_cast<scalar>(mixture_->mixtureFrozenCpMass());
    }

    scalar mutationMixture::mu() const
    {
        return static_cast<scalar>(mixture_->viscosity());
    }

    scalar mutationMixture::kappa() const
    {
        return static_cast<scalar>(mixture_->frozenThermalConductivity());
    }

    // --- Helper Functions for Wrapper ---

    scalar mutationMixture::mixtureMw() const
    {
        return static_cast<scalar>(mixture_->mixtureMw());
    }

    scalar mutationMixture::Ru() const
    {
        return static_cast<scalar>(Mutation::RU);
    }

    // --- 2-Temperature Logic ---

    scalar mutationMixture::speciesEve(const label speciei) const
    {
        // Calculate Eve = E_total(T,Tve) - E_trans(T) - E_rot(T)

        // 1. Get Total Enthalpy (H/RT)
        // Uses T_work_[0] which is T_tr
        mixture_->speciesHOverRT(T_work_[0], H_RT_work_.data());

        const double T_tr = T_work_[0];
        const double Ru = Mutation::RU;
        const double Mw = Mw_work_[speciei];
        const double R_spec = Ru / Mw;

        // 2. Convert Enthalpy to Internal Energy (e = h - RT)
        const double h_spec = H_RT_work_[speciei] * R_spec * T_tr;
        const double e_total = h_spec - R_spec * T_tr;

        // 3. Translational Energy (1.5 R T)
        const double e_trans = 1.5 * R_spec * T_tr;

        // 4. Rotational Energy (Heuristic)
        double e_rot = 0.0;

        // Safe check for atom count using species name heuristic or size
        // (Mutation++ nAtoms() is safer if available, but this fallback works)
        if (species_.size() > static_cast<label>(speciei))
        {
            const word sname = species_[speciei];

            // Assume 1-char names (N, O) are atoms -> e_rot = 0
            // Assume names > 1 char (N2, NO, NO+) are molecules
            if (sname.size() > 1)
            {
                // Linear molecule approximation (Rigid Rotor)
                e_rot = 1.0 * R_spec * T_tr;
            }
        }

        const double e_ve = e_total - e_trans - e_rot;
        return static_cast<scalar>(e_ve);
    }

} // namespace Foam