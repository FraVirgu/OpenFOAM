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
        // Ensure we use Vibrational Temperature (T_work_[1])
        double T_use = T_work_[1];

        // Get enthalpy at Tve
        mixture_->speciesHOverRT(T_use, H_RT_work_.data());

        const double Ru = Mutation::RU;
        const double Mw = Mw_work_[speciei];
        const double R_spec = Ru / Mw;

        // e = h - RT
        const double h_spec = H_RT_work_[speciei] * R_spec * T_use;
        const double e_tot = h_spec - R_spec * T_use;

        // Subtract Translational (1.5RT) and Rotational (1.0RT for molecules)
        double e_trans = 1.5 * R_spec * T_use;
        double e_rot = 0.0;

        // Simple heuristic: If name length > 1 (e.g., N2, NO), it has rotation
        if (species_[speciei].size() > 1)
        {
            e_rot = 1.0 * R_spec * T_use;
        }

        // Result is Vibrational + Electronic energy
        return static_cast<scalar>(max(0.0, e_tot - e_trans - e_rot));
    }

    // 2. NEW HELPER: Get Source Term Qve [J/m^3/s]
    scalar mutationMixture::getVibrationalSource() const
    {
        // Vector size depends on # of energy equations (usually 1 for Tve)
        std::vector<double> sources(mixture_->nEnergyEqns());

        // Get source terms (Relaxation + Chemistry)
        mixture_->energyTransferSource(sources.data());

        // Return the first component (Vib-Electronic source)
        return static_cast<scalar>(sources[0]);
    }

    // 3. NEW HELPER: Invert Energy -> Temperature (Newton-Raphson)
    scalar mutationMixture::solveTve(
        const scalar rho,
        const scalar eve_target,
        const scalarField &Y,
        const scalar Tve_guess) const
    {
        double T_new = Tve_guess;
        const int maxIter = 20;
        const double tol = 1.0e-5;

        for (int iter = 0; iter < maxIter; ++iter)
        {
            // Set state with T_tr=T_new and T_ve=T_new to calculate properites at this temp
            // (We only care about the vibrational part here)
            this->setState(101325.0, T_new, T_new, Y);

            // Calculate calculated Energy and Cv at this temperature
            double e_calc = 0.0;
            double cv_calc = 0.0;

            // Finite Difference for Cv (simplest robust method)
            double dT = 1.0;
            for (int k = 0; k < Y.size(); ++k)
                e_calc += Y[k] * speciesEve(k);

            // Perturb T
            T_work_[1] += dT;
            double e_plus = 0.0;
            for (int k = 0; k < Y.size(); ++k)
                e_plus += Y[k] * speciesEve(k);

            cv_calc = (e_plus - e_calc) / dT;

            // Newton Step
            if (std::abs(cv_calc) < 1e-10)
                break; // Avoid div by 0

            double delta = (e_calc - eve_target) / cv_calc;

            // Limit step size for stability
            delta = std::max(-500.0, std::min(500.0, delta));

            T_new -= delta;
            T_new = std::max(300.0, std::min(20000.0, T_new)); // Clamp

            if (std::abs(delta / T_new) < tol)
                return T_new;
        }
        return T_new;
    }

    // 4. NEW HELPER: Solve T (Translational)
    scalar mutationMixture::solveT(
        const scalar rho,
        const scalar e_tr_target,
        const scalarField &Y,
        const scalar T_guess) const
    {
        double T_new = T_guess;
        const int maxIter = 20;
        const double tol = 1.0e-6;
        const double Ru = Mutation::RU; // Universal Gas Constant

        for (int iter = 0; iter < maxIter; ++iter)
        {
            double T_old = T_new;

            double e_calc = 0.0;
            double cv_calc = 0.0;

            // --- Loop over species to sum up Energy and Cv ---
            for (int k = 0; k < Y.size(); ++k)
            {
                // Specific Gas Constant for this species: R = Ru / M
                double R_spec = Ru / Mw_work_[k];

                // Degrees of Freedom Factor:
                // Translation = 1.5
                // Rotation    = 1.0 (only for molecules)
                double factor = 1.5;

                // Heuristic: If name length > 1 (e.g. "N2", "NO"), it is a molecule
                if (species_[k].size() > 1)
                {
                    factor += 1.0;
                }

                double cv_spec = factor * R_spec;

                // Accumulate Mixture properties
                // e = sum( Y_i * Cv_i * T )
                e_calc += Y[k] * cv_spec * T_old;
                cv_calc += Y[k] * cv_spec;
            }

            // --- Newton-Raphson Step ---
            // f(T) = e_calc(T) - e_target = 0
            // f'(T) = cv_calc
            // delta = f(T) / f'(T)

            if (cv_calc < 1.0e-12)
                break; // Avoid division by zero

            double delta = (e_calc - e_tr_target) / cv_calc;

            T_new = T_old - delta;

            // --- Safety Clamps ---
            // Prevent T from going negative or exploding during iteration
            T_new = std::max(200.0, std::min(25000.0, T_new));

            // --- Convergence Check ---
            if (std::abs(delta / T_new) < tol)
            {
                return static_cast<scalar>(T_new);
            }
        }

        return static_cast<scalar>(T_new);
    }
} // namespace Foam