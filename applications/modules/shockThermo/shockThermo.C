#include "shockThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "fvModels.H"
#include "fvConstraints.H"

namespace Foam
{
    namespace solvers
    {
        // 1. Register the solver
        defineTypeNameAndDebug(shockThermo, 0);
        addToRunTimeSelectionTable(solver, shockThermo, fvMesh);

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        shockThermo::shockThermo(fvMesh &mesh)
            : // Initialize Base Class (Creates the thermo from dictionary)
              shockFluidThermo(mesh, autoPtr<fluidThermo>(nullptr)),

              // Initialize References
              thermo_(refCast<fluidMulticomponentThermo>(thermoPtr_())),
              Y_(thermo_.Y()),

              // Initialize Models
              reaction(combustionModel::New(thermo_, momentumTransport())),

              thermophysicalTransport(
                  fluidMulticomponentThermophysicalTransportModel::New(
                      momentumTransport(),
                      thermo_)),

              // --- Initialize Wrapper ---
              // We pass the reference 'thermo_' to your wrapper
              heThermo_(new highEnthalpyMulticomponentThermo(thermo_)),

              // Initialize Public References
              thermo(thermo_),
              Y(Y_)
        {
            // Validation
            thermo.validate(type(), "h", "e");

            // Collect fields for interpolation
            forAll(Y_, i)
            {
                fields.add(Y_[i]);
            }
            fields.add(thermo.he());
        }

        // Destructor
        shockThermo::~shockThermo() = default;

        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        void shockThermo::preSolve()
        {
            shockFluidThermo::preSolve();
        }

        void shockThermo::thermophysicalPredictor()
        {
            // 1. Convection scheme
            tmp<fv::convectionScheme<scalar>> mvConvection(
                fv::convectionScheme<scalar>::New(
                    mesh,
                    fields,
                    phi,
                    mesh.schemes().div("div(phi,Yi_h)")));

            // 2. Reaction
            reaction->correct();

            // 3. Species Transport
            forAll(Y_, i)
            {
                volScalarField &Yi = Y_[i];
                if (thermo_.solveSpecie(i))
                {
                    fvScalarMatrix YiEqn(
                        fvm::ddt(rho, Yi) + mvConvection->fvmDiv(phi, Yi) + thermophysicalTransport->divj(Yi) ==
                        reaction->R(Yi) + fvModels().source(rho, Yi));

                    YiEqn.relax();
                    fvConstraints().constrain(YiEqn);
                    YiEqn.solve("Yi");
                    fvConstraints().constrain(Yi);
                }
            }

            // 4. Normalise
            thermo_.normaliseY();

            // 5. Energy Equation setup
            volScalarField &e = thermo_.he();

            // Interpolate energy
            const surfaceScalarField e_pos(interpolate(e, pos, thermo_.T().name()));
            const surfaceScalarField e_neg(interpolate(e, neg, thermo_.T().name()));

            // Convective flux
            surfaceScalarField phiEp(
                "phiEp",
                aphiv_pos() * (rho_pos() * (e_pos + 0.5 * magSqr(U_pos())) + p_pos()) + aphiv_neg() * (rho_neg() * (e_neg + 0.5 * magSqr(U_neg())) + p_neg()) + aSf() * (p_pos() - p_neg()));

            if (mesh.moving())
            {
                phiEp += mesh.phi() * (a_pos() * p_pos() + a_neg() * p_neg());
            }

            fvScalarMatrix EEqn(
                fvm::ddt(rho, e) + fvc::div(phiEp) + fvm::ddt(rho, K) ==
                fvModels().source(rho, e));

            if (!inviscid)
            {
                const surfaceScalarField devTauDotU(
                    "devTauDotU",
                    devTau() & (a_pos() * U_pos() + a_neg() * U_neg()));
                EEqn += thermophysicalTransport->divq(e) + fvc::div(devTauDotU);
            }

            // 6. Solve Energy
            EEqn.relax();
            fvConstraints().constrain(EEqn);
            EEqn.solve();
            fvConstraints().constrain(e);

            // 7. Update Thermodynamics
            // --- HERE IS THE MAGIC ---
            // Instead of thermo_.correct(), we call the wrapper!
            heThermo_->correct();
        }

    } // namespace solvers
} // namespace Foam