#ifndef HighEnthalpyMulticomponentThermo_H
#define HighEnthalpyMulticomponentThermo_H

#include "fluidMulticomponentThermo.H"
#include "multicomponentThermo.H"
#include "psiThermo.H"
#include "basicThermo.H"
#include "fvMesh.H"

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                Class HighEnthalpyMulticomponentThermo Declaration
    \*---------------------------------------------------------------------------*/

    template <class Thermo>
    class HighEnthalpyMulticomponentThermo
        : public FluidMulticomponentThermo<Thermo>
    {
    public:
        // Constructors
        HighEnthalpyMulticomponentThermo(
            const fvMesh &mesh,
            const word &phaseName = word::null)
            : FluidMulticomponentThermo<Thermo>(mesh, phaseName)
        {
        }

        virtual ~HighEnthalpyMulticomponentThermo() = default;

        // Member Functions
        virtual void correct() override
        {

            Info << "HighEnthalpyMulticomponentThermo::correct() called"
                 << endl;

            auto &T = this->T_;
            auto &p = this->p_;

            auto &Cp = this->Cp_;
            auto &Cv = this->Cv_;
            auto &psi = this->psi_;
            auto &mu = this->mu_;
            auto &kappa = this->kappa_;

            auto Yslicer = this->Yslicer();

            forAll(T, celli)
            {
                auto composition = this->cellComposition(Yslicer, celli);

                const auto &thermoMixture =
                    this->thermoMixture(composition);

                const auto &transportMixture =
                    this->transportMixture(composition, thermoMixture);

                Cp[celli] = thermoMixture.Cp(p[celli], T[celli]);
                Cv[celli] = thermoMixture.Cv(p[celli], T[celli]);
                psi[celli] = thermoMixture.psi(p[celli], T[celli]);

                mu[celli] = transportMixture.mu(p[celli], T[celli]);
                kappa[celli] = transportMixture.kappa(p[celli], T[celli]);
            }

            // Boundary conditions
            Cp.correctBoundaryConditions();
            Cv.correctBoundaryConditions();
            psi.correctBoundaryConditions();
            mu.correctBoundaryConditions();
            kappa.correctBoundaryConditions();
        }
    };

} // End namespace Foam

#endif

// ************************************************************************* //
