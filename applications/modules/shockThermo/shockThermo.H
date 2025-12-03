#ifndef shockThermo_H
#define shockThermo_H

#include "shockFluidThermo.H"
#include "fluidMulticomponentThermo.H"
#include "combustionModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "multivariateScheme.H"

// INCLUDE YOUR WRAPPER HEADER
#include "highEnthalpyMulticomponentThermo.H"

namespace Foam
{
    namespace solvers
    {
        /*---------------------------------------------------------------------------*\
                                  Class shockThermo Declaration
        \*---------------------------------------------------------------------------*/

        class shockThermo
            : public shockFluidThermo
        {

        private:
            // Thermophysical properties
            fluidMulticomponentThermo &thermo_;

            // Composition
            PtrList<volScalarField> &Y_;

            // Reactions
            autoPtr<combustionModel> reaction;

            // Fields for interpolation
            multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

            // Thermophysical transport
            autoPtr<fluidMulticomponentThermophysicalTransportModel> thermophysicalTransport;

            // --- ADDED: The High Enthalpy Wrapper ---
            autoPtr<highEnthalpyMulticomponentThermo> heThermo_;

        public:
            // Public Data references (read-only access for others)
            const fluidMulticomponentThermo &thermo;
            const PtrList<volScalarField> &Y;

            // Runtime type information
            TypeName("shockThermo");

            // Constructors

            //- Construct from region mesh
            explicit shockThermo(fvMesh &mesh);

            //- Disallow default bitwise copy construction
            shockThermo(const shockThermo &) = delete;

            //- Destructor
            virtual ~shockThermo();

            // Member Functions

            //- Called at the start of the time-step
            virtual void preSolve();

            //- Solve Energy Equation
            virtual void thermophysicalPredictor();

            // Member Operators
            void operator=(const shockThermo &) = delete;
        };
    } // End namespace solvers
} // End namespace Foam

#endif