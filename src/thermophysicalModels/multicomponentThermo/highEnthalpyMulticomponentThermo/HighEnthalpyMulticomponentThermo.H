#ifndef HighEnthalpyMulticomponentThermo_H
#define HighEnthalpyMulticomponentThermo_H

#include "fluidMulticomponentThermo.H"
#include "multicomponentThermo.H"
#include "psiThermo.H"
#include "basicThermo.H"
#include "fvMesh.H"

// Temporarily undefine OpenFOAM's log macro to avoid collision with Mutation++'s Log class
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

#include "mutationMixture.H"

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
    };

} // End namespace Foam

#endif

// ************************************************************************* //
