#ifndef psiThermoMixture_H
#define psiThermoMixture_H

#include "PsiThermo.H"

namespace Foam
{

    template <class Thermo>
    class psiThermoMixture : public PsiThermo<Thermo>
    {
    public:
        // Typedef for convenience
        using BaseThermo = PsiThermo<Thermo>;

        // Constructor
        psiThermoMixture(const fvMesh &mesh, const word &phaseName)
            : PsiThermo<Thermo>(mesh, phaseName) {}

        // Destructor
        virtual ~psiThermoMixture() {}

        // Overload correct() instead of calculate()
        virtual void correct()
        {
            // Call base class correct() if exists
            BaseThermo::correct();

            Info << "psiThermoMixture::correct() called" << endl;
            // Add any mixture-specific correction here
        }
    };

} // End namespace Foam

#endif
