#ifndef PsiThermoMixture_H
#define PsiThermoMixture_H

#include "psiThermoMixture.H"

namespace Foam
{

    template <class Thermo>
    class PsiThermoMixture : public psiThermoMixture<Thermo>
    {
    public:
        // Typedef for base
        using Base = psiThermoMixture<Thermo>;

        // Constructor
        PsiThermoMixture(const fvMesh &mesh, const word &phaseName)
            : Base(mesh, phaseName) {}

        // Destructor
        virtual ~PsiThermoMixture() {}

        // Optionally, override correct() again
        virtual void correct() override
        {
            Base::correct();
            Info << "PsiThermoMixture::correct() called" << endl;
        }
    };

} // End namespace Foam

#ifdef NoRepository
#include "PsiThermoMixture.C"
#endif

#endif
