#include "PsiThermoMixture.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

    // Example of runtime selection table for a specific Thermo type
    typedef PsiThermoMixture<MyThermoType> MyPsiThermoMixtureType;

    addToRunTimeSelectionTable(
        basicThermo,
        MyPsiThermoMixtureType,
        fvMesh);

} // End namespace Foam
