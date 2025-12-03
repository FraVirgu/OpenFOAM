#include "psiThermoMixture.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

    // Example of runtime selection table for a specific Thermo type
    // You need to replace MyThermoType with your actual type
    typedef psiThermoMixture<MyThermoType> MyPsiThermoMixture;

    addToRunTimeSelectionTable(
        basicThermo,
        MyPsiThermoMixture,
        fvMesh);

} // End namespace Foam
