#ifndef mutationMixture_H
#define mutationMixture_H

#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "speciesTable.H"

#include <vector>

// Forward Declaration
namespace Mutation
{
    class Mixture;
}

namespace Foam
{

    class mutationMixture
    {
        // Private Data
        const fvMesh &mesh_;
        speciesTable species_;
        PtrList<volScalarField> Y_;

        Mutation::Mixture *mixture_;
        word mechanismName_;

        // Working arrays (using std::vector for memory safety)
        mutable std::vector<double> Y_work_;
        mutable std::vector<double> T_work_;
        mutable std::vector<double> rho_work_;
        mutable std::vector<double> Mw_work_;
        mutable std::vector<double> E_work_;
        mutable std::vector<double> H_RT_work_;

    public:
        TypeName("mutationMixture");

        // Constructor
        mutationMixture(
            const dictionary &thermoDict,
            const fvMesh &mesh,
            const word &phaseName);

        // Destructor
        virtual ~mutationMixture();

        // Accessors
        const fvMesh &mesh() const { return mesh_; }
        const speciesTable &species() const { return species_; }
        const PtrList<volScalarField> &Y() const { return Y_; }
        PtrList<volScalarField> &Y() { return Y_; }

        // Return raw pointer if needed, though wrappers below are preferred
        Mutation::Mixture &mixture() { return *mixture_; }

        // --- Core Physics Function ---
        void setState(
            const scalar p,
            const scalar T,
            const scalar Tve,
            const scalarField &Y) const;

        // --- Property Accessors ---
        scalar H() const;
        scalar Cp() const;
        scalar mu() const;
        scalar kappa() const;

        // --- Helper Functions (Crucial for your Thermo wrapper) ---
        // Return the mixture mean molecular weight [kg/mol]
        scalar mixtureMw() const;

        // Return the Universal Gas Constant [J/(mol K)]
        scalar Ru() const;

        // --- 2-Temperature Specific ---
        // Vibrational + electronic energy for species (mass basis J/kg)
        scalar speciesEve(const label speciei) const;
    };

} // namespace Foam

#endif