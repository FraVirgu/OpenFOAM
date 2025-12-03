#ifndef highEnthalpyMulticomponentThermo_H
#define highEnthalpyMulticomponentThermo_H

#include "fluidMulticomponentThermo.H"
#include "mutationMixture.H"
#include "volFields.H"
#include "IOdictionary.H"

namespace Foam
{
    class highEnthalpyMulticomponentThermo
    {
    private:
        // 1. Reference to the base thermo (fields, mesh)
        // We use a Reference (&) because shockFluidThermo owns the object.
        fluidMulticomponentThermo &wrappedThermo_;

        // 2. Reference to the specific Mutation++ interface (physics)
        mutationMixture &mutationMix_;

        // 3. The Vibrational Temperature Field
        volScalarField Tve_;

        // Helper function
        scalar readEnthalpyCorrection() const;

    public:
        // Define TypeName for debug switches
        ClassName("highEnthalpyMulticomponentThermo");

        // Constructor: Takes BOTH the thermo and the mutation mixture
        highEnthalpyMulticomponentThermo(
            fluidMulticomponentThermo &wrapped,
            mutationMixture &mix);

        // Destructor
        ~highEnthalpyMulticomponentThermo() = default;

        // Accessors
        fluidMulticomponentThermo &wrappedThermo() { return wrappedThermo_; }
        mutationMixture &mutation() { return mutationMix_; }

        // Field Access / Forwarding
        const volScalarField &T() const;
        volScalarField &T();

        const volScalarField &p() const;
        volScalarField &p();

        const PtrList<volScalarField> &Y() const;
        PtrList<volScalarField> &Y();

        const volScalarField &Cp() const;
        const volScalarField &Cv() const;

        // This accesses the Internal Energy field 'e' in shockThermo
        const volScalarField &he() const;
        volScalarField &he();

        // Accessor for the Vibrational Temperature
        const volScalarField &Tve() const { return Tve_; }
        volScalarField &Tve() { return Tve_; }

        // The Magic Function: Updates Physics using T_tr AND T_ve
        void correct();
    };
}

#endif