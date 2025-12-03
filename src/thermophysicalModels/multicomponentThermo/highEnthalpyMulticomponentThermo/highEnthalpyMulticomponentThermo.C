#include "highEnthalpyMulticomponentThermo.H"
#include "dimensionedScalar.H"
#include "IOdictionary.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(highEnthalpyMulticomponentThermo, 0);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    highEnthalpyMulticomponentThermo::highEnthalpyMulticomponentThermo(
        fluidMulticomponentThermo &wrapped,
        mutationMixture &mix)
        : wrappedThermo_(wrapped),
          mutationMix_(mix),
          Tve_(
              IOobject(
                  "Tve",
                  wrapped.T().time().name(),
                  wrapped.T().mesh(),
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE),
              wrapped.T().mesh())
    {
        Info << "highEnthalpyMulticomponentThermo: Linked to "
             << wrappedThermo_.thermoName()
             << " and Mutation++ mechanism." << endl;
    }

    // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    // Forwarding
    const volScalarField &highEnthalpyMulticomponentThermo::T() const { return wrappedThermo_.T(); }
    volScalarField &highEnthalpyMulticomponentThermo::T() { return wrappedThermo_.T(); }

    const volScalarField &highEnthalpyMulticomponentThermo::p() const { return wrappedThermo_.p(); }
    volScalarField &highEnthalpyMulticomponentThermo::p() { return wrappedThermo_.p(); }

    const PtrList<volScalarField> &highEnthalpyMulticomponentThermo::Y() const { return wrappedThermo_.Y(); }
    PtrList<volScalarField> &highEnthalpyMulticomponentThermo::Y() { return wrappedThermo_.Y(); }

    const volScalarField &highEnthalpyMulticomponentThermo::Cp() const { return wrappedThermo_.Cp(); }
    const volScalarField &highEnthalpyMulticomponentThermo::Cv() const { return wrappedThermo_.Cv(); }

    const volScalarField &highEnthalpyMulticomponentThermo::he() const { return wrappedThermo_.he(); }
    volScalarField &highEnthalpyMulticomponentThermo::he() { return wrappedThermo_.he(); }

    scalar highEnthalpyMulticomponentThermo::readEnthalpyCorrection() const
    {
        return 0.0;
    }

    // ------------------------------------------------------------------------
    // THE 2-TEMPERATURE UPDATE LOOP
    // ------------------------------------------------------------------------
    void highEnthalpyMulticomponentThermo::correct()
    {
        // 1. Run standard update (Calculates T from current energy, updates transport)
        wrappedThermo_.correct();

        // 2. OVERWRITE Energy with Non-Equilibrium Physics (T_tr != T_ve)
        Info << "HighEnthalpy: Updating 2-T Physics via Mutation++..." << endl;

        const volScalarField &p = wrappedThermo_.p();
        const volScalarField &Ttr = wrappedThermo_.T();
        const PtrList<volScalarField> &Y = wrappedThermo_.Y();

        // Access the Energy Field (This is 'e' internal energy in shockThermo)
        volScalarField &e_field = wrappedThermo_.he();

        // NOTE: We do NOT update Cp, mu, alpha here.
        // In the Wrapper pattern, these are read-only from the base thermo.
        // We rely on the base models (Janaf/Sutherland) for transport,
        // but we enforce strict Energy Conservation using Mutation++.

        forAll(Ttr, celli)
        {
            // Gather Data
            scalar P_c = p[celli];
            scalar T_c = Ttr[celli];
            scalar Tve_c = Tve_[celli];

            // Build local mass fractions
            scalarField Y_c(Y.size());
            forAll(Y, i) Y_c[i] = Y[i][celli];

            // CALL MUTATION++
            // Sets the state using T_tr AND T_ve
            mutationMix_.setState(P_c, T_c, Tve_c, Y_c);

            // GET ENTHALPY [J/kg]
            double H_mut = mutationMix_.H();

            // CONVERT TO INTERNAL ENERGY [J/kg]
            // e = h - p/rho = h - R_mixture * T
            // R_mixture = Ru / M_mixture

            double Mw_mix = mutationMix_.mixtureMw(); // Calls your new helper
            double Ru = mutationMix_.Ru();            // Calls your new helper
            double R_mix = Ru / Mw_mix;               // J/(kg K)

            // Update the OpenFOAM field
            e_field[celli] = H_mut - (R_mix * T_c);
        }

        // 3. Correct Boundaries
        e_field.correctBoundaryConditions();
        Tve_.correctBoundaryConditions();
    }
}