/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::highEnthalpyMulticomponentThermo

Description
    Class for multi-component high-enthalpy fluid thermodynamic properties
    based on compressibility.

See also
    Foam::basicThermo

SourceFiles
    highEnthalpyMulticomponentThermo.C

\*---------------------------------------------------------------------------*/

#ifndef highEnthalpyMulticomponentThermo_H
#define highEnthalpyMulticomponentThermo_H

#include "HighEnthalpyMulticomponentThermo.H"
#include "psiThermo.H"
#include "fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                          Class highEnthalpyMulticomponentThermo Declaration
    \*---------------------------------------------------------------------------*/

    class highEnthalpyMulticomponentThermo
        : virtual public psiThermo,
          virtual public fluidMulticomponentThermo
    {
    public:
        // Public Classes

        //- Forward declare the composite class
        class composite;

        virtual volScalarField &Tve() = 0;
        virtual const volScalarField &Tve() const = 0;
        virtual void correctTve(const volScalarField &newTve) = 0;

        // Public Typedefs

        //- The derived type
        template <class MixtureType>
        using DerivedThermoType =
            HighEnthalpyMulticomponentThermo<
                MulticomponentThermo<
                    PsiThermo<
                        BasicThermo<MixtureType, composite>>>>;

        //- Runtime type information
        TypeName("highEnthalpyMulticomponentThermo");

        //- Declare run-time constructor selection tables
        declareRunTimeSelectionTable(
            autoPtr,
            highEnthalpyMulticomponentThermo,
            fvMesh,
            (const fvMesh &mesh, const word &phaseName),
            (mesh, phaseName));

        // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<highEnthalpyMulticomponentThermo> New(
            const fvMesh &,
            const word &phaseName = word::null);

        //- The derived name
        static word derivedThermoName()
        {
            return "highEnthalpyThermo";
        };

        //- Destructor
        virtual ~highEnthalpyMulticomponentThermo();
    };

    /*---------------------------------------------------------------------------*\
                    Class highEnthalpyMulticomponentThermo::composite Declaration
    \*---------------------------------------------------------------------------*/

    class highEnthalpyMulticomponentThermo::composite
        : public basicThermo::implementation,
          public fluidThermo::implementation,
          public psiThermo::implementation,
          public multicomponentThermo::implementation,
          public highEnthalpyMulticomponentThermo
    {

    private:
        //- Mutation mixture for two-temperature relaxation
        autoPtr<mutationMixture> mutationMixPtr_;

        //- Time step for energy relaxation (read from dict or default)
        scalar relaxationTimeStep_;

    public:
        volScalarField Tve_;

        // Constructors

        //- Construct from dictionary, mesh and phase name
        template <class MixtureType>
        composite(
            const dictionary &dict,
            const MixtureType &mixture,
            const fvMesh &mesh,
            const word &phaseName)
            : basicThermo::implementation(dict, mesh, phaseName),
              fluidThermo::implementation(dict, mesh, phaseName),
              psiThermo::implementation(dict, mesh, phaseName),
              multicomponentThermo::implementation(
                  dict,
                  mixture.specieNames(),
                  mesh,
                  phaseName),
              Tve_(
                  IOobject(
                      "Tve",
                      mesh.time().name(),
                      mesh,
                      IOobject::MUST_READ,
                      IOobject::AUTO_WRITE),
                  mesh)
        {

            Info << "highEnthalpyMulticomponentThermo::composite - "
                 << "constructor"
                 << endl;

            // Read mixture name from thermophysicalProperties dictionary
            word mixtureName = "air_5"; // default
            if (dict.found("mixture"))
            {
                mixtureName = word(dict.lookup("mixture"));
            }

            Info << "Initializing mutationMixture with mechanism: "
                 << mixtureName << endl;

            mutationMixPtr_.reset(
                new mutationMixture(mixtureName.c_str()));

            // Read relaxation time step (default: 1e-8)
            relaxationTimeStep_ = dict.lookupOrDefault<scalar>("relaxationTimeStep", 1.0e-8);

            Info << "Relaxation time step: " << relaxationTimeStep_ << " s" << endl;
        }

        volScalarField &Tve() override
        {
            return Tve_;
        }

        const volScalarField &Tve() const override
        {
            return Tve_;
        }

        void correctTve(const volScalarField &newTve) override
        {
            // 1) Basic safety checks (optional but useful)
            if (newTve.size() != Tve_.size())
            {
                FatalErrorInFunction
                    << "Tve size mismatch: newTve=" << newTve.size()
                    << " Tve_=" << Tve_.size()
                    << exit(FatalError);
            }

            if (newTve.dimensions() != Tve_.dimensions())
            {
                FatalErrorInFunction
                    << "Tve dimensions mismatch: newTve=" << newTve.dimensions()
                    << " Tve_=" << Tve_.dimensions()
                    << exit(FatalError);
            }

            // 2) Copy INTERNAL values (raw)
            forAll(Tve_.internalFieldRef(), i)
            {
                Tve_.internalFieldRef()[i] = newTve.internalField()[i];
            }

            // 3) Copy BOUNDARY values (raw)
            forAll(Tve_.boundaryFieldRef(), patchi)
            {
                forAll(Tve_.boundaryFieldRef()[patchi], facei)
                {
                    Tve_.boundaryFieldRef()[patchi][facei] =
                        newTve.boundaryField()[patchi][facei];
                }
            }

            // 4) Enforce BC consistency
            Tve_.correctBoundaryConditions();
        }

        // Override correct() to use mutationMixture
        virtual void correct() override
        {
            Info << "highEnthalpyMulticomponentThermo::composite::correct() called" << endl;

            volScalarField &T = this->T_;
            volScalarField &Tve = this->Tve_;
            volScalarField &p = this->p_;

            // rho() returns tmp<volScalarField> in many implementations
            tmp<volScalarField> trho = this->rho();
            const volScalarField &rho = trho();

            // Species count (volScalarField list)
            const label nSpecies = this->Y().size();

            // Universal gas constant
            static constexpr double Ru = 8.31446261815324;

            // Build a slicer for Y fields (OpenFOAM pattern you already use)
            auto Yslicer = this->Yslicer();

            forAll(T, celli)
            {
                // Grab mass fractions for this cell
                std::vector<double> Y_cell(nSpecies, 0.0);
                for (label i = 0; i < nSpecies; ++i)
                    Y_cell[i] = this->Y()[i][celli];

                // Compute Rmix = sum(Yi * Ru/Mwi)
                double Rmix = 0.0;
                for (label i = 0; i < nSpecies; ++i)
                {
                    if (Y_cell[i] <= 0.0)
                        continue;
                    const scalar Mw = this->Wi(i).value(); // kg/mol
                    if (Mw > 0.0)
                        Rmix += Y_cell[i] * Ru / Mw;
                }

                if (Rmix <= SMALL || rho[celli] <= SMALL)
                {
                    // fallback
                    Tve[celli] = T[celli];
                    continue;
                }

                // Build energies from current temperatures
                // (translation-only model consistent with your wrapper)
                double Et = rho[celli] * 1.5 * Rmix * T[celli];

                // IMPORTANT: you previously used a "simplified" Ev.
                // Better: store Ev as a field later.
                // For now, keep a conservative "Ev = 0" or use previous cell Ev if you add it.
                double Ev = 0.0;

                double Ttr = T[celli];
                double Tv = Tve[celli];

                if (mutationMixPtr_.valid())
                {
                    mutationMixPtr_->step(
                        relaxationTimeStep_,
                        rho[celli],
                        Y_cell,
                        Et,
                        Ev,
                        Ttr,
                        Tv);

                    T[celli] = Ttr;
                    Tve[celli] = Tv;

                    // Optional: keep pressure consistent
                    p[celli] = rho[celli] * Rmix * T[celli];
                }
                else
                {
                    // fallback
                    T[celli] = p[celli] / (rho[celli] * Rmix);
                    Tve[celli] = T[celli];
                }
            }

            T.correctBoundaryConditions();
            Tve.correctBoundaryConditions();
            p.correctBoundaryConditions();
        }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
