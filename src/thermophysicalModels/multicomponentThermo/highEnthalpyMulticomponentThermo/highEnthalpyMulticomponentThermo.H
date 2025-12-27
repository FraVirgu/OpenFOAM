/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
\*---------------------------------------------------------------------------*/

#ifndef highEnthalpyMulticomponentThermo_H
#define highEnthalpyMulticomponentThermo_H

#include "HighEnthalpyMulticomponentThermo.H"
#include "psiThermo.H"
#include "fluidMulticomponentThermo.H"

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
        class composite;

        virtual volScalarField &Tve() = 0;
        virtual const volScalarField &Tve() const = 0;
        virtual void correctTve(const volScalarField &newTve) = 0;
        virtual void correct_he() = 0;

        template <class MixtureType>
        using DerivedThermoType =
            HighEnthalpyMulticomponentThermo<
                MulticomponentThermo<
                    PsiThermo<
                        BasicThermo<MixtureType, composite>>>>;

        TypeName("highEnthalpyMulticomponentThermo");

        declareRunTimeSelectionTable(
            autoPtr,
            highEnthalpyMulticomponentThermo,
            fvMesh,
            (const fvMesh &mesh, const word &phaseName),
            (mesh, phaseName));

        static autoPtr<highEnthalpyMulticomponentThermo> New(
            const fvMesh &,
            const word &phaseName = word::null);

        static word derivedThermoName() { return "highEnthalpyThermo"; };

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
        autoPtr<mutationMixture> mutationMixPtr_;
        scalar relaxationTimeStep_;

    public:
        volScalarField Tve_;
        volScalarField Et_;
        volScalarField Ev_;

        // Constructors
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
                  mesh),
              Et_(
                  IOobject("Et", mesh.time().name(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                  mesh,
                  dimensionedScalar(dimEnergy / dimVolume, 0.0)),
              Ev_(
                  IOobject("Ev", mesh.time().name(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                  mesh,
                  dimensionedScalar(dimEnergy / dimVolume, 0.0))
        {
            Info << "highEnthalpyMulticomponentThermo::composite - constructor" << endl;

            word mixtureName = "air_5";
            if (dict.found("mixture"))
                mixtureName = word(dict.lookup("mixture"));

            Info << "Initializing mutationMixture with mechanism: " << mixtureName << endl;
            mutationMixPtr_.reset(new mutationMixture(mixtureName.c_str()));

            relaxationTimeStep_ = dict.lookupOrDefault<scalar>("relaxationTimeStep", 1.0e-8);
            Info << "Relaxation time step: " << relaxationTimeStep_ << " s" << endl;

            // ------------------------------------------------------------
            // Initialization Logic (Using Mutation++ for Mw to avoid Link Error)
            // ------------------------------------------------------------
            if (max(mag(Et_)).value() <= SMALL)
            {
                Info << "Initializing Et/Ev fields from Temperature..." << endl;

                volScalarField &T = this->T_;
                tmp<volScalarField> trho = this->rho();
                const volScalarField &rho = trho();
                const label nSpecies = this->Y().size();
                static constexpr double Ru = 8.31446261815324;

                forAll(Et_, celli)
                {
                    double Rmix = 0.0;
                    for (label i = 0; i < nSpecies; ++i)
                    {
                        double Yi = this->Y()[i][celli];
                        if (Yi > SMALL)
                        {
                            // FIX: Use Mutation++ Mw (kg/mol) instead of OpenFOAM Wi
                            // This bypasses the linker error and matches Ru units.
                            scalar Mw = mutationMixPtr_->speciesMw(i);
                            Rmix += Yi * Ru / Mw;
                        }
                    }

                    if (Rmix > SMALL && rho[celli] > SMALL)
                        Et_[celli] = rho[celli] * 1.5 * Rmix * T[celli];
                    else
                        Et_[celli] = 1.0;

                    Ev_[celli] = 0.0;
                }
            }
        }

        volScalarField &Tve() override { return Tve_; }
        const volScalarField &Tve() const override { return Tve_; }

        void correctTve(const volScalarField &newTve) override
        {
            Tve_ = newTve;
            Tve_.correctBoundaryConditions();
        }

        virtual void correct_he() override
        {
            Info << "highEnthalpyMulticomponentThermo::composite::correct() called" << endl;

            volScalarField &T = this->T_;
            volScalarField &Tve = this->Tve_;
            volScalarField &p = this->p_;
            volScalarField &e = this->he();

            tmp<volScalarField> trho = this->rho();
            const volScalarField &rho = trho();

            const label nSpecies = this->Y().size();
            static constexpr double Ru = 8.31446261815324;

            forAll(T, celli)
            {
                std::vector<double> Y_cell(nSpecies, 0.0);
                for (label i = 0; i < nSpecies; ++i)
                    Y_cell[i] = this->Y()[i][celli];

                double Rmix = 0.0;
                for (label i = 0; i < nSpecies; ++i)
                {
                    if (Y_cell[i] <= 0.0)
                        continue;
                    // FIX: Use Mutation++ Mw here too
                    const scalar Mw = mutationMixPtr_->speciesMw(i);
                    if (Mw > 0.0)
                        Rmix += Y_cell[i] * Ru / Mw;
                }

                if (Rmix <= SMALL || rho[celli] <= SMALL)
                {
                    Tve[celli] = T[celli];
                    continue;
                }

                double Etot_new = rho[celli] * e[celli];
                double Et_prev = Et_[celli];
                double Ev_prev = Ev_[celli];

                if (Et_prev + Ev_prev <= SMALL)
                {
                    Et_prev = rho[celli] * 1.5 * Rmix * T[celli];
                    Ev_prev = 0.0;
                }

                double scale = 1.0;
                if ((Et_prev + Ev_prev) > SMALL)
                {
                    scale = Etot_new / (Et_prev + Ev_prev);
                }

                double Et_input = Et_prev * scale;
                double Ev_input = Ev_prev * scale;

                double Ttr = T[celli];
                double Tv = Tve[celli];

                if (mutationMixPtr_.valid())
                {
                    mutationMixPtr_->step(
                        relaxationTimeStep_,
                        rho[celli],
                        Y_cell,
                        Et_input,
                        Ev_input,
                        Ttr,
                        Tv);

                    T[celli] = Ttr;
                    Tve[celli] = Tv;
                    p[celli] = rho[celli] * Rmix * T[celli];
                    e[celli] = (Et_input + Ev_input) / rho[celli];

                    Et_[celli] = Et_input;
                    Ev_[celli] = Ev_input;
                }
                else
                {
                    T[celli] = p[celli] / (rho[celli] * Rmix);
                    Tve[celli] = T[celli];
                }
            }

            T.correctBoundaryConditions();
            Tve.correctBoundaryConditions();
            p.correctBoundaryConditions();
            Et_.correctBoundaryConditions();
            Ev_.correctBoundaryConditions();
        }
    };
} // End namespace Foam

#endif