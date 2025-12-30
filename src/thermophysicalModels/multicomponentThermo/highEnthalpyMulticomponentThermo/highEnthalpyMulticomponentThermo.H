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
#include <cmath>

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
        List<label> ofToMut_;

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
            const wordList &ofNames = mixture.specieNames(); // OpenFOAM species names
            ofToMut_.setSize(ofNames.size(), -1);

            forAll(ofNames, iOF)
            {
                const std::string n = ofNames[iOF];
                const int iMut = mutationMixPtr_->speciesIndex(n); // your wrapper has this in the test

                if (iMut < 0)
                {
                    WarningInFunction
                        << "OpenFOAM species " << ofNames[iOF]
                        << " not found in Mutation++ mechanism " << mixtureName
                        << " -> it will be ignored in VT step." << nl;
                }
                else
                {
                    ofToMut_[iOF] = iMut;
                }
            }

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
            Info << "highEnthalpyMulticomponentThermo::composite::correct_he() called" << endl;

            volScalarField &T = this->T_;
            volScalarField &Tve = this->Tve_;
            volScalarField &p = this->p_;
            volScalarField &e = this->he();

            tmp<volScalarField> trho = this->rho();
            const volScalarField &rho = trho();

            const label nsOF = this->Y().size();
            const label nsMut = mutationMixPtr_->nSpecies();

            auto finite = [](double x)
            { return std::isfinite(x); };

            forAll(T, celli)
            {
                if (!mutationMixPtr_.valid())
                {
                    // Fallback: just keep 1T
                    Tve[celli] = T[celli];
                    continue;
                }

                if (rho[celli] <= SMALL || !finite(rho[celli]))
                {
                    Tve[celli] = T[celli];
                    continue;
                }

                // --- Build Y in Mutation++ ordering
                std::vector<double> Y_mut(nsMut, 0.0);

                for (label iOF = 0; iOF < nsOF; ++iOF)
                {
                    const label iMut = ofToMut_[iOF];
                    if (iMut < 0)
                        continue;

                    const double y = this->Y()[iOF][celli];
                    if (finite(y))
                        Y_mut[iMut] += y;
                }

                // --- Clip/renormalize Y_mut
                double sumY = 0.0;
                for (auto &y : Y_mut)
                {
                    if (!finite(y))
                        y = 0.0;
                    y = std::max(y, 0.0);
                    sumY += y;
                }

                if (sumY <= SMALL || !finite(sumY))
                {
                    Tve[celli] = T[celli];
                    continue;
                }

                for (auto &y : Y_mut)
                    y /= sumY;

                // --- Compute Rmix using Mutation++ MWs (Mutation ordering)
                static constexpr double Ru = 8.31446261815324; // J/mol/K
                double Rmix = 0.0;

                for (label iMut = 0; iMut < nsMut; ++iMut)
                {
                    const double Yi = Y_mut[iMut];
                    if (Yi <= SMALL)
                        continue;

                    const double Mw = mutationMixPtr_->speciesMw(iMut); // kg/mol
                    if (Mw > SMALL && finite(Mw))
                    {
                        Rmix += Yi * Ru / Mw; // J/kg/K
                    }
                }

                if (Rmix <= SMALL || !finite(Rmix))
                {
                    Tve[celli] = T[celli];
                    continue;
                }

                // --- Energy bookkeeping (Et/Ev are volumetric J/m3 in your wrapper)
                const double Etot_new = rho[celli] * e[celli];

                double Et_prev = Et_[celli];
                double Ev_prev = Ev_[celli];

                // If Et/Ev not initialized sensibly, rebuild Et from T
                // (avoids the Et_prev=1 -> huge scale behavior)
                if (
                    !finite(Et_prev) || !finite(Ev_prev) ||
                    (Et_prev + Ev_prev) <= 0.1 * std::max(Etot_new, 0.0) || // heuristic
                    (Et_prev + Ev_prev) <= SMALL)
                {
                    const double Ttr0 = std::max(double(T[celli]), 300.0);
                    Et_prev = rho[celli] * 1.5 * Rmix * Ttr0;
                    Ev_prev = 0.0;
                }

                double scale = 1.0;
                if ((Et_prev + Ev_prev) > SMALL && finite(Etot_new) && Etot_new > SMALL)
                {
                    scale = Etot_new / (Et_prev + Ev_prev);
                }
                if (!finite(scale) || scale <= 0.0)
                    scale = 1.0;

                // --- Temperatures
                double Ttr = double(T[celli]);
                double Tv = double(Tve[celli]);

                double Ev_input = mutationMixPtr_->EvFromTv(Tv, rho[celli], Y_mut);
                double Et_input = rho[celli] * e[celli] - Ev_input;

                if (Et_input < 0)
                {
                    Et_input = rho[celli] * e[celli];
                    Ev_input = 0;
                }

                // Clamp (prevents MW correlation singularities)
                const double Tmin = 50.0;
                Ttr = std::max(Ttr, Tmin);
                Tv = std::max(Tv, Tmin);

                // --- Final hard check before calling Mutation++
                if (
                    !finite(Ttr) || !finite(Tv) ||
                    !finite(Et_input) || !finite(Ev_input))
                {
                    FatalErrorInFunction
                        << "Non-finite state before Mutation++ step at cell " << celli << nl
                        << "rho=" << rho[celli] << " Rmix=" << Rmix << nl
                        << "Ttr=" << Ttr << " Tv=" << Tv << nl
                        << "Et=" << Et_input << " Ev=" << Ev_input << nl
                        << exit(FatalError);
                }

                Info << "DBG cell " << celli << nl << "rho=" << rho[celli] << " finite=" << finite(rho[celli]) << nl << "Ttr=" << Ttr << " finite=" << finite(Ttr) << nl << "Tv=" << Tv << " finite=" << finite(Tv) << nl << "Rmix=" << Rmix << " finite=" << finite(Rmix) << nl << "Etot_new=" << Etot_new << " finite=" << finite(Etot_new) << nl << "Et_prev=" << Et_prev << " finite=" << finite(Et_prev) << nl << "Ev_prev=" << Ev_prev << " finite=" << finite(Ev_prev) << nl << "scale=" << scale << " finite=" << finite(scale) << nl << "Et_input=" << Et_input << " finite=" << finite(Et_input) << nl << "Ev_input=" << Ev_input << " finite=" << finite(Ev_input) << nl << endl;
                const scalar dtSolver = this->mesh().time().deltaTValue();
                Info << "Using solver dt = " << dtSolver << endl;
                // --- Call Mutation++ with properly ordered composition
                mutationMixPtr_->step(
                    dtSolver,
                    rho[celli],
                    Y_mut,
                    Et_input,
                    Ev_input,
                    Ttr,
                    Tv);

                // --- Write back
                T[celli] = Ttr;
                Tve[celli] = Tv;
                p[celli] = rho[celli] * Rmix * Ttr;
                e[celli] = (Et_input + Ev_input) / rho[celli];

                Et_[celli] = Et_input;
                Ev_[celli] = Ev_input;
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