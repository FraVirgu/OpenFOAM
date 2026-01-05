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

            // Use the actual solver time-step for coupling
            const scalar dtSolver = this->mesh().time().deltaTValue();

            static constexpr double Ru = 8.31446261815324; // J/mol/K
            const double Tmin = 50.0;

            forAll(T, celli)
            {
                if (!mutationMixPtr_.valid())
                {
                    // Fallback: 1-T
                    Tve[celli] = T[celli];
                    continue;
                }

                if (rho[celli] <= SMALL || !finite(rho[celli]))
                {
                    Tve[celli] = T[celli];
                    continue;
                }

                // ------------------------------------------------------------
                // 1) Build Y in Mutation++ ordering (name-mapped via ofToMut_)
                // ------------------------------------------------------------
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

                // Clip + renormalize
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

                // ------------------------------------------------------------
                // 2) Compute Rmix (Mutation MWs, Mutation ordering)
                // ------------------------------------------------------------
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

                // ------------------------------------------------------------
                // 3) Total energy density from OpenFOAM (authoritative)
                // ------------------------------------------------------------
                const double Etot = rho[celli] * e[celli];

                if (!finite(Etot) || Etot <= SMALL)
                {
                    Tve[celli] = T[celli];
                    continue;
                }

                // ------------------------------------------------------------
                // 4) Use stored Et_/Ev_ as the evolving state (DO NOT rebuild every step)
                // ------------------------------------------------------------
                double Et_local = Et_[celli];
                double Ev_local = Ev_[celli];

                // Initialize if missing/invalid
                if (!finite(Et_local) || !finite(Ev_local) || (Et_local + Ev_local) <= SMALL)
                {
                    double Tv0 = std::max(double(Tve[celli]), 300.0);

                    Ev_local = mutationMixPtr_->EvFromTv(Tv0, rho[celli], Y_mut);
                    if (!finite(Ev_local) || Ev_local < 0.0)
                        Ev_local = 0.0;

                    // Keep Ev below total energy
                    Ev_local = std::min(Ev_local, 0.99 * Etot);
                    Et_local = Etot - Ev_local;
                }

                // Positivity
                Et_local = std::max(Et_local, 0.0);
                Ev_local = std::max(Ev_local, 0.0);

                // Enforce conservation to match OpenFOAM energy
                const double sumE0 = Et_local + Ev_local;
                if (sumE0 > SMALL && finite(sumE0))
                {
                    const double fac = Etot / sumE0;
                    Et_local *= fac;
                    Ev_local *= fac;
                }
                else
                {
                    // fallback: all translational
                    Et_local = Etot;
                    Ev_local = 0.0;
                }

                // ------------------------------------------------------------
                // 5) Temperatures passed to Mutation++ (caller-provided)
                // ------------------------------------------------------------
                double Ttr = std::max(double(T[celli]), Tmin);
                double Tv = std::max(double(Tve[celli]), Tmin);

                if (!finite(Ttr) || !finite(Tv))
                {
                    Tve[celli] = T[celli];
                    continue;
                }

                // ------------------------------------------------------------
                // 6) Advance VT exchange using Mutation++ (updates Et_local/Ev_local)
                // ------------------------------------------------------------
                mutationMixPtr_->step(dtSolver, rho[celli], Y_mut, Et_local, Ev_local, Ttr, Tv);

                // After update: enforce positivity and conservation again
                Et_local = std::max(Et_local, 0.0);
                Ev_local = std::max(Ev_local, 0.0);

                const double sumE1 = Et_local + Ev_local;
                if (sumE1 > SMALL && finite(sumE1))
                {
                    const double fac = Etot / sumE1;
                    Et_local *= fac;
                    Ev_local *= fac;
                }
                else
                {
                    Et_local = Etot;
                    Ev_local = 0.0;
                }

                // ------------------------------------------------------------
                // 7) Write back fields
                // ------------------------------------------------------------
                Et_[celli] = Et_local;
                Ev_[celli] = Ev_local;

                // Keep e consistent with Et+Ev (still total internal energy / rho)
                e[celli] = (Et_local + Ev_local) / rho[celli];

                // OPTIONAL (recommended): update Tve from Ev so it evolves and you can see it
                // If your wrapper does not expose invertTv publicly, add a public TvFromEv wrapper.
                // --- Update Tv from Ev (use previous Tve as initial guess)
                double Tv_new = mutationMixPtr_->invertTv(
                    Ev_local,
                    rho[celli],
                    Y_mut,
                    double(Tve[celli]) // Tv_init
                );

                if (finite(Tv_new) && Tv_new > Tmin)
                {
                    Tve[celli] = Tv_new;
                }

                // --- Update Ttr from Et (use updated Tv, previous T as initial guess)
                double Ttr_new = mutationMixPtr_->invertTtr(
                    Et_local,
                    rho[celli],
                    Y_mut,
                    double(Tve[celli]), // Tv_fixed
                    double(T[celli])    // Ttr_init
                );

                if (finite(Ttr_new) && Ttr_new > Tmin)
                {
                    T[celli] = Ttr_new;
                }

                // --- Now compute pressure using UPDATED T
                p[celli] = rho[celli] * Rmix * T[celli];

                // Debug one cell if you want
                if (celli == 0)
                {
                    Info << "DBG cell0: Et=" << Et_local << " Ev=" << Ev_local
                         << " Etot=" << Etot << " T=" << T[celli] << " Tve=" << Tve[celli]
                         << " dt=" << dtSolver << nl;
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