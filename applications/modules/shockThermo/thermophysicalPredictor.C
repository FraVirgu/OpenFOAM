/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "shockThermo.H"
#include "fvmDdt.H"
#include "fvcDiv.H"
#include "fvcDdt.H"
#include <cmath>

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockThermo::thermophysicalPredictor()
{

    // add support to multi-specie chemistry
    tmp<fv::convectionScheme<scalar>> mvConvection(
        fv::convectionScheme<scalar>::New(
            mesh,
            fields,
            phi,
            mesh.schemes().div("div(phi,Yi_h)")));

    reaction->correct();

    forAll(Y, i)
    {
        volScalarField &Yi = Y_[i];

        if (thermo_.solveSpecie(i))
        {
            fvScalarMatrix YiEqn(
                fvm::ddt(rho, Yi) + mvConvection->fvmDiv(phi, Yi) + thermophysicalTransport->divj(Yi) ==
                reaction->R(Yi) + fvModels().source(rho, Yi));

            YiEqn.relax();

            fvConstraints().constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints().constrain(Yi);
        }
        else
        {
            Yi.correctBoundaryConditions();
        }
    }

    thermo_.normaliseY();

    if (thermo_.he().name() != "e")
    {
        FatalErrorInFunction()
            << "sensible energy e is required as primary variable"
            << exit(FatalError);
    }

    //- ------------------------------------------------------------------------

    // solve the equation of energy.

    // IMPORTANT: I cannot use the function call
    // shockFluid::thermophysicalPredictor(), as it causes inconsistencies
    // between thermo classes. The thermo.correct() function must be the one
    // defined in the derived class; shockFluid::thermophysicalPredictor() is
    // pasted here below.

    volScalarField &e = thermo_.he();

    const surfaceScalarField e_pos(interpolate(e, pos, thermo.T().name()));
    const surfaceScalarField e_neg(interpolate(e, neg, thermo.T().name()));

    surfaceScalarField phiEp(
        "phiEp",
        aphiv_pos() * (rho_pos() * (e_pos + 0.5 * magSqr(U_pos())) + p_pos()) + aphiv_neg() * (rho_neg() * (e_neg + 0.5 * magSqr(U_neg())) + p_neg()) + aSf() * (p_pos() - p_neg()));

    // Make flux for pressure-work absolute
    if (mesh.moving())
    {
        phiEp += mesh.phi() * (a_pos() * p_pos() + a_neg() * p_neg());
    }

    //- for high enthalpy flows, e = e_rt + e_ve.

    fvScalarMatrix EEqn(
        fvm::ddt(rho, e) + fvc::div(phiEp) + fvc::ddt(rho, K) ==
        fvModels().source(rho, e));

    if (!inviscid)
    {
        const surfaceScalarField devTauDotU(
            "devTauDotU",
            devTau() & (a_pos() * U_pos() + a_neg() * U_neg()));

        EEqn += thermophysicalTransport->divq(e) + fvc::div(devTauDotU);
    }

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(e);
    // ------------------------------------------------------------------------
    // 3. VIBRATIONAL-ELECTRONIC ENERGY (Using local mixture)
    // ------------------------------------------------------------------------
    if (heThermoPtr_)
    {

        volScalarField &Tve = heThermoPtr_->Tve();

        mutationMixture &mix = *mutationMixPtr_;

        // 1. Get references to the fields we need for the State
        const volScalarField &p = thermo_.p();
        const volScalarField &T = thermo_.T();

        // FIX 1: Use 'heThermo_->Tve()' instead of 'thermo_.Tve()'
        // 'thermo_' is the standard OpenFOAM class (doesn't know about Tve).
        // 'heThermo_' is your custom wrapper (knows about Tve).

        // 2. Initialize the fields we want to calculate
        volScalarField eve(
            IOobject(
                "eve",
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("0", dimEnergy / dimMass, 0.0));

        volScalarField Qve(
            IOobject(
                "Qve",
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("0", dimEnergy / dimVolume / dimTime, 0.0));

        // 3. SET BOUNDARY CONDITIONS FOR eve AND Qve
        forAll(eve.boundaryFieldRef(), patchi)
        {
            const fvPatch &patch = mesh.boundary()[patchi];
            const word &TpatchType = thermo_.T().boundaryField()[patchi].type();

            // --- eve BC
            eve.boundaryFieldRef().set(
                patchi,
                fvPatchField<scalar>::New(
                    TpatchType, // copy BC type from T
                    patch,
                    eve));

            // --- Qve BC (zeroGradient is safest for a source term)
            Qve.boundaryFieldRef().set(
                patchi,
                fvPatchField<scalar>::New(
                    "zeroGradient",
                    patch,
                    Qve));
        }

        volScalarField kByCp(
            IOobject(
                "kByCp",
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("0", dimMass / (dimLength * dimTime), 0.0) // kg/(m s)
        );

        forAll(kByCp.boundaryFieldRef(), patchi)
        {
            const fvPatch &patch = mesh.boundary()[patchi];
            const word pType = patch.type();

            if (pType == "empty")
            {
                kByCp.boundaryFieldRef().set(
                    patchi,
                    fvPatchField<scalar>::New("empty", patch, kByCp));
            }
            else
            {
                kByCp.boundaryFieldRef().set(
                    patchi,
                    fvPatchField<scalar>::New("zeroGradient", patch, kByCp));
            }
        }

        kByCp.correctBoundaryConditions();

        // 3. Reusable array
        scalarField Y_cell(Y.size());

        forAll(eve, celli)
        {
            // --- Read primitive values
            scalar rhoVal = rho[celli];
            scalar TVal = T[celli];
            scalar TveVal = Tve[celli];

            // --- Check finiteness FIRST (max() does not fix NaNs)
            if (!std::isfinite(rhoVal) || !std::isfinite(TVal) || !std::isfinite(TveVal))
            {
                FatalErrorInFunction
                    << "Non-finite state at cell " << celli
                    << " rho=" << rhoVal << " T=" << TVal << " Tve=" << TveVal
                    << exit(FatalError);
            }

            // --- Clamp + renormalize composition for Mutation++ (NO exact zeros)
            const scalar Ymin = scalar(1e-30); // tiny, but nonzero
            scalar sumY = 0.0;
            scalar minY = GREAT;

            forAll(Y, k)
            {
                scalar yk = Y[k][celli];
                if (!std::isfinite(yk))
                    yk = 0.0;

                minY = min(minY, yk);

                yk = max(yk, Ymin); // <-- key change: enforce strictly positive
                Y_cell[k] = yk;
                sumY += yk;
            }

            // renormalize
            forAll(Y, k) Y_cell[k] /= sumY;

            // --- Physical clamps (finite already guaranteed)
            scalar rho_safe = max(rhoVal, scalar(1e-6)); // safer than 1e-8 for relaxation models
            scalar T_safe = max(TVal, scalar(200.0));
            scalar Tve_safe = max(TveVal, scalar(200.0));

            // --- Push state into Mutation++
            mix.setState(rho_safe, T_safe, Tve_safe, Y_cell);

            // You MUST have already called mix.setState(...) for this cell
            scalar kappa = mix.kappa(); // [W/m/K]
            scalar Cp = mix.Cp();       // [J/kg/K]

            // Safety clamp
            Cp = max(Cp, scalar(1e-6));

            kByCp[celli] = kappa / Cp; // [kg/(m s)]

            // --- Compute mixture vibrational energy
            scalar cellEve = 0.0;
            forAll(Y, k)
            {
                cellEve += Y_cell[k] * mix.speciesEve(k); // use sanitized Y_cell here too
            }
            eve[celli] = cellEve;

            // --- Source term (this is where you currently crash)
            Qve[celli] = mix.getVibrationalSource();
        }

        eve.correctBoundaryConditions();
        Qve.correctBoundaryConditions();

        // ------------------------------------------------------------------------
        // 5. SOLVE THE EQUATION
        // ------------------------------------------------------------------------

        fvScalarMatrix EveEqn(
            fvm::ddt(rho, eve) + fvc::div(phi, eve) - fvm::laplacian(kByCp, eve) ==
            Qve);

        EveEqn.relax();
        EveEqn.solve();

        // --- 4. UPDATE TEMPERATURES (Inversion Step) ---

        forAll(T, celli)
        {
            // --------------------------------------------------
            // 1. Read conserved quantities
            // --------------------------------------------------
            const scalar e_total = thermo_.he()[celli]; // total specific energy
            const scalar eve_new = eve[celli];          // vibrational energy

            // --------------------------------------------------
            // 2. Build a SAFE composition vector for Mutation++
            // --------------------------------------------------
            const scalar Ymin = 1e-20;
            scalar sumY = 0.0;

            forAll(Y, k)
            {
                scalar yk = Y[k][celli];

                if (!std::isfinite(yk))
                    yk = 0.0;

                yk = max(yk, Ymin);
                Y_cell[k] = yk;
                sumY += yk;
            }

            forAll(Y_cell, k)
                Y_cell[k] /= sumY;

            // --------------------------------------------------
            // 3. Safe primitive values
            // --------------------------------------------------
            const scalar rho_safe =
                max(rho[celli], scalar(1e-6));

            const scalar T_old =
                max(thermo_.T()[celli], scalar(200.0));

            const scalar Tve_old =
                max(Tve[celli], scalar(200.0));

            // --------------------------------------------------
            // 4. IMPORTANT: refresh Mutation++ state
            // --------------------------------------------------
            mix.setState(rho_safe, T_old, Tve_old, Y_cell);

            // --------------------------------------------------
            // 5. Invert vibrational temperature from eve
            // --------------------------------------------------
            Tve[celli] =
                mix.solveTve(
                    rho_safe,
                    eve_new,
                    Y_cell,
                    Tve_old);

            // --------------------------------------------------
            // 6. Invert translational temperature from (e − eve)
            // --------------------------------------------------
            const scalar e_tr_target = e_total - eve_new;

            thermo_.T()[celli] =
                mix.solveT(
                    rho_safe,
                    e_tr_target,
                    Y_cell,
                    T_old);
        }

        Tve.correctBoundaryConditions();
        thermo_.T().correctBoundaryConditions();
        // Make sure the thermo object keeps your updated Tve
        heThermoPtr_->correctTve(Tve);
    }

    thermo_.correct();
}

// ************************************************************************* //
