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

    if (!heThermoPtr_)
    {
        // ------------------------------------------------------------------------
        // 3. VIBRATIONAL-ELECTRONIC ENERGY (Using local mixture)
        // ------------------------------------------------------------------------

        // 1. Get references to the fields we need for the State
        const volScalarField &p = thermo_.p();
        const volScalarField &T = thermo_.T();

        // FIX 1: Use 'heThermo_->Tve()' instead of 'thermo_.Tve()'
        // 'thermo_' is the standard OpenFOAM class (doesn't know about Tve).
        // 'heThermo_' is your custom wrapper (knows about Tve).

        volScalarField &Tve = heThermoPtr_->Tve();

        mutationMixture &mix = *mutationMixPtr_;

        // 2. Initialize the fields we want to calculate
        volScalarField eve(
            IOobject("eve", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("0", dimEnergy / dimMass, 0.0));

        volScalarField Qve(
            IOobject("Qve", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("0", dimEnergy / dimMass / dimTime, 0.0));

        // 3. Reusable array
        scalarField Y_cell(Y.size());

        // 4. MAIN LOOP
        forAll(eve, celli)
        {
            forAll(Y, k)
            {
                Y_cell[k] = Y[k][celli];
            }

            mix.setState(
                p[celli],
                T[celli],
                Tve[celli],
                Y_cell);

            scalar cellEve = 0.0;
            forAll(Y, k)
            {
                cellEve += Y[k][celli] * mix.speciesEve(k);
            }
            eve[celli] = cellEve;

            // Get Source Term [J/m^3/s]
            Qve[celli] = mix.getVibrationalSource();
        }

        eve.correctBoundaryConditions();
        Qve.correctBoundaryConditions();

        // ------------------------------------------------------------------------
        // 5. SOLVE THE EQUATION
        // ------------------------------------------------------------------------

        fvScalarMatrix EveEqn(
            fvm::ddt(rho, eve) + fvc::div(phi, eve) - fvm::laplacian(thermophysicalTransport->alpha(), eve) ==
            Qve);

        EveEqn.relax();
        EveEqn.solve();

        // --- 4. UPDATE TEMPERATURES (Inversion Step) ---

        forAll(T, celli)
        {
            for (label k = 0; k < Y.size(); ++k)
                Y_cell[k] = Y[k][celli];

            scalar e_total_new = thermo_.he()[celli];
            scalar eve_new = eve[celli];

            // A. Update Tve
            scalar Tve_old = Tve[celli];
            Tve[celli] = mix.solveTve(rho[celli], eve_new, Y_cell, Tve_old);

            // B. Update T (Translational)
            // e_trans_rot = e_total - e_vib
            scalar e_tr_target = e_total_new - eve_new;
            scalar T_old = thermo_.T()[celli];

            thermo_.T()[celli] = mix.solveT(rho[celli], e_tr_target, Y_cell, T_old);
        }

        Tve.correctBoundaryConditions();
        thermo_.T().correctBoundaryConditions();
    }

    // FIX 3: Use 'heThermo_->correct()' to ensure T and Tve are updated together
    thermo_.correct();
}

// ************************************************************************* //
