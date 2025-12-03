/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "shockFluidThermo.H"
#include "fvMeshStitcher.H"
#include "localEulerDdtScheme.H"
#include "hydrostaticInitialisation.H"
#include "fvcMeshPhi.H"
#include "fvcVolumeIntegrate.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace solvers
    {
        defineTypeNameAndDebug(shockFluidThermo, 0);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::shockFluidThermo::shockFluidThermo(fvMesh &mesh, autoPtr<fluidThermo> thermoPtr)
    : shockFluid(mesh),

      thermoPtr_(fluidThermo::New(mesh)),

      thermo_(thermoPtr_()),

      p_(thermo_.p()),

      rho_(
          IOobject(
              "rho",
              runTime.name(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          thermo_.renameRho()),

      U_(
          IOobject(
              "U",
              runTime.name(),
              mesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
          mesh),

      phi_(
          IOobject(
              "phi",
              runTime.name(),
              mesh,
              IOobject::READ_IF_PRESENT,
              IOobject::AUTO_WRITE),
          linearInterpolate(rho_ * U_) & mesh.Sf()),

      K("K", 0.5 * magSqr(U_)),

      inviscid(
          max(thermo_.mu().primitiveField()) > 0
              ? false
              : true),

      momentumTransport(
          inviscid
              ? autoPtr<compressibleMomentumTransportModel>(nullptr)
              : compressible::momentumTransportModel::New(
                    rho_,
                    U_,
                    phi_,
                    thermo_)),

      thermophysicalTransport(
          inviscid
              ? autoPtr<fluidThermoThermophysicalTransportModel>(nullptr)
              : fluidThermoThermophysicalTransportModel::New(
                    momentumTransport(),
                    thermo_)),

      fluxScheme(
          mesh.schemes().dict().lookupOrDefault<word>("fluxScheme", "Kurganov")),

      thermo(thermo_),
      p(p_),
      rho(rho_),
      U(U_),
      phi(phi_)
{
    thermo.validate(type(), "e");

    if (momentumTransport.valid())
    {
        momentumTransport->validate();
        mesh.schemes().setFluxRequired(U.name());
    }

    fluxPredictor();

    if (transient())
    {
        const surfaceScalarField amaxSf(
            max(mag(aphiv_pos()), mag(aphiv_neg())));

        correctCoNum(amaxSf);
    }
    else if (LTS)
    {
        Info << "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>(
            new volScalarField(
                IOobject(
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE),
                mesh,
                dimensionedScalar(dimless / dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName));
    }
}

Foam::solvers::shockFluidThermo::~shockFluidThermo() = default;

// ************************************************************************* //
