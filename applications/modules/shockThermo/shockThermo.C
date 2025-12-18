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
#include "fvMeshStitcher.H"
#include "localEulerDdtScheme.H"
#include "hydrostaticInitialisation.H"
#include "fvcMeshPhi.H"
#include "fvcVolumeIntegrate.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"
#include "addToRunTimeSelectionTable.H"

#include "psiThermo.H"

#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace solvers
    {
        defineTypeNameAndDebug(shockThermo, 0);
        addToRunTimeSelectionTable(solver, shockThermo, fvMesh);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::shockThermo::shockThermo(fvMesh &mesh)
    : shockFluid(mesh, fluidThermo::New(mesh)),
      thermo_(refCast<fluidMulticomponentThermo>(shockFluid::thermo_)),

      Y_(thermo_.Y()),

      reaction(combustionModel::New(thermo_, momentumTransport())),

      thermophysicalTransport(
          fluidMulticomponentThermophysicalTransportModel::New(
              momentumTransport(),
              thermo_)),
      thermo(thermo_),
      Y(Y_)

{
    thermo.validate(type(), "h", "e");

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::shockThermo::~shockThermo()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockThermo::preSolve()
{
    // cast the thermo class to support multicomponent fluid.
    refCast<fluidMulticomponentThermo>(shockFluid::thermo_);

    {
        const surfaceScalarField amaxSf(
            max(mag(aphiv_pos()), mag(aphiv_neg())));

        if (transient())
        {
            correctCoNum(amaxSf);
        }
        else if (LTS)
        {
            setRDeltaT(amaxSf);
        }
    }

    fvModels().preUpdateMesh();

    if (mesh.topoChanging() || mesh.stitcher().stitches())
    {
        pos.clear();
        neg.clear();

        clearTemporaryFields();
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}

// ************************************************************************* //
