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

Class
    Foam::solvers::shockThermo

Description
    Solver module for density-based solution of compressible flow

    Solver module for EXPLICIT density-based solution of compressible REACTIVE
    flows with specie transport.
   
SourceFiles
    shockThermo.C

Author:
    Federico Piscaglia
    Dept. of Aerospace Science and Technology (DAER)
    Politecnico di Milano
    

\*---------------------------------------------------------------------------*/

#ifndef shockThermo_H
#define shockThermo_H

#include "shockFluid.H"
#include "fluidMulticomponentThermo.H"
#include "combustionModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "multivariateScheme.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class shockThermo Declaration
\*---------------------------------------------------------------------------*/

class shockThermo
:
    public shockFluid
{

private:

    // Thermophysical properties

        fluidMulticomponentThermo& thermo_;

    // Composition

        PtrList<volScalarField>& Y_;

    // Reactions

        autoPtr<combustionModel> reaction;
        
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    // Thermophysical transport

        autoPtr<fluidMulticomponentThermophysicalTransportModel>
            thermophysicalTransport;

    // Private Member Functions

        //- Interpolate field vf according to direction dir
        template<class Type>
        inline tmp<SurfaceField<Type>> interpolate
        (
            const VolField<Type>& vf,
            const surfaceScalarField& dir,
            const word& reconFieldName = word::null
        )
        {
            tmp<SurfaceField<Type>> tsf
            (
                fvc::interpolate
                (
                    vf,
                    dir,
                    "reconstruct("
                  + (reconFieldName != word::null ? reconFieldName : vf.name())
                  + ')'
                )
            );

            SurfaceField<Type>& sf = tsf.ref();

            sf.rename(vf.name() + '_' + dir.name());

            return tsf;
        }

        //- Set rDeltaT for LTS
        void setRDeltaT(const surfaceScalarField& amaxSf);
 
public:

    // Public Data

        //- Reference to the fluid thermophysical properties
        const fluidMulticomponentThermo& thermo;

        //- Reference to the composition
        const PtrList<volScalarField>& Y;

 
    //- Runtime type information
        TypeName("shockThermo");

    // Constructors

        //- Construct from region mesh
        shockThermo(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        shockThermo(const shockThermo&) = delete;

    //- Destructor
        virtual ~shockThermo();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        // //- Called at the start of the PIMPLE loop to move the mesh
        // virtual void moveMesh();

        // //- Called at the start of the PIMPLE loop
        // virtual void prePredictor();

        // //- Construct and optionally solve the momentum equation
        // virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        // //- Construct and solve the pressure equation in the PISO loop
        // virtual void pressureCorrector();

        // //- Correct the momentum and thermophysical transport modelling
        // virtual void postCorrector();

        // //- Called after the PIMPLE loop at the end of the time-step
        // virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const shockThermo&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
