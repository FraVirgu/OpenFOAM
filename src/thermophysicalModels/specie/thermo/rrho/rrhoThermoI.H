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

\*---------------------------------------------------------------------------*/

#include "janafThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::rrhoThermo<EquationOfState>::rrhoThermo
(
    const EquationOfState& st,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon,
    const typename rrhoThermo<EquationOfState>::coeffArray& highCpCoeffs,
    const typename rrhoThermo<EquationOfState>::coeffArray& lowCpCoeffs,
    const bool convertCoeffs
)
:
janafThermo<EquationOfState>(st,Tlow,Thigh,Tcommon,highCpCoeffs,lowCpCoeffs,convertCoeffs)
{

}


template<class EquationOfState>
inline Foam::rrhoThermo<EquationOfState>::rrhoThermo
(
    const word& name,
    const rrhoThermo& jt
)
:
    janafThermo<EquationOfState>(name, jt)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class EquationOfState>
inline Foam::scalar Foam::rrhoThermo<EquationOfState>::Cpve
(
    const scalar p,
    const scalar T
) const
{
    Info<< "This is Cp vibro-electronic" << endl;
    return this->Cp(p,T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //






// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
  

// ************************************************************************* //
