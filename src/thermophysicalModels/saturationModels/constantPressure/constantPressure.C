/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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

#include "saturationModels.H"
#include "constantPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(constantPressure, 0);
    addToRunTimeSelectionTable
    (
        saturationPressureModel,
        constantPressure,
        dictionary
    );

    static const dimensionedScalar oneP(dimPressure, 1);
    static const dimensionedScalar zeroPbyT(dimPressure/dimTemperature, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
Foam::tmp<FieldType>
Foam::saturationModels::constantPressure::pSat(const FieldType& T) const
{
    return evaluate(T, "pSat", pSat_);
}


template<class FieldType>
Foam::tmp<FieldType>
Foam::saturationModels::constantPressure::pSatPrime(const FieldType& T) const
{
    return evaluate(T, "pSatPrime", zeroPbyT);
}


template<class FieldType>
Foam::tmp<FieldType>
Foam::saturationModels::constantPressure::lnPSat(const FieldType& T) const
{
    return evaluate(T, "lnPSat", log(pSat_/oneP));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::constantPressure::constantPressure
(
    const dictionary& dict
)
:
    saturationPressureModel(),
    pSat_("value", dimPressure, dict)
{}


Foam::saturationModels::constantPressure::constantPressure
(
    const dimensionedScalar& pSat
)
:
    saturationPressureModel(),
    pSat_(pSat)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::constantPressure::~constantPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

IMPLEMENT_PSAT(saturationModels::constantPressure, scalarField);


IMPLEMENT_PSAT(saturationModels::constantPressure, volScalarField::Internal);


IMPLEMENT_PSAT(saturationModels::constantPressure, volScalarField);


// ************************************************************************* //
