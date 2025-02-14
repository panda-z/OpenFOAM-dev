/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "function1Temperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(function1Temperature, 0);
    addToRunTimeSelectionTable
    (
        saturationTemperatureModel,
        function1Temperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::function1Temperature::function1Temperature
(
    const dictionary& dict
)
:
    saturationTemperatureModel(),
    function_
    (
        Function1<scalar>::New("function", dimPressure, dimTemperature, dict)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::function1Temperature::~function1Temperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationModels::function1Temperature::Tsat
(
    const volScalarField::Internal& p
) const
{
    tmp<volScalarField::Internal> tTsat
    (
        volScalarField::Internal::New
        (
            "Tsat",
            p.mesh(),
            dimensionedScalar(dimTemperature, 0)
        )
    );

    volScalarField::Internal& Tsat = tTsat.ref();

    Tsat.primitiveFieldRef() = function_->value(p.primitiveField());

    return tTsat;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationModels::function1Temperature::TsatPrime
(
    const volScalarField::Internal& p
) const
{
    const scalar dp(rootSmall*p.average().value());

    tmp<volScalarField::Internal> tTsatPrime
    (
        volScalarField::Internal::New
        (
            "TsatPrime",
            p.mesh(),
            dimensionedScalar(dimTemperature/dimPressure, 0)
        )
    );

    volScalarField::Internal& TsatPrime = tTsatPrime.ref();

    TsatPrime.primitiveFieldRef() =
        (
            function_->value(p.primitiveField() + dp/2)
          - function_->value(p.primitiveField() - dp/2)
        )/dp;

    return tTsatPrime;
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::function1Temperature::Tsat
(
    const volScalarField& p
) const
{
    tmp<volScalarField> tTsat
    (
        volScalarField::New
        (
            "Tsat",
            p.mesh(),
            dimensionedScalar(dimTemperature, 0)
        )
    );

    volScalarField& Tsat = tTsat.ref();

    Tsat.primitiveFieldRef() = function_->value(p.primitiveField());

    volScalarField::Boundary& TsatBf = Tsat.boundaryFieldRef();

    forAll(Tsat.boundaryField(), patchi)
    {
        const scalarField& pp = p.boundaryField()[patchi];

        TsatBf[patchi] = function_->value(pp);
    }

    return tTsat;
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::function1Temperature::TsatPrime
(
    const volScalarField& p
) const
{
    const scalar dp(rootSmall*p.average().value());

    tmp<volScalarField> tTsatPrime
    (
        volScalarField::New
        (
            "TsatPrime",
            p.mesh(),
            dimensionedScalar(dimTemperature/dimPressure, 0)
        )
    );

    volScalarField& TsatPrime = tTsatPrime.ref();

    TsatPrime.primitiveFieldRef() =
        (
            function_->value(p.primitiveField() + dp/2)
          - function_->value(p.primitiveField() - dp/2)
        )/dp;

    volScalarField::Boundary& TsatPrimeBf = TsatPrime.boundaryFieldRef();

    forAll(TsatPrime.boundaryField(), patchi)
    {
        const scalarField& pp = p.boundaryField()[patchi];

        TsatPrimeBf[patchi] =
            (
                function_->value(pp + dp/2)
              - function_->value(pp - dp/2)
            )/dp;
    }

    return tTsatPrime;
}


// ************************************************************************* //
