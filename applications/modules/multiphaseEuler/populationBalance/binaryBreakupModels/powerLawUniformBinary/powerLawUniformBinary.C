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

#include "powerLawUniformBinary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace binaryBreakupModels
{
    defineTypeNameAndDebug(powerLawUniformBinary, 0);
    addToRunTimeSelectionTable
    (
        binaryBreakupModel,
        powerLawUniformBinary,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::binaryBreakupModels::powerLawUniformBinary::
powerLawUniformBinary
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binaryBreakupModel(popBal, dict),
    power_(dict.lookup<scalar>("power"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::binaryBreakupModels::powerLawUniformBinary::
addToBinaryBreakupRate
(
    volScalarField::Internal& binaryBreakupRate,
    const label i,
    const label j
)
{
    const dimensionedScalar& vj = popBal_.vs()[j];

    binaryBreakupRate.primitiveFieldRef() +=
        pow(vj.value(), power_)*2/vj.value();
}


// ************************************************************************* //
