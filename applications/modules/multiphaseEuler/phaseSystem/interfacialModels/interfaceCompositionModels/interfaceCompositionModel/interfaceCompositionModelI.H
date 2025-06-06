/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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

#include "interfaceCompositionModel.H"
#include "rhoFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::sidedPhaseInterface&
Foam::interfaceCompositionModel::interface() const
{
    return interface_;
}


const Foam::hashedWordList& Foam::interfaceCompositionModel::species() const
{
    return species_;
}


const Foam::rhoFluidMulticomponentThermo&
Foam::interfaceCompositionModel::thermo() const
{
    return thermo_;
}


const Foam::rhoThermo& Foam::interfaceCompositionModel::otherThermo() const
{
    return otherThermo_;
}


const Foam::rhoFluidMulticomponentThermo&
Foam::interfaceCompositionModel::otherMulticomponentThermo() const
{
    return refCast<const rhoFluidMulticomponentThermo>(otherThermo_);
}


// ************************************************************************* //
