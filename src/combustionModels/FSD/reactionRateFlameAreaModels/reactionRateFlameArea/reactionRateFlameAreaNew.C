/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "reactionRateFlameArea.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::reactionRateFlameArea> Foam::reactionRateFlameArea::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
{
    const word reactionRateFlameAreaType
    (
        dict.lookup("reactionRateFlameArea")
    );

    Info<< "Selecting reaction rate flame area correlation "
        << reactionRateFlameAreaType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reactionRateFlameAreaType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Unknown reactionRateFlameArea type "
            << reactionRateFlameAreaType << endl << endl
            << "Valid reaction rate flame area types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    const label tempOpen = reactionRateFlameAreaType.find('<');

    const word modelType = reactionRateFlameAreaType(0, tempOpen);

    Info<< incrIndent;

    autoPtr<reactionRateFlameArea> modelPtr
    (
        cstrIter()
        (
            modelType,
            dict,
            dict.optionalSubDict(modelType + "Coeffs"),
            mesh,
            combModel
        )
    );

    Info<< decrIndent;

    return modelPtr;
}


// ************************************************************************* //
