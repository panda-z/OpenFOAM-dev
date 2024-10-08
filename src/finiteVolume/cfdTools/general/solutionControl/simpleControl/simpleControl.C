/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "simpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::simpleControl::read()
{
    return
        fluidSolutionControl::read()
     && singleRegionConvergenceControl::read();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleControl::simpleControl(fvMesh& mesh, const word& algorithmName)
:
    fluidSolutionControl(mesh, algorithmName),
    singleRegionConvergenceControl
    (
        static_cast<singleRegionSolutionControl&>(*this)
    )
{
    printResidualControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleControl::~simpleControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::simpleControl::run(Time& time)
{
    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.run();
}


bool Foam::simpleControl::loop(Time& time)
{
    if (!endIfConverged(time))
    {
        storePrevIterFields();
    }

    return time.loop();
}


// ************************************************************************* //
