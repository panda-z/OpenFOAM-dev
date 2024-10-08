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

#include "manual.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(manual, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        manual,
        decomposer
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        manual,
        distributor
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::manual::manual
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict),
    decompDataFile_(methodDict.lookup("dataFile"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::manual::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    labelIOList finalDecomp
    (
        IOobject
        (
            decompDataFile_,
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );

    // check if the final decomposition is OK

    if (finalDecomp.size() != points.size())
    {
        FatalErrorInFunction
            << "Size of decomposition list does not correspond "
            << "to the number of points.  Size: "
            << finalDecomp.size() << " Number of points: "
            << points.size()
            << ".\n" << "Manual decomposition data read from file "
            << decompDataFile_ << "." << endl
            << exit(FatalError);
    }

    if (min(finalDecomp) < 0 || max(finalDecomp) > nProcessors_ - 1)
    {
        FatalErrorInFunction
            << "According to the decomposition, cells assigned to "
            << "impossible processor numbers.  Min processor = "
            << min(finalDecomp) << " Max processor = " << max(finalDecomp)
            << ".\n" << "Manual decomposition data read from file "
            << decompDataFile_ << "." << endl
            << exit(FatalError);
    }

    return finalDecomp;
}


// ************************************************************************* //
