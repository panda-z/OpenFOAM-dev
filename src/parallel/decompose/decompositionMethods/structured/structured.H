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

Class
    Foam::decompositionMethods::structured

Description
    Decomposition by walking out decomposition of patch cells mesh.

SourceFiles
    structured.C

\*---------------------------------------------------------------------------*/

#ifndef structured_H
#define structured_H

#include "decompositionMethod.H"

namespace Foam
{
namespace decompositionMethods
{

/*---------------------------------------------------------------------------*\
                      Class structured Declaration
\*---------------------------------------------------------------------------*/

class structured
:
    public decompositionMethod
{
    // Private Data

        wordReList patches_;

        autoPtr<decompositionMethod> method_;


public:

    //- Runtime type information
    TypeName("structured");


    // Constructors

        //- Construct given the decomposition dictionary
        structured
        (
            const dictionary& decompositionDict,
            const dictionary& methodDict
        );

        //- Disallow default bitwise copy construction
        structured(const structured&) = delete;


    //- Destructor
    virtual ~structured()
    {}


    // Member Functions

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const structured&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
