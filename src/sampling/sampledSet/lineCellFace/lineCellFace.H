/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::sampledSets::lineCellFace

Description
    Face-intersections along a line, plus cell-samples at the mid-points
    in-between

Usage
    \table
        Property    | Description                            | Req'd? | Default
        start       | The start point of the line            | yes    |
        end         | The end point of the line              | yes    |
        axis        | The coordinate axis that is written    | yes    |
    \endtable

    Example specification:
    \verbatim
    {
        type        lineCellFace;
        start       (0.55 0.6 0.5);
        end         (0.55 -0.3 -0.1);
        axis        x;
    }
    \endverbatim

SourceFiles
    lineCellFace.C

\*---------------------------------------------------------------------------*/

#ifndef lineCellFace_H
#define lineCellFace_H

#include "lineCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{

/*---------------------------------------------------------------------------*\
                        Class lineCellFace Declaration
\*---------------------------------------------------------------------------*/

class lineCellFace
:
    public sampledSet
{
    // Private Data

        //- Start point
        const point start_;

        //- End point
        const point end_;

        //- If this line intersects a processor boundary, should both the owner
        //  and neighbour patch faces be stored? In general, for
        //  post-processing no. For some specific programmatic use cases, yes.
        const bool storeBothProcFaces_;


    // Private Member Functions

        //- Calculate all the sampling points
        virtual bool calcSamples
        (
            DynamicList<point>& samplingPositions,
            DynamicList<scalar>& samplingDistances,
            DynamicList<label>& samplingSegments,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces
        ) const;


public:

    //- Runtime type information
    TypeName("lineCellFace");


    // Constructors

        //- Construct from dictionary
        lineCellFace
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from components
        lineCellFace
        (
            const word& name,
            const polyMesh& mesh,
            const word& axis,
            const point& start,
            const point& end
        );


    //- Destructor
    virtual ~lineCellFace();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
