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
    Foam::sampledSets::points

Description
    Specified point samples. Optionally ordered into a continuous path.
    Ordering is an optimisation; it enables tracking from one point to the
    next. If ordering is off, each point is searched for individually.

Usage
    \table
        Property    | Description                            | Req'd? | Default
        points      | The points to sample                   | yes    |
        ordered     | Are the points in order?               | yes    |
        axis        | The coordinate axis that is written    | yes    |
    \endtable

    Example specification:
    \verbatim
    {
        type        points;
        points
        (
            (0.95 0 0.25)
            (0.951251 0 0.250119)
            (0.952468 0 0.250473)
            (0.953618 0 0.251057)
            (0.954669 0 0.251859)
            (0.95559 0 0.252865)
            (0.956353 0 0.254057)
            (0.956931 0 0.255413)
            (0.9573 0 0.256908)
            (0.957441 0 0.258513)
        );
        ordered     yes;
        axis        x;
    }
    \endverbatim

SourceFiles
    points.C

\*---------------------------------------------------------------------------*/

#ifndef points_H
#define points_H

#include "sampledSet.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{

/*---------------------------------------------------------------------------*\
                           Class points Declaration
\*---------------------------------------------------------------------------*/

class points
:
    public sampledSet
{
    // Private Data

        //- Sampling points
        const List<point> points_;

        //- Do the points form an ordered sequence?
        const bool ordered_;


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
    TypeName("points");


    // Static Member Functions

        //- Calculate all the sampling points
        static void calcSamples
        (
            const polyMesh& mesh,
            const pointField& points,
            DynamicList<point>& samplingPositons,
            DynamicList<scalar>& samplingDistances,
            DynamicList<label>& samplingSegments,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces
        );


    // Constructors

        //- Construct from dictionary
        points
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~points();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
