/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::pointEdgePoint

Description
    Holds information regarding nearest wall point. Used in PointEdgeWave.
    (so not standard FaceCellWave)
    To be used in wall distance calculation.

SourceFiles
    pointEdgePointI.H
    pointEdgePoint.C

\*---------------------------------------------------------------------------*/

#ifndef pointEdgePoint_H
#define pointEdgePoint_H

#include "point.H"
#include "label.H"
#include "scalar.H"
#include "tensor.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;


// Forward declaration of friend functions and operators

class pointEdgePoint;

Istream& operator>>(Istream&, pointEdgePoint&);
Ostream& operator<<(Ostream&, const pointEdgePoint&);


/*---------------------------------------------------------------------------*\
                       Class pointEdgePoint Declaration
\*---------------------------------------------------------------------------*/

class pointEdgePoint
{
    // Private Data

        //- Position of nearest wall center
        point origin_;

        //- Normal distance (squared) from point to origin
        scalar distSqr_;


    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer pt. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const pointEdgePoint& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Constructors

        //- Construct null
        inline pointEdgePoint();

        //- Construct from origin, distance
        inline pointEdgePoint(const point&, const scalar);


    // Member Functions

        // Access

            inline const point& origin() const;

            inline scalar distSqr() const;


        // Needed by PointEdgeWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of edge on point
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const label edgeI,
                const pointEdgePoint& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const pointEdgePoint& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge.
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgeI,
                const label pointi,
                const pointEdgePoint& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const pointEdgePoint&, TrackingData& td) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const pointEdgePoint&) const;
        inline bool operator!=(const pointEdgePoint&) const;


    // IOstream Operators

        inline friend Ostream& operator<<(Ostream&, const pointEdgePoint&);
        inline friend Istream& operator>>(Istream&, pointEdgePoint&);
};


//- Data associated with pointEdgePoint type are contiguous
template<>
inline bool contiguous<pointEdgePoint>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointEdgePointI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
