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
    Foam::trackedParticle

Description
    Particle class that marks cells it passes through. Used to mark cells
    visited by feature edges.

SourceFiles
    trackedParticle.C

\*---------------------------------------------------------------------------*/

#ifndef trackedParticle_H
#define trackedParticle_H

#include "particle.H"
#include "Cloud.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class trackedParticle;

Ostream& operator<<(Ostream&, const trackedParticle&);


/*---------------------------------------------------------------------------*\
                       Class trackedParticle Declaration
\*---------------------------------------------------------------------------*/

class trackedParticle
:
    public particle
{
    // Private Data

        //- Start point to track from
        point start_;

        //- End point to track to
        point end_;

        //- Level of this particle
        label level_;

        //- Passive label (used to store feature edge mesh)
        label i_;

        //- Passive label (used to store feature edge point)
        label j_;

        //- Passive label (used to store feature edge label)
        label k_;


public:

    friend class lagrangian::Cloud<trackedParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
    public:

        const scalar maxTrackLen_;

        labelList& maxLevel_;

        List<PackedBoolList>& featureEdgeVisited_;


        // Constructors

            trackingData
            (
                lagrangian::Cloud<trackedParticle>& cloud,
                const scalar maxTrackLen,
                labelList& maxLevel,
                List<PackedBoolList>& featureEdgeVisited
            )
            :
                particle::trackingData(cloud),
                maxTrackLen_(maxTrackLen),
                maxLevel_(maxLevel),
                featureEdgeVisited_(featureEdgeVisited)
            {}
    };



    // Constructors

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology
        trackedParticle
        (
            const meshSearch& searchEngine,
            const vector& position,
            const label celli,
            label& nLocateBoundaryHits,
            const point& end,
            const label level,
            const label i,
            const label j,
            const label k
        );

        //- Construct from Istream
        trackedParticle(Istream& is, bool readFields = true);

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new trackedParticle(*this));
        }

        //- Construct from Istream and return
        static autoPtr<trackedParticle> New(Istream& is)
        {
            return autoPtr<trackedParticle>(new trackedParticle(is));
        }


    // Member Functions

        //- Point to track from
        point& start()
        {
            return start_;
        }

        //- Point to track to
        point& end()
        {
            return end_;
        }

        //- Transported label
        label i() const
        {
            return i_;
        }

        //- Transported label
        label& i()
        {
            return i_;
        }

        //- Transported label
        label j() const
        {
            return j_;
        }

        //- Transported label
        label& j()
        {
            return j_;
        }

        //- Transported label
        label k() const
        {
            return k_;
        }

        //- Transported label
        label& k()
        {
            return k_;
        }



        // Tracking

            //- Track all particles to their end point
            bool move(lagrangian::Cloud<trackedParticle>&, trackingData&);

            //- Overridable function to handle the particle hitting a wedge
            void hitWedgePatch
            (
                lagrangian::Cloud<trackedParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            void hitSymmetryPlanePatch
            (
                lagrangian::Cloud<trackedParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            void hitSymmetryPatch
            (
                lagrangian::Cloud<trackedParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a cyclic
            void hitCyclicPatch
            (
                lagrangian::Cloud<trackedParticle>&,
                trackingData&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                lagrangian::Cloud<trackedParticle>&,
                trackingData&
            );

            //- Do corrections to the particle and tracking data following a
            //  transfer between processors
            void correctAfterParallelTransfer
            (
                lagrangian::Cloud<trackedParticle>&,
                trackingData&
            );


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const trackedParticle&);
};


template<>
inline bool contiguous<trackedParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
