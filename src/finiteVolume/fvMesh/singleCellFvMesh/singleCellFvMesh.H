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
    Foam::singleCellFvMesh

Description
    fvMesh as subset of other mesh. Consists of one cell and all original
    boundary faces. Useful when manipulating boundary data. Single internal
    cell only needed to be able to manipulate in a standard way.

SourceFiles
    singleCellFvMesh.C
    singleCellFvMeshInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef singleCellFvMesh_H
#define singleCellFvMesh_H

#include "generalFieldMapper.H"
#include "fvMesh.H"
#include "labelListIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class singleCellFvMesh Declaration
\*---------------------------------------------------------------------------*/

class singleCellFvMesh
:
    public fvMesh
{
    // Private Data

        const labelListIOList patchFaceAgglomeration_;

        //- From patch faces back to agglomeration or fine mesh
        labelListIOList patchFaceMap_;

        //- From fine mesh faces to coarse mesh
        labelIOList reverseFaceMap_;

        //- From coarse points back to original mesh
        labelIOList pointMap_;

        //- From fine points to coarse mesh
        labelIOList reversePointMap_;


    // Private Member Functions

        //- Calculate agglomerated mesh
        void agglomerateMesh(const fvMesh&, const labelListList&);


public:

    // Public Classes

        //- Patch field mapper class for agglomerated meshes
        class agglomPatchFieldMapper
        :
            public generalFieldMapper
        {
            // Private Data

                //- Interpolated addressing
                const labelListList& addressing_;

                //- Interpolation weights
                const scalarListList& weights_;


        public:

            // Constructors

                //- Construct given addressing
                agglomPatchFieldMapper
                (
                    const labelListList& addressing,
                    const scalarListList& weights
                )
                :
                    addressing_(addressing),
                    weights_(weights)
                {}


            // Member Functions

                //- Is the mapping direct?
                virtual bool direct() const
                {
                    return false;
                }

                //- Return interpolated addressing
                virtual const labelListList& addressing() const
                {
                    return addressing_;
                }

                //- Return interpolation weights
                virtual const scalarListList& weights() const
                {
                    return weights_;
                }
        };


    // Constructors

        //- Construct from fvMesh and no agglomeration
        singleCellFvMesh(const IOobject& io, const fvMesh&);

        //- Construct from fvMesh and agglomeration of boundary faces.
        //  agglomeration is per patch, per patch face index the agglomeration
        //  the face goes into.
        singleCellFvMesh
        (
            const IOobject& io,
            const fvMesh&,
            const labelListList& patchFaceAgglomeration
        );

        //- Read from IOobject
        singleCellFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        singleCellFvMesh(const singleCellFvMesh&) = delete;


    // Member Functions

        bool agglomerate() const
        {
            return patchFaceAgglomeration_.size() > 0;
        }

        //- From patchFace on this back to original mesh or agglomeration
        const labelListList& patchFaceMap() const
        {
            return patchFaceMap_;
        }

        //- From point on this back to original mesh
        const labelList& pointMap() const
        {
            return pointMap_;
        }

        //- From face on original mesh to face on this
        const labelList& reverseFaceMap() const
        {
            return reverseFaceMap_;
        }

        //- From point on original mesh to point on this (or -1 for removed
        //  points)
        const labelList& reversePointMap() const
        {
            return reversePointMap_;
        }

        //- Map volField. Internal field set to average, patch fields straight
        //  copies.
        template<class Type>
        tmp<VolField<Type>>
        interpolate
        (
            const VolField<Type>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singleCellFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleCellFvMeshInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
