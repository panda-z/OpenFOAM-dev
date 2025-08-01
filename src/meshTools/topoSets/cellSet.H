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
    Foam::cellSet

Description
    A collection of cell labels.

SourceFiles
    cellSet.C

\*---------------------------------------------------------------------------*/

#ifndef cellSet_H
#define cellSet_H

#include "topoSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellSet Declaration
\*---------------------------------------------------------------------------*/

class cellSet
:
    public topoSet
{
public:

    //- Runtime type information
    TypeName("cellSet");


    // Constructors

        //- Construct from IOobject. No checking.
        cellSet(const IOobject& obj);

        //- Construct from polyMesh and name. Checks for valid cell ids.
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            readOption r=MUST_READ,
            writeOption w=NO_WRITE
        );

        //- Construct empty from size of labelHashSet
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const label sizes,
            writeOption w=NO_WRITE
        );

        //- Construct from existing set
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet&,
            writeOption w=NO_WRITE
        );

        //- Construct from labelHashSet
        cellSet
        (
            const polyMesh& mesh,
            const word& name,
            const labelHashSet&,
            writeOption w=NO_WRITE
        );

        //- Disallow default bitwise copy construction
        cellSet(const cellSet&) = delete;


    //- Destructor
    virtual ~cellSet();


    // Member Functions

        //- Sync cellSet across coupled patches.
        virtual void sync(const polyMesh& mesh)
        {}

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
