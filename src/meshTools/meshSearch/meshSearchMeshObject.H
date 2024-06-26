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
    Foam::meshSearchMeshObject

Description
    DemandDrivenMeshObject wrapper around meshSearch(mesh).

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef meshSearchMeshObject_H
#define meshSearchMeshObject_H

#include "DemandDrivenMeshObject.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshSearchMeshObject Declaration
\*---------------------------------------------------------------------------*/

class meshSearchMeshObject
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        DeletableMeshObject,
        meshSearchMeshObject
    >,
    public meshSearch
{

protected:

    friend class DemandDrivenMeshObject
    <
        polyMesh,
        DeletableMeshObject,
        meshSearchMeshObject
    >;

    // Protected Constructors

        //- Constructor from polyMesh
        explicit meshSearchMeshObject(const polyMesh& mesh);


public:

    // Declare name of the class and its debug switch
    TypeName("meshSearchMeshObject");

    //- Destructor
    virtual ~meshSearchMeshObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
