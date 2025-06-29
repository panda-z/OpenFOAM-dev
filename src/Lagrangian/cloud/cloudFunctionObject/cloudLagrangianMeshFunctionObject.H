/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionObjects::cloudLagrangianMeshFunctionObject

Description
    Base class for function objects that refer to a LagrangianMesh and a cloud.
    Used, for example, to construct, evolve and write derived Lagrangian
    properties.

SourceFiles
    cloudLagrangianMeshFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef cloudLagrangianMeshFunctionObject_functionObject_H
#define cloudLagrangianMeshFunctionObject_functionObject_H

#include "LagrangianMeshFunctionObject.H"
#include "cloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
              Class cloudLagrangianMeshFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class cloudLagrangianMeshFunctionObject
:
    public LagrangianMeshFunctionObject,
    public cloudFunctionObject
{
public:

    //- Runtime type information
    TypeName("cloudLagrangianMeshFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        cloudLagrangianMeshFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudLagrangianMeshFunctionObject
        (
            const cloudLagrangianMeshFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~cloudLagrangianMeshFunctionObject();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudLagrangianMeshFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
