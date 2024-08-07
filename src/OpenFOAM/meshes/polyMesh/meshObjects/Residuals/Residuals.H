/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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
    Foam::Residuals

Description
    DemandDrivenMeshObject to store the solver performance residuals of all
    the fields of the type it is instantiated on.

SourceFiles
    Residuals.C

\*---------------------------------------------------------------------------*/

#ifndef Residuals_H
#define Residuals_H

#include "DemandDrivenMeshObject.H"
#include "polyMesh.H"
#include "solverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Residuals Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Residuals
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        TopoChangeableMeshObject,
        Residuals<Type>
    >,
    public HashTable<DynamicList<SolverPerformance<Type>>>
{
    // Private Data

        //- Previously used time-index, used for reset between iterations
        mutable label prevTimeIndex_;


protected:

    friend class DemandDrivenMeshObject
    <
        polyMesh,
        TopoChangeableMeshObject,
        Residuals<Type>
    >;

    // Protected Constructors

        //- Construct for given mesh
        explicit Residuals(const polyMesh& mesh);


public:

    //- Runtime type information
    TypeName("residuals");


    // Constructors

        //- Disallow default bitwise copy construction
        Residuals(const Residuals<Type>&) = delete;


    // Member Functions

        //- Return the list of field names of the particular type
        //  for which residuals are stored
        static List<word> fieldNames(const polyMesh& mesh);

        //- Return true if residuals for the given field are stored
        static bool found(const polyMesh& mesh, const word& fieldName);

        //- Return the list of solver performance residuals for the given field
        static const DynamicList<SolverPerformance<Type>>& field
        (
            const polyMesh& mesh,
            const word& fieldName
        );

        //- Append the given solver performance residuals
        //  in the corresponding list
        static void append
        (
            const polyMesh& mesh,
            const SolverPerformance<Type>&
        );

        //- Update for mesh motion
        virtual bool movePoints();

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap& map);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap& map);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Residuals<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Residuals.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
