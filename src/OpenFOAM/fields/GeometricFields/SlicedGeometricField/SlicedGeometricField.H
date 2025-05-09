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
    Foam::SlicedGeometricField

Description
    Specialisation of GeometricField which holds slices of given complete
    fields in a form that they act as a GeometricField.

    The destructor is wrapped to avoid deallocation of the storage of the
    complete fields when this is destroyed.

    SlicedGeometricField can only be instantiated with a valid form of
    SlicedPatchField to handle the slicing and storage deallocation of the
    boundary field.

SourceFiles
    SlicedGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef SlicedGeometricField_H
#define SlicedGeometricField_H

#include "GeometricField.H"
#include "SlicedDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class SlicedGeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class SlicedGeometricField
:
    public GeometricField<Type, GeoMesh>
{
public:

    // Public Typedefs

        //- Type of mesh on which this SlicedGeometricField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of the internal field from which this SlicedGeometricField is
        //  derived
        typedef SlicedDimensionedField<Type, GeoMesh> Internal;

        //- Type of the patch field of which the Boundary is composed
        typedef typename GeoMesh::template PatchField<Type> Patch;

        //- Type of the sliced patch field of which the Boundary is composed
        typedef typename Patch::Sliced SlicedPatch;


private:

    // Private Member Functions

        //- Slice the given field and a create a PtrList of sliced patch fields
        //  from which the boundary field is built
        tmp<FieldField<GeoMesh::template PatchField, Type>> slicedBoundaryField
        (
            const Mesh& mesh,
            const Field<Type>& completeField,
            const bool preserveCouples,
            const bool preserveProcessorOnly = false
        );

        //- Slice the given field and a create a PtrList of sliced patch fields
        //  from which the boundary field is built
        tmp<FieldField<GeoMesh::template PatchField, Type>> slicedBoundaryField
        (
            const Mesh& mesh,
            const FieldField<GeoMesh::template PatchField, Type>& bField,
            const bool preserveCouples
        );


public:

    // Constructors

        //- Construct from components and field to slice
        SlicedGeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const Field<Type>& completeField,
            const bool preserveCouples=true
        );

        //- Construct from components and separate fields to slice for the
        //  internal field and boundary field
        SlicedGeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const Field<Type>& completeIField,
            const Field<Type>& completeBField,
            const bool preserveCouples=true,
            const bool preserveProcessorOnly = false
        );

        //- Construct from GeometricField. Reuses full internal and
        //  patch fields except on couples (preserveCouples=true).
        SlicedGeometricField
        (
            const IOobject&,
            const GeometricField<Type, GeoMesh>&,
            const bool preserveCouples=true
        );

        //- Copy constructor
        SlicedGeometricField(const SlicedGeometricField<Type, GeoMesh>&);

        //- Clone
        tmp<SlicedGeometricField<Type, GeoMesh>> clone() const;


    //- Destructor
    ~SlicedGeometricField();


    // Member Functions

        //- Splice the sliced field and return the complete field
        tmp<Field<Type>> splice() const;

        //- Correct boundary field
        void correctBoundaryConditions();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SlicedGeometricField&) = delete;

        //- Disallow standard assignment to GeometricField,
        //  == assignment is allowed.
        void operator=(const GeometricField<Type, GeoMesh>&) = delete;

        //- Disallow standard assignment to tmp<GeometricField>,
        //  == assignment is allowed.
        void operator=(const tmp<GeometricField<Type, GeoMesh>>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SlicedGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
