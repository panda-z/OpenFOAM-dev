/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

Description
    Perform mapping of finite volume fields required by stitching.

\*---------------------------------------------------------------------------*/

#include "fvMeshStitcher.H"
#include "conformedFvPatchField.H"
#include "conformedFvsPatchField.H"
#include "nonConformalErrorFvPatch.H"
#include "setSizeFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class GeoField>
void Foam::fvMeshStitcher::resizePatchFields()
{
    UPtrList<GeoField<Type>> fields(mesh_.fields<GeoField<Type>>());

    forAll(fields, i)
    {
        forAll(mesh_.boundary(), patchi)
        {
            typename GeoField<Type>::Patch& pf =
                fields[i].boundaryFieldRefNoStoreOldTimes()[patchi];

            if (isA<nonConformalFvPatch>(pf.patch()))
            {
                pf.map(pf, setSizeFieldMapper(pf.patch().size()));
            }
        }
    }
}


template<template<class> class GeoField>
void Foam::fvMeshStitcher::resizePatchFields()
{
    #define ResizePatchFields(Type, nullArg) \
        resizePatchFields<Type, GeoField>();
    FOR_ALL_FIELD_TYPES(ResizePatchFields);
    #undef ResizePatchFields
}


template<class Type>
void Foam::fvMeshStitcher::preConformSurfaceFields()
{
    UPtrList<SurfaceField<Type>> fields(mesh_.fields<SurfaceField<Type>>());

    forAll(fields, i)
    {
        conformedFvsPatchField<Type>::conform
        (
            fields[i].boundaryFieldRefNoStoreOldTimes()
        );
    }
}


template<class Type>
void Foam::fvMeshStitcher::preConformVolFields()
{
    UPtrList<VolField<Type>> fields(mesh_.fields<VolField<Type>>());

    forAll(fields, i)
    {
        conformedFvPatchField<Type>::conform
        (
            fields[i].boundaryFieldRefNoStoreOldTimes()
        );
    }
}


template<class Type>
void Foam::fvMeshStitcher::postUnconformSurfaceFields()
{
    if (mesh_.topoChanged())
    {
        UPtrList<SurfaceField<Type>> curFields
        (
            mesh_.curFields<SurfaceField<Type>>()
        );

        forAll(curFields, i)
        {
            curFields[i].clearOldTimes();
        }
    }

    UPtrList<SurfaceField<Type>> fields(mesh_.fields<SurfaceField<Type>>());

    forAll(fields, i)
    {
        conformedFvsPatchField<Type>::unconform
        (
            fields[i].boundaryFieldRefNoStoreOldTimes()
        );
    }
}


template<class Type>
void Foam::fvMeshStitcher::postUnconformVolFields()
{
    UPtrList<VolField<Type>> fields(mesh_.fields<VolField<Type>>());

    forAll(fields, i)
    {
        conformedFvPatchField<Type>::unconform
        (
            fields[i].boundaryFieldRefNoStoreOldTimes()
        );
    }
}


template<class Type>
void Foam::fvMeshStitcher::postUnconformEvaluateVolFields()
{
    auto evaluate = [](const typename VolField<Type>::Patch& pf)
    {
        return
            (
                isA<nonConformalFvPatch>(pf.patch())
             && pf.type() == pf.patch().patch().type()
             && polyPatch::constraintType(pf.patch().patch().type())
            )
         || isA<nonConformalErrorFvPatch>(pf.patch());
    };

    UPtrList<VolField<Type>> fields(mesh_.fields<VolField<Type>>());

    forAll(fields, i)
    {
        const label nReq = Pstream::nRequests();

        forAll(mesh_.boundary(), patchi)
        {
            typename VolField<Type>::Patch& pf =
                fields[i].boundaryFieldRefNoStoreOldTimes()[patchi];

            if (evaluate(pf))
            {
                pf.initEvaluate(Pstream::defaultCommsType);
            }
        }

        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(mesh_.boundary(), patchi)
        {
            typename VolField<Type>::Patch& pf =
                fields[i].boundaryFieldRefNoStoreOldTimes()[patchi];

            if (evaluate(pf))
            {
                pf.evaluate(Pstream::defaultCommsType);
            }
        }
    }
}


// ************************************************************************* //
