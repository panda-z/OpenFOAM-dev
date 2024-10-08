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

\*---------------------------------------------------------------------------*/

#include "singleCellFvMesh.H"
#include "calculatedFvPatchFields.H"
#include "identityFieldMapper.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>> singleCellFvMesh::interpolate
(
    const VolField<Type>& vf
) const
{
    // 1. Create the complete field with dummy patch fields
    PtrList<fvPatchField<Type>> patchFields(vf.boundaryField().size());

    forAll(patchFields, patchi)
    {
        patchFields.set
        (
            patchi,
            fvPatchField<Type>::New
            (
                calculatedFvPatchField<Type>::typeName,
                boundary()[patchi],
                DimensionedField<Type, volMesh>::null()
            )
        );
    }

    // Create the complete field from the pieces
    tmp<VolField<Type>> tresF
    (
        new VolField<Type>
        (
            IOobject
            (
                vf.name(),
                time().name(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this,
            vf.dimensions(),
            Field<Type>(1, gAverage(vf)),
            patchFields
        )
    );
    VolField<Type>& resF = tresF.ref();


    // 2. Change the fvPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    typename VolField<Type>::
        Boundary& bf = resF.boundaryFieldRef();

    if (agglomerate())
    {
        forAll(vf.boundaryField(), patchi)
        {
            const labelList& agglom = patchFaceAgglomeration_[patchi];
            label nAgglom = max(agglom)+1;

            // Use inverse of agglomeration. This is from agglomeration to
            // original (fine) mesh patch face.
            labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
            inplaceReorder(patchFaceMap_[patchi], coarseToFine);
            scalarListList coarseWeights(nAgglom);
            forAll(coarseToFine, coarseI)
            {
                const labelList& fineFaces = coarseToFine[coarseI];
                coarseWeights[coarseI] = scalarList
                (
                    fineFaces.size(),
                    1.0/fineFaces.size()
                );
            }

            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[patchi],
                    boundary()[patchi],
                    resF(),
                    agglomPatchFieldMapper(coarseToFine, coarseWeights)
                )
            );
        }
    }
    else
    {
        forAll(vf.boundaryField(), patchi)
        {
            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[patchi],
                    boundary()[patchi],
                    resF(),
                    identityFieldMapper()
                )
            );
        }
    }

    return tresF;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
