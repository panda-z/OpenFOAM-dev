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

\*---------------------------------------------------------------------------*/

#include "fvFieldDecomposer.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::fvFieldDecomposer::completePatchID
(
    const label proci,
    const label procPatchi
) const
{
    const fvPatch& procPatch = procMeshes_[proci].boundary()[procPatchi];

    if (procPatchi < completeMesh_.boundary().size())
    {
        return procPatchi;
    }
    else if (isA<processorCyclicFvPatch>(procPatch))
    {
        return
            refCast<const processorCyclicFvPatch>(procPatch)
           .referPatchIndex();
    }
    else
    {
        return -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const labelUList& addressing
)
:
    labelList(mag(addressing) - 1),
    forwardFieldMapper(static_cast<const labelList&>(*this))
{}


Foam::fvFieldDecomposer::fvFieldDecomposer
(
    const fvMesh& completeMesh,
    const PtrList<fvMesh>& procMeshes,
    const labelListList& faceProcAddressing,
    const labelListList& cellProcAddressing,
    const PtrList<surfaceLabelField::Boundary>& faceProcAddressingBf
)
:
    completeMesh_(completeMesh),
    procMeshes_(procMeshes),
    faceProcAddressing_(faceProcAddressing),
    cellProcAddressing_(cellProcAddressing),
    faceProcAddressingBf_(faceProcAddressingBf),
    patchFieldDecomposers_(procMeshes_.size())
{
    forAll(procMeshes_, proci)
    {
        patchFieldDecomposers_.set
        (
            proci,
            new PtrList<patchFieldDecomposer>
            (
                procMeshes_[proci].boundary().size()
            )
        );

        forAll(procMeshes_[proci].boundary(), procPatchi)
        {
            const label completePatchi = completePatchID(proci, procPatchi);

            if (completePatchi >= 0)
            {
                patchFieldDecomposers_[proci].set
                (
                    procPatchi,
                    new patchFieldDecomposer
                    (
                        faceProcAddressingBf[proci][completePatchi]
                    )
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvFieldDecomposer::~fvFieldDecomposer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvFieldDecomposer::decomposes(const IOobjectList& objects)
{
    bool result = false;

    #define DO_FV_FIELDS_TYPE(Type, nullArg)                                   \
        result = result                                                        \
         || !objects.lookupClass(VolField<Type>::Internal::typeName).empty()   \
         || !objects.lookupClass(VolField<Type>::typeName).empty()             \
         || !objects.lookupClass(SurfaceField<Type>::typeName).empty();
    FOR_ALL_FIELD_TYPES(DO_FV_FIELDS_TYPE)
    #undef DO_FV_FIELDS_TYPE

    return result;
}


// ************************************************************************* //
