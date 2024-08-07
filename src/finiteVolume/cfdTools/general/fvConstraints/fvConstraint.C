/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "fvConstraint.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvConstraint, 0);
    defineRunTimeSelectionTable(fvConstraint, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::fvConstraint::constrainType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    return false;
}


template<class Type>
bool Foam::fvConstraint::constrainType(VolField<Type>& field) const
{
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvConstraint::fvConstraint
(
    const word& name,
    const word& constraintType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    constraintType_(constraintType),
    mesh_(mesh)
{
    Info<< incrIndent << indent << "Name: " << name_
        << endl << decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvConstraint> Foam::fvConstraint::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word constraintType(dict.lookup("type"));

    Info<< indent
        << "Selecting finite volume constraint type " << constraintType << endl;

    if
    (
        !dictionaryConstructorTablePtr_
     || dictionaryConstructorTablePtr_->find(constraintType)
        == dictionaryConstructorTablePtr_->end()
    )
    {
        if
        (
           !libs.open
            (
                dict,
                "libs",
                dictionaryConstructorTablePtr_
            )
        )
        {
            libs.open("lib" + constraintType.remove(':') + ".so", false);
        }

        if (!dictionaryConstructorTablePtr_)
        {
            FatalErrorInFunction
                << "Unknown constraint type "
                << constraintType << nl << nl
                << "Table of fvConstraints is empty"
                << exit(FatalError);
        }
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(constraintType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown fvConstraint " << constraintType << nl << nl
            << "Valid fvConstraints are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<fvConstraint>
    (
        cstrIter()(name, constraintType, mesh, dict)
    );
}


Foam::fvConstraint::~fvConstraint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fvConstraint::constrainedFields() const
{
    return wordList::null();
}


bool Foam::fvConstraint::constrainsField(const word& fieldName) const
{
    return findIndex(constrainedFields(), fieldName) != -1;
}


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_CONSTRAINT_CONSTRAIN, fvConstraint);


FOR_ALL_FIELD_TYPES(IMPLEMENT_FV_CONSTRAINT_CONSTRAIN_FIELD, fvConstraint);


bool Foam::fvConstraint::read(const dictionary& dict)
{
    return true;
}


// ************************************************************************* //
