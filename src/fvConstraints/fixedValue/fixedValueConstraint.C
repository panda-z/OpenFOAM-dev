/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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

#include "fixedValueConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvcSurfaceIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fixedValueConstraint, 0);

    addToRunTimeSelectionTable
    (
        fvConstraint,
        fixedValueConstraint,
        dictionary
    );
    addBackwardCompatibleToRunTimeSelectionTable
    (
        fvConstraint,
        fixedValueConstraint,
        dictionary,
        fixedValueConstraint,
        "fixedValueConstraint"
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::fixedValueConstraint::readCoeffs(const dictionary& dict)
{
    fieldValues_.clear();
    forAllConstIter(dictionary, dict.subDict("fieldValues"), iter)
    {
        fieldValues_.set
        (
            iter().keyword(),
            new unknownTypeFunction1
            (
                iter().keyword(),
                mesh().time().userUnits(),
                dict.subDict("fieldValues")
            )
        );
    }

    fraction_ =
        dict.found("fraction")
      ? Function1<scalar>::New
        (
            "fraction",
            mesh().time().userUnits(),
            unitFraction,
            dict
        )
      : autoPtr<Function1<scalar>>();
}


template<class Type>
bool Foam::fv::fixedValueConstraint::constrainType
(
    fvMatrix<Type>& eqn,
    const word& fieldName
) const
{
    // Set the value units for the function
    fieldValues_[fieldName]->template setValueUnits<Type>
    (
        eqn.psi().dimensions()
    );

    const scalar t = mesh().time().value();

    const List<Type> values
    (
        zone_.nCells(),
        fieldValues_[fieldName]->value<Type>(t)
    );

    if (fraction_.valid())
    {
        eqn.setValues
        (
            zone_.zone(),
            values,
            scalarList(zone_.nCells(), fraction_->value(t))
        );
    }
    else
    {
        eqn.setValues(zone_.zone(), values);
    }

    return zone_.nCells();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedValueConstraint::fixedValueConstraint
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    zone_(mesh, coeffs(dict))
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::fixedValueConstraint::constrainedFields() const
{
    return fieldValues_.toc();
}


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FV_CONSTRAINT_CONSTRAIN,
    fv::fixedValueConstraint
);


bool Foam::fv::fixedValueConstraint::movePoints()
{
    zone_.movePoints();
    return true;
}


void Foam::fv::fixedValueConstraint::topoChange(const polyTopoChangeMap& map)
{
    zone_.topoChange(map);
}


void Foam::fv::fixedValueConstraint::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::fixedValueConstraint::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
}


bool Foam::fv::fixedValueConstraint::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        zone_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
