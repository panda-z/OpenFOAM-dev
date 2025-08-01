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

#include "histogram.H"
#include "volFields.H"
#include "setWriter.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(histogram, 0);
    addToRunTimeSelectionTable(functionObject, histogram, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::histogram::histogram
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    zone_(mesh(), dict),
    file_(obr_, name)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::histogram::~histogram()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::histogram::read(const dictionary& dict)
{
    dict.lookup("field") >> fieldName_;
    dict.lookup("max") >> max_;
    min_ = dict.lookupOrDefault<scalar>("min", 0);
    dict.lookup("nBins") >> nBins_;

    formatterPtr_ = setWriter::New(dict.lookup("setFormat"), dict);

    return true;
}


Foam::wordList Foam::functionObjects::histogram::fields() const
{
    return wordList(fieldName_);
}


bool Foam::functionObjects::histogram::execute()
{
    return true;
}


bool Foam::functionObjects::histogram::write()
{
    Log << type() << " " << name() << " write:" << nl;

    const volScalarField& field
    (
        obr_.lookupObject<volScalarField>(fieldName_)
    );

    // Calculate the mid-points of bins for the graph axis
    scalarField xBin(nBins_);
    const scalar delta = (max_- min_)/nBins_;
    scalar x = min_ + 0.5*delta;
    forAll(xBin, i)
    {
        xBin[i] = x;
        x += delta;
    }

    scalarField volFrac(nBins_, 0);
    const scalarField& V = mesh_.V();

    if (zone_.all())
    {
        forAll(field, celli)
        {
            const label bini = (field[celli] - min_)/delta;
            if (bini >= 0 && bini < nBins_)
            {
                volFrac[bini] += V[celli];
            }
        }
    }
    else
    {
        const labelList& zoneCells = zone_.zone();
        forAll(zoneCells, i)
        {
            const label celli = zoneCells[i];
            const label bini = (field[celli] - min_)/delta;
            if (bini >= 0 && bini < nBins_)
            {
                volFrac[bini] += V[celli];
            }
        }
    }

    Pstream::listCombineGather(volFrac, plusEqOp<scalar>());

    if (Pstream::master())
    {
        const scalar sumVol = sum(volFrac);

        if (sumVol > small)
        {
            volFrac /= sumVol;

            formatterPtr_().write
            (
                file_.baseTimeDir(),
                typeName,
                coordSet(true, fieldName_, xBin),
                "v/vTotal",
                volFrac
            );
        }
    }

    return true;
}


void Foam::functionObjects::histogram::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &this->mesh())
    {
        zone_.movePoints();
    }
}


void Foam::functionObjects::histogram::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.topoChange(map);
    }
}


void Foam::functionObjects::histogram::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.mapMesh(map);
    }
}


void Foam::functionObjects::histogram::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh())
    {
        zone_.distribute(map);
    }
}


// ************************************************************************* //
