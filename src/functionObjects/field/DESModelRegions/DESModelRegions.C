/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "DESModelRegions.H"
#include "Time.H"
#include "fvMesh.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "DESModelBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(DESModelRegions, 0);
    addToRunTimeSelectionTable(functionObject, DESModelRegions, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void Foam::functionObjects::DESModelRegions::writeFileHeader(const label i)
{
    writeHeader(file(), "DES model region coverage (% volume)");

    writeCommented(file(), "Time");
    writeTabbed(file(), "LES");
    writeTabbed(file(), "RAS");

    file().endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::DESModelRegions::DESModelRegions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log)
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::DESModelRegions::~DESModelRegions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::DESModelRegions::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    return true;
}


bool Foam::functionObjects::DESModelRegions::execute()
{
    if (mesh_.foundObject<DESModelBase>(turbulenceModel::propertiesName))
    {
        const DESModelBase& model =
            mesh_.lookupObject<DESModelBase>
            (
                turbulenceModel::propertiesName
            );

        tmp<volScalarField> tLESRegion = model.LESRegion();

        word name(type());

        return store(name, tLESRegion);
    }
    else
    {
        Log << "    No DES turbulence model found in database" << nl
            << endl;

        return false;
    }
}


bool Foam::functionObjects::DESModelRegions::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& LESRegion = obr_.lookupObject<volScalarField>(type());

    scalar prc =
        gSum(LESRegion.primitiveField()*mesh_.V())
        /gSum(mesh_.V())*100.0;

    if (Pstream::master())
    {
        file() << obr_.time().value()
            << tab << prc
            << tab << 100.0 - prc
            << endl;
    }

    Log << "    LES = " << prc << " % (volume)" << nl
        << "    RAS = " << 100.0 - prc << " % (volume)" << nl
        << endl;

    return true;
}


// ************************************************************************* //
