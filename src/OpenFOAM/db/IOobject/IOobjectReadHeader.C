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

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::IOobject::headerOk()
{
    return typeHeaderOk<IOobject>(false);
}


bool Foam::IOobject::readHeader(Istream& is)
{
    if (IOobject::debug)
    {
        InfoInFunction << "Reading header for file " << is.name() << endl;
    }

    // Check Istream not already bad
    if (!is.good())
    {
        if (rOpt_ == MUST_READ || rOpt_ == MUST_READ_IF_MODIFIED)
        {
            FatalIOErrorInFunction(is)
                << " stream not open for reading essential object from file "
                << is.name()
                << exit(FatalIOError);
        }

        if (IOobject::debug)
        {
            SeriousIOErrorInFunction(is)
                << " stream not open for reading from file "
                << is.name() << endl;
        }

        return false;
    }

    token firstToken(is);

    if
    (
        is.good()
     && firstToken.isWord()
     && firstToken.wordToken() == foamFile
    )
    {
        dictionary headerDict(is);

        IOstream::versionNumber ver = IOstream::currentVersion;
        headerDict.readIfPresent("version", ver);
        is.version(ver);

        is.format(headerDict.lookup("format"));

        headerClassName_ = word(headerDict.lookup("class"));

        const word headerObject(headerDict.lookup("object"));
        if (IOobject::debug && headerObject != name())
        {
            IOWarningInFunction(is)
                << " object renamed from "
                << name() << " to " << headerObject
                << " for file " << is.name() << endl;
        }

        // The note entry is optional
        headerDict.readIfPresent("note", note_);
    }
    else
    {
        if (IOobject::debug)
        {
            IOWarningInFunction(is)
                << "First token could not be read "
                   "or is not the keyword 'FoamFile'"
                << nl << nl << "Check header is of the form:" << nl << endl;

            writeHeader(Info);
        }

        return false;
    }

    // Check stream is still OK
    if (is.good())
    {
        objState_ = GOOD;
    }
    else
    {
        if (rOpt_ == MUST_READ || rOpt_ == MUST_READ_IF_MODIFIED)
        {
            FatalIOErrorInFunction(is)
                << " stream failure while reading header"
                << " on line " << is.lineNumber()
                << " of file " << is.name()
                << " for essential object" << name()
                << exit(FatalIOError);
        }

        if (IOobject::debug)
        {
            InfoInFunction
                << "Stream failure while reading header"
                << " on line " << is.lineNumber()
                << " of file " << is.name() << endl;
        }

        objState_ = BAD;

        return false;
    }

    if (IOobject::debug)
    {
        Info<< " .... read" << endl;
    }

    return true;
}


// ************************************************************************* //
