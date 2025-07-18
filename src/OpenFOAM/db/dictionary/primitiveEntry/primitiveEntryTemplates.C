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

#include "primitiveEntry.H"
#include "dictionary.H"
#include "IStringStream.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const T& t
)
:
    entry(key),
    ITstream(key, tokenList(10))
{
    OStringStream os;
    os  << t << token::END_STATEMENT;
    readEntry(dictionary::null, IStringStream(os.str())());
}


template<class T>
Foam::primitiveEntry::primitiveEntry
(
    const keyType& key,
    const T& t,
    const label startLineNumber,
    const label endLineNumber
)
:
    entry(key, startLineNumber),
    ITstream(key, tokenList(10))
{
    OStringStream os;
    os  << t << token::END_STATEMENT;
    IStringStream iss(os.str());

    if (endLineNumber != -1)
    {
        iss.lineNumber() = endLineNumber;
    }
    else
    {
        iss.lineNumber() = startLineNumber;
    }

    readEntry(dictionary::null, iss());
}


// ************************************************************************* //
