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
    Foam::IStringStream

Description
    Input from memory buffer stream.

SourceFiles
    StringStreamsPrint.C

\*---------------------------------------------------------------------------*/

#ifndef IStringStream_H
#define IStringStream_H

#include "ISstream.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IStringStream Declaration
\*---------------------------------------------------------------------------*/

class IStringStream
:
    public ISstream
{

public:

    // Constructors

        //- Construct from string
        IStringStream
        (
            const string& buffer,
            const streamFormat format = ASCII,
            const versionNumber version = currentVersion
        )
        :
            ISstream
            (
                *(new std::istringstream(buffer)),
                "IStringStream.sourceFile",
                format,
                version
            )
        {}


        //- Construct from string and name
        IStringStream
        (
            const string& name,
            const string& buffer,
            const streamFormat format = ASCII,
            const versionNumber version = currentVersion
        )
        :
            ISstream
            (
                *(new std::istringstream(buffer)),
                name,
                format,
                version
            )
        {}


        //- Construct from char*
        IStringStream
        (
            const char* buffer,
            const streamFormat format = ASCII,
            const versionNumber version = currentVersion
        )
        :
            ISstream
            (
                *(new std::istringstream(buffer)),
                "IStringStream.sourceFile",
                format,
                version
            )
        {}


    //- Destructor
    ~IStringStream()
    {
        delete &dynamic_cast<std::istringstream&>(stdStream());
    }


    // Member Functions

        // Access

            //- Return the string
            string str() const
            {
                return dynamic_cast<const std::istringstream&>
                (
                    stdStream()
                ).str();
            }


        // Print

            //- Print description to Ostream
            void print(Ostream&) const;


    // Member Operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary:
        //  e.g. thing thisThing(IFstream("thingFileName")());
        Istream& operator()() const
        {
            return const_cast<IStringStream&>(*this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
