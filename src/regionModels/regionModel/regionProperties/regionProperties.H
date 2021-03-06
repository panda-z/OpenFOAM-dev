/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::regionProperties

Description
    Simple class to hold region information for coupled region simulations.

    Gives per physics ('fluid', 'solid') the names of the regions. There
    is no assumption on this level that one region should only have one
    set of physics.

SourceFiles
    regionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef regionProperties_H
#define regionProperties_H

#include "Time.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regionProperties Declaration
\*---------------------------------------------------------------------------*/

class regionProperties
:
    public HashTable<wordList>
{
public:

    // Constructors

        //- Construct from components
        regionProperties(const Time& runTime);

        //- Disallow default bitwise copy construction
        regionProperties(const regionProperties&) = delete;


    //- Destructor
    ~regionProperties();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const regionProperties&) = delete;
};


const word& regionDir(const word& regionName);

wordList selectRegionNames(const argList& args, const Time& runTime);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
