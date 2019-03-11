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

Class
    Foam::functionObjects::DESModelRegions

Description
    This function object writes out an indicator field for DES turbulence
    calculations, that is:
    - 0 for RAS regions
    - 1 for LES regions

    Example of function object specification:
    \verbatim
    DESModelRegions1
    {
        type           DESModelRegions;
        libs           ("libfieldFunctionObject.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: DESModelRegions | yes      |
    \endtable

See also
    Foam::functionObject

SourceFiles
    DESModelRegions.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_DESModelRegions_H
#define functionObjects_DESModelRegions_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class DESModelRegions Declaration
\*---------------------------------------------------------------------------*/

class DESModelRegions
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private Member Functions


        //- Disallow default bitwise copy construct
        DESModelRegions(const DESModelRegions&);

        //- Disallow default bitwise assignment
        void operator=(const DESModelRegions&);


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("DESModelRegions");


    // Constructors

        //- Construct from Time and dictionary
        DESModelRegions
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~DESModelRegions();


    // Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary&);

        //- Calculate the DESModelRegions indicator
        virtual bool execute();

        //- Write the DESModelRegions field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
