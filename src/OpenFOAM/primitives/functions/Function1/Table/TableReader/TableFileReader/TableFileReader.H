/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2024 OpenFOAM Foundation
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
    Foam::TableFileReader

Description
    Base class to read table data for tables

SourceFiles
    TableFileReader.C

\*---------------------------------------------------------------------------*/

#ifndef TableFileReader_H
#define TableFileReader_H

#include "TableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class TableFileReader Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TableFileReader
:
    public TableReader<Type>
{
    // Private Data

        //- File name for table
        const fileName fName_;

        //- Unit conversions
        const autoPtr<Function1s::unitConversions> unitsPtr_;


    // Private Member Functions

        //- Read the units, if any
        autoPtr<Function1s::unitConversions> readUnits
        (
            const Function1s::unitConversions& defaultUnits,
            const dictionary& dict
        ) const;

        //- Read a 1D table
        virtual void read(ISstream&, List<Tuple2<scalar, Type>>&) const = 0;


protected:

    // Protected Member Functions

        //- Read a 1D table
        void read
        (
            const Function1s::unitConversions& defaultUnits,
            const dictionary& dict,
            List<Tuple2<scalar, Type>>&
        ) const;


public:

    // Constructors

        //- Construct from dictionary
        TableFileReader
        (
            const Function1s::unitConversions& defaultUnits,
            const dictionary& dict
        );

        //- Copy construct
        TableFileReader(const TableFileReader<Type>&);


    //- Destructor
    virtual ~TableFileReader();


    // Member Functions

        //- Read values
        virtual List<Tuple2<scalar, Type>> read
        (
            const Function1s::unitConversions& units,
            const dictionary& dict,
            const word& valuesKeyword=word::null
        ) const;

        //- Write additional information
        virtual void write
        (
            Ostream& os,
            const Function1s::unitConversions& units,
            const List<Tuple2<scalar, Type>>& table,
            const word& valuesKeyword=word::null
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TableFileReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
