/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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
    Foam::Function1s::quarterCosineRamp

Description
    Quarter-cosine ramp function starting from 0 and increasing to 1 from \c
    start over the \c duration and remaining at 1 thereafter.

See also
    Foam::Function1s::Ramp

SourceFiles
    quarterCosineRamp.C

\*---------------------------------------------------------------------------*/

#ifndef quarterCosineRamp_H
#define quarterCosineRamp_H

#include "Ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                      Class quarterCosineRamp Declaration
\*---------------------------------------------------------------------------*/

class quarterCosineRamp
:
    public Ramp<quarterCosineRamp>
{

public:

    // Runtime type information
    TypeName("quarterCosineRamp");


    // Constructors

        //- Construct from name and dictionary
        quarterCosineRamp
        (
            const word& name,
            const unitConversions& units,
            const dictionary& dict
        );


    //- Destructor
    virtual ~quarterCosineRamp();


    // Member Functions

        //- Return value for time t
        virtual inline scalar value(const scalar t) const;

        //- Return the integral between times t1 and t2
        virtual inline scalar integral
        (
            const scalar t1,
            const scalar t2
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const quarterCosineRamp&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "quarterCosineRampI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
