/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
    Foam::interfaceCompositionModels::saturated

Description
    Model which uses a saturation pressure model for a single species to
    calculate the interface composition.

SourceFiles
    saturated.C

\*---------------------------------------------------------------------------*/

#ifndef saturated_H
#define saturated_H

#include "interfaceCompositionModel.H"
#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceCompositionModels
{

/*---------------------------------------------------------------------------*\
                          Class saturated Declaration
\*---------------------------------------------------------------------------*/

class saturated
:
    public interfaceCompositionModel
{
protected:

    // Private Data

        //- saturated species name
        word saturatedName_;

        //- saturated species index
        label saturatedIndex_;

        //- Saturation pressure model
        autoPtr<saturationPressureModel> saturationModel_;


    // Private Member Functions

        //- Constant of proportionality between partial pressure and mass
        //  fraction
        tmp<volScalarField> wRatioByP() const;


public:

    //- Runtime type information
    TypeName("saturated");


    // Constructors

        //- Construct from a dictionary and an interface
        saturated
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~saturated();


    // Member Functions

        //- Update the composition
        virtual void update(const volScalarField& Tf);

        //- The interface species fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- The interface species fraction derivative w.r.t. temperature
        virtual tmp<volScalarField> YfPrime
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceCompositionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
