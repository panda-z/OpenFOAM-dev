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
    Foam::saturationTemperatureModel

Description
    Model to describe the dependence of saturation temperature on pressure

SourceFiles
    saturationTemperatureModel.C
    saturationTemperatureModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef saturationTemperatureModel_H
#define saturationTemperatureModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DEFINE_TSAT(FieldType, Modifier)                                       \
                                                                               \
    virtual tmp<FieldType> Tsat(const FieldType& p) const Modifier;            \
                                                                               \
    virtual tmp<FieldType> TsatPrime(const FieldType& p) const Modifier;       \


#define IMPLEMENT_TSAT(ModelType, FieldType)                                   \
                                                                               \
    Foam::tmp<Foam::FieldType>                                                 \
    Foam::ModelType::Tsat(const FieldType& p) const                            \
    {                                                                          \
        return Tsat<FieldType>(p);                                             \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::FieldType>                                                 \
    Foam::ModelType::TsatPrime(const FieldType& p) const                       \
    {                                                                          \
        return TsatPrime<FieldType>(p);                                        \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class saturationTemperatureModel Declaration
\*---------------------------------------------------------------------------*/

class saturationTemperatureModel
{
private:

    // Private Member Functions

        //- Saturation temperature
        template<class FieldType>
        tmp<FieldType> Tsat(const FieldType& p) const;

        //- Saturation temperature derivative w.r.t. pressure
        template<class FieldType>
        tmp<FieldType> TsatPrime(const FieldType& p) const;


public:

    //- Runtime type information
    TypeName("saturationTemperatureModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        saturationTemperatureModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Default construct
        saturationTemperatureModel();


    // Selectors

        //- Select with dictionary
        static autoPtr<saturationTemperatureModel> New
        (
            const dictionary& dict
        );

        //- Select with name within a dictionary
        static autoPtr<saturationTemperatureModel> New
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~saturationTemperatureModel();


    // Member Functions

        //- Saturation temperature for scalarField
        DEFINE_TSAT(scalarField, = 0);

        //- Saturation temperature for volScalarField::Internal
        DEFINE_TSAT(volScalarField::Internal, );

        //- Saturation temperature for volScalarField
        DEFINE_TSAT(volScalarField, );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
