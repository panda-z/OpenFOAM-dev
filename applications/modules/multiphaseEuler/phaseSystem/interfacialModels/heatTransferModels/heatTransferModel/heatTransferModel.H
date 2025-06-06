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

Class
    Foam::heatTransferModel

Description
    Model for heat transfer between phases

SourceFiles
    heatTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferModel_H
#define heatTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "BlendedInterfacialModel.H"
#include "SidedInterfacialModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class heatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;


public:

    //- Runtime type information
    TypeName("heatTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface,
                const bool registerObject
            ),
            (dict, interface, registerObject)
        );


    // Static Data Members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        //- Construct from a dictionary and an interface
        heatTransferModel
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~heatTransferModel();


    // Selectors

        static autoPtr<heatTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true,
            const bool registerObject=true
        );


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        tmp<volScalarField> K() const;

        //- The heat transfer function K used in the enthalpy equation
        //    ddt(alpha1*rho1*ha) + ... = ... K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ... = ... K*(Tb - Ta)
        //  with a specified residual volume fraction
        virtual tmp<volScalarField> K(const scalar residualAlpha) const = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                   Class blendedHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class blendedHeatTransferModel
:
    public BlendedInterfacialModel<heatTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<heatTransferModel>::
            BlendedInterfacialModel;


    // Selectors

        inline static autoPtr<blendedHeatTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const dictionary& blendingDict
        )
        {
            return autoPtr<blendedHeatTransferModel>
            (
                new blendedHeatTransferModel(dict, interface, blendingDict)
            );
        }


    // Member Functions

        //- Return the heat transfer coefficient
        tmp<volScalarField> K() const;

        //- Return the heat transfer coefficient
        tmp<volScalarField> K(const scalar residualAlpha) const;
};


/*---------------------------------------------------------------------------*\
                 Class sidedBlendedHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class sidedBlendedHeatTransferModel
:
    public SidedInterfacialModel<blendedHeatTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            SidedInterfacialModel<blendedHeatTransferModel>::
            SidedInterfacialModel;


    // Selectors

        inline static autoPtr<sidedBlendedHeatTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const dictionary& blendingDict
        )
        {
            return autoPtr<sidedBlendedHeatTransferModel>
            (
                new sidedBlendedHeatTransferModel(dict, interface, blendingDict)
            );
        }
};


/*---------------------------------------------------------------------------*\
                   Class sidedHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class sidedHeatTransferModel
:
    public SidedInterfacialModel<heatTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            SidedInterfacialModel<heatTransferModel>::
            SidedInterfacialModel;


    // Selectors

        inline static autoPtr<sidedHeatTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        )
        {
            return autoPtr<sidedHeatTransferModel>
            (
                new sidedHeatTransferModel(dict, interface)
            );
        }


    // Member Functions

        //- Return the heat transfer coefficient
        tmp<volScalarField> KinThe(const phaseModel& phase) const;

        //- Return the heat transfer coefficient
        tmp<volScalarField> KinThe
        (
            const phaseModel& phase,
            const scalar residualAlpha
        ) const;
};


/*---------------------------------------------------------------------------*\
                 Class blendedSidedHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class blendedSidedHeatTransferModel
:
    public BlendedInterfacialModel<sidedHeatTransferModel>
{
public:

    // Constructors

        //- Inherit base class constructors
        using
            BlendedInterfacialModel<sidedHeatTransferModel>::
            BlendedInterfacialModel;


    // Selectors

        inline static autoPtr<blendedSidedHeatTransferModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const dictionary& blendingDict
        )
        {
            return autoPtr<blendedSidedHeatTransferModel>
            (
                new blendedSidedHeatTransferModel(dict, interface, blendingDict)
            );
        }


    // Member Functions

        //- Return the heat transfer coefficient
        tmp<volScalarField> KinThe(const phaseModel& phase) const;

        //- Return the heat transfer coefficient
        tmp<volScalarField> KinThe
        (
            const phaseModel& phase,
            const scalar residualAlpha
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
