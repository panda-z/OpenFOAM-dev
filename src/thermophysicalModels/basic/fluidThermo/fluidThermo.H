/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::fluidThermo

Description
    Base-class for fluid thermodynamic properties.

See also
    Foam::basicThermo

SourceFiles
    fluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidThermo_H
#define fluidThermo_H

#include "basicThermo.H"
#include "viscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fluidThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo
:
    virtual public basicThermo,
    public viscosity
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("fluidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<fluidThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~fluidThermo();


    // Member Functions

        // Thermodynamic state

            //- Pressure [Pa]
            virtual const volScalarField& p() const = 0;

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p() = 0;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Derived thermodynamic properties

            //- Rename the thermodynamic density field if stored and return
            //  This is used by solvers which create a separate continuity rho
            //  [kg/m^3]
            virtual tmp<volScalarField> renameRho() = 0;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho) = 0;


        // Transport state

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const = 0;


        // Derived transport properties

            //- Kinematic viscosity of mixture [m^2/s]
            tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            tmp<scalarField> nu(const label patchi) const;
};


/*---------------------------------------------------------------------------*\
                 Class fluidThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo::implementation
:
    virtual public fluidThermo
{
protected:

    // Protected data

        // Fields

            //- Pressure [Pa]
            volScalarField& p_;

            //- Compressibility [s^2/m^2]
            volScalarField psi_;

            //- Dynamic viscosity [kg/m/s]
            volScalarField mu_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        implementation(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Thermodynamic state

            //- Pressure [Pa]
            virtual const volScalarField& p() const;

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;


        // Transport state

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
