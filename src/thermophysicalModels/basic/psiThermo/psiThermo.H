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
    Foam::psiThermo

Description
    Base-class for fluid thermodynamic properties based on compressibility.

See also
    Foam::basicThermo

SourceFiles
    psiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiThermo_H
#define psiThermo_H

#include "PsiThermo.H"
#include "pureThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class psiThermo Declaration
\*---------------------------------------------------------------------------*/

class psiThermo
:
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    // Public Typedefs

        //- The derived type
        template<class MixtureType>
        using DerivedThermoType =
            PsiThermo<BasicThermo<MixtureType, composite>>;

        //- The derived name
        static const word derivedThermoName;


    //- Runtime type information
    TypeName("psiThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<psiThermo> New
        (
            const fvMesh& mesh,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~psiThermo();


    // Member Functions

        // Derived thermodynamic properties

            //- Density [kg/m^3] - uses current value of pressure
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Return the thermodynamic density field [kg/m^3]
            //  This is used by solvers which create a separate continuity rho
            virtual tmp<volScalarField> renameRho();

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution.
            //  For psiThermo does nothing.
            virtual void correctRho(const volScalarField& deltaRho);
};


/*---------------------------------------------------------------------------*\
                  Class psiThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class psiThermo::implementation
:
    virtual public psiThermo
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        implementation(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Derived thermodynamic properties

            //- Density [kg/m^3] - uses current value of pressure
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


/*---------------------------------------------------------------------------*\
                    Class psiThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class psiThermo::composite
:
    public basicThermo::implementation,
    public pureThermo,
    public fluidThermo::implementation,
    public psiThermo::implementation
{
public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        composite
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& phaseName
        )
        :
            basicThermo::implementation(dict, mesh, phaseName),
            fluidThermo::implementation(dict, mesh, phaseName),
            psiThermo::implementation(dict, mesh, phaseName)
        {}

        //- Construct from dictionary, mesh and phase name
        template<class MixtureType>
        composite
        (
            const dictionary& dict,
            const MixtureType& mixture,
            const fvMesh& mesh,
            const word& phaseName
        )
        :
            basicThermo::implementation(dict, mesh, phaseName),
            fluidThermo::implementation(dict, mesh, phaseName),
            psiThermo::implementation(dict, mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
