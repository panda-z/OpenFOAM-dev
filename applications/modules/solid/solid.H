/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    Foam::solvers::solid

Description
    Solver module for thermal transport in solid domains and regions for
    conjugate heat transfer, HVAC and similar applications, with optional mesh
    motion and mesh topology changes.

SourceFiles
    solid.C

\*---------------------------------------------------------------------------*/

#ifndef solid_H
#define solid_H

#include "solver.H"
#include "solidThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class solid Declaration
\*---------------------------------------------------------------------------*/

class solid
:
    public solver
{

protected:

    // Control parameters

        scalar maxDeltaT_;


    // Thermophysical properties

        autoPtr<solidThermo> thermoPtr_;

        solidThermo& thermo_;

        volScalarField& T_;


    // Thermophysical transport

        //- Pointer to the solid thermophysical transport model
        autoPtr<solidThermophysicalTransportModel> thermophysicalTransport;


    // Protected Member Functions

        //- Return true if the solver's dependencies have been modified
        virtual bool dependenciesModified() const;

        //- Read controls
        virtual bool read();


public:

    // Public Data

        //- Reference to the solid thermophysical properties
        const solidThermo& thermo;

        //- Reference to the temperature field
        const volScalarField& T;


    //- Runtime type information
    TypeName("solid");


    // Constructors

        //- Construct from region mesh and thermophysical properties
        solid(fvMesh& mesh, autoPtr<solidThermo>);

        //- Construct from region mesh
        solid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        solid(const solid&) = delete;


    //- Destructor
    virtual ~solid();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the beginning of the PIMPLE loop
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum transport
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport
        virtual void thermophysicalTransportCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
