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
    Foam::CollidingCloud

Description
    Adds collisions to clouds

SourceFiles
    CollidingCloudI.H
    CollidingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef CollidingCloud_H
#define CollidingCloud_H

#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "dimensionedTypes.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CollisionModel;


/*---------------------------------------------------------------------------*\
                     Class CollidingCloudName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(CollidingCloud);


/*---------------------------------------------------------------------------*\
                       Class CollidingCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CollidingCloud
:
    public CloudType,
    public CollidingCloudName
{
public:

    // Public Typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef CollidingCloud<CloudType> collidingCloudType;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<CollidingCloud<CloudType>> cloudCopyPtr_;


protected:

    // Protected data

        //- Thermo parcel constant properties
        typename parcelType::constantProperties constProps_;


        // References to the cloud sub-models

            //- Collision model
            autoPtr<CollisionModel<CollidingCloud<CloudType>>>
                collisionModel_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(CollidingCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier fields
        CollidingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            const bool readFields = true
        );

        //- Construct given carrier fields and thermo
        CollidingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const fluidThermo& carrierThermo,
            const bool readFields = true
        );

        //- Copy constructor with new name
        CollidingCloud
        (
            CollidingCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        CollidingCloud
        (
            const fvMesh& mesh,
            const word& name,
            const CollidingCloud<CloudType>& c
        );

        //- Disallow default bitwise copy construction
        CollidingCloud(const CollidingCloud&) = delete;

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<lagrangian::Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<lagrangian::Cloud<parcelType>>
            (
                new CollidingCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<lagrangian::Cloud<parcelType>> cloneBare
        (
            const word& name
        ) const
        {
            return autoPtr<lagrangian::Cloud<parcelType>>
            (
                new CollidingCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~CollidingCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const CollidingCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;


            // Sub-models

                //- Return const access to the collision model
                inline const CollisionModel<CollidingCloud<CloudType>>&
                    collision() const;

                //- Return reference to the collision model
                inline CollisionModel<CollidingCloud<CloudType>>&
                    collision();

        // Check

            //- Total rotational kinetic energy in the system
            inline scalar rotationalKineticEnergyOfSystem() const;


        // Cloud evolution functions

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the cloud
            void evolve();

            //- Particle motion
            template<class TrackCloudType>
            void motion
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );


        // I-O

            //- Print cloud information
            void info();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CollidingCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CollidingCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CollidingCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
