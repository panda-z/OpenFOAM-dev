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
    Foam::Cloud

Description
    Base cloud calls templated on particle type

SourceFiles
    Cloud.C
    CloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef Cloud_H
#define Cloud_H

#include "cloud.H"
#include "IDLList.H"
#include "IOField.H"
#include "CompactIOField.H"
#include "polyMesh.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParticleType>
class IOPosition;

namespace lagrangian
{
    template<class ParticleType>
    class Cloud;
}

template<class ParticleType>
Ostream& operator<<(Ostream&, const lagrangian::Cloud<ParticleType>&);

namespace lagrangian
{

/*---------------------------------------------------------------------------*\
                            Class Cloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class Cloud
:
    public cloud,
    public IDLList<ParticleType>
{
    // Private Data

        //- Reference to the mesh
        const polyMesh& pMesh_;

        //- Map from patch index to the neighbouring processor index
        labelList patchNbrProc_;

        //- Map from patch index to the corresponding patch index on the
        //  neighbouring processor
        labelList patchNbrProcPatch_;

        //- Map from patch index to connected non-conformal cyclics
        labelListList patchNonConformalCyclicPatches_;

        //- Temporary storage for the global particle positions
        mutable autoPtr<vectorField> globalPositionsPtr_;

        //- Time index
        mutable label timeIndex_;


    // Private Member Functions

        //- Initialise cloud on IO constructor
        void initCloud(const bool checkClass);

        //- Read cloud properties dictionary
        void readCloudUniformProperties();

        //- Write cloud properties dictionary
        void writeCloudUniformProperties() const;

        //- Map from patch index to the neighbouring processor index
        static labelList patchNbrProc(const polyMesh&);

        //- Map from patch index to the corresponding patch index on the
        //  neighbouring processor
        static labelList patchNbrProcPatch(const polyMesh&);

        //- Map from patch index to connected non-conformal cyclics
        static labelListList patchNonConformalCyclicPatches(const polyMesh&);

        //- Store rays necessary for non conformal cyclic transfer
        void storeRays() const;


public:

    friend class particle;
    template<class ParticleT>
    friend class IOPosition;

    typedef ParticleType particleType;

    typedef typename IDLList<ParticleType>::value_type value_type;
    typedef typename IDLList<ParticleType>::iterator iterator;
    typedef typename IDLList<ParticleType>::const_iterator const_iterator;


    // Static data

        //- Name of cloud properties dictionary
        static word cloudPropertiesName;


    // Constructors

        //- Construct from mesh and a list of particles
        Cloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<ParticleType>& particles
        );

        //- Construct from mesh by reading from file with given cloud instance
        //  Optionally disable checking of class name for post-processing
        Cloud
        (
            const polyMesh& pMesh,
            const word& cloudName,
            const bool checkClass = true
        );


    // Member Functions

        // Access

            //- Return the polyMesh reference
            const polyMesh& pMesh() const
            {
                return pMesh_;
            }

            //- Map from patch index to the neighbouring processor index
            const labelList& patchNbrProc() const
            {
                return patchNbrProc_;
            }

            //- Map from patch index to the corresponding patch index on the
            //  neighbouring processor
            const labelList& patchNbrProcPatch() const
            {
                return patchNbrProcPatch_;
            }

            //- Return map from patch index to connected non-conformal cyclics
            const labelListList& patchNonConformalCyclicPatches() const
            {
                return patchNonConformalCyclicPatches_;
            }

            //- Return the number of particles in the cloud
            label size() const
            {
                return IDLList<ParticleType>::size();
            };

            //- Return true to cache per-cell CPU load
            //  Overridden in derived clouds, defaults to false
            bool cpuLoad() const
            {
                return false;
            }


            // Iterators

                const const_iterator begin() const
                {
                    return IDLList<ParticleType>::begin();
                };

                const const_iterator cbegin() const
                {
                    return IDLList<ParticleType>::cbegin();
                };

                const const_iterator end() const
                {
                    return IDLList<ParticleType>::end();
                };

                const const_iterator cend() const
                {
                    return IDLList<ParticleType>::cend();
                };

                iterator begin()
                {
                    return IDLList<ParticleType>::begin();
                };

                iterator end()
                {
                    return IDLList<ParticleType>::end();
                };


        // Edit

            void clear()
            {
                IDLList<ParticleType>::clear();
            };

            //- Transfer particle to cloud
            void addParticle(ParticleType* pPtr);

            //- Remove particle from cloud and delete
            void deleteParticle(ParticleType&);

            //- Remove lost particles from cloud and delete
            void deleteLostParticles();

            //- Reset the particles
            void cloudReset(const Cloud<ParticleType>& c);

            //- Change the particles' state from the end of the previous time
            //  step to the start of the next time step
            void changeTimeStep();

            //- Move the particles
            template<class TrackCloudType>
            void move
            (
                TrackCloudType& cloud,
                typename ParticleType::trackingData& td
            );


        // Mapping

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // Read

            //- Helper to construct IOobject for field and current time.
            IOobject fieldIOobject
            (
                const word& fieldName,
                const IOobject::readOption r
            ) const;

            //- Check lagrangian data field
            template<class DataType>
            void checkFieldIOobject
            (
                const Cloud<ParticleType>& c,
                const IOField<DataType>& data
            ) const;

            //- Check lagrangian data fieldfield
            template<class DataType>
            void checkFieldFieldIOobject
            (
                const Cloud<ParticleType>& c,
                const CompactIOField<Field<DataType>>& data
            ) const;


        // Write

            //- Write the field data for the cloud of particles Dummy at
            //  this level.
            virtual void writeFields() const;

            //- Write using given format, version and compression.
            //  Only writes the cloud file if the Cloud isn't empty
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write = true
            ) const;

            //- Write positions to \<cloudName\>_positions.obj file
            void writePositions() const;

            //- Call this before a topology change. Stores the particles global
            //  positions in the database for use during mapping.
            void storeGlobalPositions() const;


    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const Cloud<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Cloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
