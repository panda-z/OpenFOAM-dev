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
    Foam::functionObjects::streamlines

Description
    Generates streamline data by sampling a set of user-specified fields along
    a particle track, transported by a user-specified velocity field.

    Example of function object specification:
    \verbatim
    streamlines1
    {
        type            streamlines;
        libs            ("libfieldFunctionObjects.so");

        writeControl    writeTime;

        setFormat       vtk;
        U               U;
        direction       both;

        fields
        (
            U
            p
        );

        lifeTime        10000;
        trackLength     1e-3;
        nSubCycle       5;

        seedSampleSet
        {
            type        lineUniform;
            axis        xyz;
            start       (-0.0205 0.0001 0.00001);
            end         (-0.0205 0.0005 0.00001);
            nPoints     100;
        }
    }
    \endverbatim

Usage
    \table
        Property      | Description                | Required   | Default value
        type          | Type name: streamlines     | yes        |
        setFormat     | Output data type           | yes        |
        U             | Tracking velocity field name | no       | U
        direction     | Direction in which to track | yes       |
        outside       | Track outside of periodic meshes | no   | no
        fields        | Fields to sample           | yes        |
        writeTime     | Write the flow time along the streamlines | no | no
        lifetime      | Maximum number of particle tracking steps | yes |
        trackLength   | Tracking segment length    | no         |
        nSubCycle     | Number of tracking steps per cell | no  |
        cloudName     | Cloud name to use          | no         |
        seedSampleSet | Seeding method (see below) | yes        |
    \endtable

    Where the \c seedSampleSet \c type is typically one of
    \plaintable
        lineUniform    | uniform particle seeding along a line
        sphereRandom   | random particle seeding within a sphere
        boundaryRandom | random particle seeding on a number of patches
        points         | a specified set of locations
    \endplaintable

    Note:
        When specifying the track resolution, the \c trackLength OR \c nSubCycle
        option should be used

See also
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::sampledSet

SourceFiles
    streamlines.C

\*---------------------------------------------------------------------------*/

#ifndef streamlines_functionObject_H
#define streamlines_functionObject_H

#include "fvMeshFunctionObject.H"
#include "setWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class sampledSet;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class streamlines Declaration
\*---------------------------------------------------------------------------*/

class streamlines
:
    public fvMeshFunctionObject
{
public:

    // Public data types

        //- Track direction enumerations
        enum class trackDirection
        {
            forward,
            backward,
            both
        };

        //- Track direction enumeration names
        static const NamedEnum<trackDirection, 3> trackDirectionNames_;


private:

    // Private Data

        //- List of fields to sample
        wordList fields_;

        //- Field to transport particle with
        word UName_;

        //- Interpolation scheme to use
        word interpolationScheme_;

        //- The direction in which to track
        trackDirection trackDirection_;

        //- Whether or not to track outside of the mesh in periodic geometries
        Switch trackOutside_;

        //- Maximum lifetime (= number of cells) of particle
        label lifeTime_;

        //- Number of subcycling steps
        label nSubCycle_;

        //- Track length
        scalar trackLength_;

        //- Optional specified name of particles
        word cloudName_;

        //- Write the streamlines ages
        Switch writeAge_;

        //- Seed set engine
        autoPtr<sampledSet> sampledSetPtr_;

        //- File writer
        autoPtr<setWriter> formatterPtr_;


public:

    //- Runtime type information
    TypeName("streamlines");


    // Constructors

        //- Construct from Time and dictionary
        streamlines
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        streamlines(const streamlines&) = delete;


    //- Destructor
    virtual ~streamlines();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Do nothing
        virtual bool execute();

        //- Calculate and write the streamlines
        virtual bool write();

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const streamlines&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
