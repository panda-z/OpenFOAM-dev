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
    Foam::ManualInjection

Description
    Manual injection.

    User specifies:
      - Parcel positions in file \c positionsFile
      - Initial parcel velocity

    Properties:
      - Parcel diameters obtained by distribution model
      - All parcels introduced at SOI

SourceFiles
    ManualInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ManualInjection_H
#define ManualInjection_H

#include "InjectionModel.H"
#include "distribution.H"
#include "Switch.H"
#include "GlobalIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ManualInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ManualInjection
:
    public InjectionModel<CloudType>
{
    // Private Data

        //- Name of file containing positions data
        const word positionsFile_;

        //- Field of parcel positions
        GlobalIOField<vector> positions_;

        //- Field of parcel diameters
        scalarList diameters_;

        //- List of coordinates of injector positions
        List<barycentric> injectorCoordinates_;

        //- List of cell labels corresponding to injector positions
        labelList injectorCells_;

        //- List of tetFace labels corresponding to injector positions
        labelList injectorTetFaces_;

        //- List of tetPt labels corresponding to injector positions
        labelList injectorTetPts_;

        //- Mass to inject
        const scalar massTotal_;

        //- Initial parcel velocity
        const vector U0_;

        //- Parcel size distribution model
        const autoPtr<distribution> sizeDistribution_;

        //- Flag to suppress errors if particle injection site is out-of-bounds
        Switch ignoreOutOfBounds_;


public:

    //- Runtime type information
    TypeName("manualInjection");


    // Constructors

        //- Construct from dictionary
        ManualInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ManualInjection(const ManualInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new ManualInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ManualInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void topoChange();

        //- Return the end-of-injection time
        virtual scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual scalar nParcelsToInject(const scalar time0, const scalar time1);

        //- Parcel mass to introduce relative to SOI
        virtual scalar massToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const meshSearch& searchEngine,
                const label parcelI,
                const label nParcels,
                const scalar time,
                barycentric& coordinates,
                label& celli,
                label& tetFacei,
                label& tetPti,
                label& facei
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType::trackingData& td,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ManualInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
