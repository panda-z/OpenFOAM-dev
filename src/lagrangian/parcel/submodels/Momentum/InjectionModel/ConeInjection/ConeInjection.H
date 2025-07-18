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
    Foam::ConeInjection

Description
    This injector injects particles in a number of cones. The user specifies a
    position and a direction to inject at, and two angles to inject between.
    Optionally, this injector can introduce particles over a disc, instead of
    at a point, in which case inner and outer diameters of the disc are also
    specified.

    The velocity is specified either as constant, or it is calculated from an
    injection pressure, or it is calculated from the injector mass flow rate
    and a discharge coefficient; i.e.:

    Constant velocity:
    \f[
    U = U_{constant}
    \f]

    Pressure driven velocity:
    \f[
    U = \sqrt{2(p_{injection} - p)/\rho}
    \f]

    Flow rate and discharge:
    \f[
    U = \dot{m}/(\rho A C_{discharge})
    \f]

Usage
    \table
    Property        | Description                    | Required | Default
    position        | The injection position         | yes      |
    direction       | The injection direction        | yes      |
    thetaInner      | The inner cone angle           | yes      |
    thetaOuter      | The outer cone angle           | yes      |
    injectionMethod | Inject at a point or on a disc | no       | point
    dInner          | The inner disc diameter        |\\
                                               if disc or flowRateAndDischarge |
    dInner          | The outer disc diameter        |\\
                                               if disc or flowRateAndDischarge |
    flowType        | Inject with constantVelocity, pressureDrivenVelocity \\
                                 or flowRateAndDischarge | no | constantVelocity
    Umag            | The injection velocity         | if constantVelocity |
    Pinj            | The injection pressure         |\\
                                                     if pressureDrivenVelocity |
    Cd              | The discharge coefficient      | if flowRateAndDischarge |
    \endtable

    Example specification:

    \verbatim
    injectionModels
    {
        model1
        {
            type            coneInjection;

            // Times
            SOI             0;
            duration        1;

            // Quantities
            parcelsPerSecond 1000000;

            // - Inject parcels with a fixed number of particles
            nParticle       1;

            //// - Or, inject a specified mass of particles
            //massTotal       6.0e-6;
            //uniformParcelSize volume;

            // Sizes
            sizeDistribution
            {
                type        fixedValue;
                fixedValueDistribution
                {
                    value   0.0025;
                }
            }

            // Geometry
            position        (-0.15 -0.1 0);
            direction       (1 0 0);
            thetaInner      0;
            thetaOuter      45;

            // - Inject at a point
            injectionMethod point;

            //// - Or, inject over a disc:
            //injectionMethod disc;
            //dInner          0;
            //dOuter          0.05;

            // Velocity

            // - Inject with constant velocity
            flowType        constantVelocity;
            Umag            1;

            //// - Or, inject with flow rate and discharge coefficient
            ////   This also requires massTotal, dInner and dOuter
            //flowType        flowRateAndDischarge;
            //Cd              0.9;

            //// - Or, inject at a pressure
            //flowType        pressureDrivenVelocity;
            //Pinj            10e5;
        }
    }
    \endverbatim

SourceFiles
    ConeInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ConeInjection_H
#define ConeInjection_H

#include "InjectionModel.H"
#include "distribution.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ConeInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConeInjection
:
    public InjectionModel<CloudType>
{
public:

    //- Injection method enumeration
    enum injectionMethod
    {
        imPoint,
        imDisc
    };

    //- Flow type enumeration
    enum flowType
    {
        ftConstantVelocity,
        ftPressureDrivenVelocity,
        ftFlowRateAndDischarge
    };


private:

    // Private Data

        //- Point/disc injection method
        injectionMethod injectionMethod_;

        //- Flow type
        flowType flowType_;

        //- Position of the injector
        const autoPtr<Function1<vector>> position_;

        //- Centreline direction in which to inject
        const autoPtr<Function1<vector>> direction_;

        //- Coordinates corresponding to the injector position
        barycentric injectorCoordinates_;

        //- Cell label corresponding to the injector position
        label injectorCell_;

        //- Tet-face label corresponding to the injector position
        label injectorTetFace_;

        //- Tet-point label corresponding to the injector position
        label injectorTetPt_;

        //- Injection duration [s]
        const scalar duration_;

        //- Mass flow rate relative to SOI []
        const autoPtr<Function1<scalar>> massFlowRate_;

        //- Number of parcels to introduce per second
        const autoPtr<Function1<scalar>> parcelsPerSecond_;

        //- Inner half-cone angle relative to SOI [deg]
        const autoPtr<Function1<scalar>> thetaInner_;

        //- Outer half-cone angle relative to SOI [deg]
        const autoPtr<Function1<scalar>> thetaOuter_;

        //- Parcel size distribution model
        const autoPtr<distribution> sizeDistribution_;


        // Disc geometry

            //- The inner disc diameter [m]
            scalar dInner_;

            //- The outer disc diameter [m]
            scalar dOuter_;


        // Velocity model coefficients

            //- Parcel velocity [m/s]
            autoPtr<Function1<scalar>> Umag_;

            //- Discharge coefficient []
            autoPtr<Function1<scalar>> Cd_;

            //- Injection pressure [Pa]
            autoPtr<Function1<scalar>> Pinj_;


    // Private Member Functions

        //- Set the injection type
        void setInjectionMethod();

        //- Set the injection flow type
        void setFlowType();


public:

    //- Runtime type information
    TypeName("coneInjection");


    // Constructors

        //- Construct from dictionary
        ConeInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ConeInjection(const ConeInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new ConeInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ConeInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void topoChange();

        //- Return the end-of-injection time
        scalar timeEnd() const;

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
    #include "ConeInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
