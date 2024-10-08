/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
LESeddyViscosity<BasicMomentumTransportModel>::LESeddyViscosity
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    eddyViscosity<LESModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    Ck_("Ck", this->coeffDict(), 0.094),
    Ce_("Ce", this->coeffDict(), 1.048)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool LESeddyViscosity<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<LESModel<BasicMomentumTransportModel>>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
        Ce_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
LESeddyViscosity<BasicMomentumTransportModel>::epsilon() const
{
    tmp<volScalarField> tk(this->k());

    return volScalarField::New
    (
        this->groupName("epsilon"),
        Ce_*tk()*sqrt(tk())/this->delta()
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
LESeddyViscosity<BasicMomentumTransportModel>::omega() const
{
    tmp<volScalarField> tk(this->k());

    return volScalarField::New
    (
        this->groupName("omega"),
        (Ce_/Ck_)*sqrt(tk())/this->delta()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
