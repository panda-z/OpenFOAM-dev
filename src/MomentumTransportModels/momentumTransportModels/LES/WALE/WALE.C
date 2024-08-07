/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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

#include "WALE.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> WALE<BasicMomentumTransportModel>::Sd
(
    const volTensorField& gradU
) const
{
    return dev(symm(gradU & gradU));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> WALE<BasicMomentumTransportModel>::k
(
    const volTensorField& gradU
) const
{
    volScalarField magSqrSd(magSqr(Sd(gradU)));

    return volScalarField::New
    (
        this->groupName("k"),
        sqr(sqr(Cw_)*this->delta()/Ck_)*
        (
            pow3(magSqrSd)
           /(
               sqr
               (
                   pow(magSqr(symm(gradU)), 5.0/2.0)
                 + pow(magSqrSd, 5.0/4.0)
               )
             + dimensionedScalar
               (
                   "small",
                   dimensionSet(0, 0, -10, 0, 0),
                   small
               )
           )
        )
    );
}


template<class BasicMomentumTransportModel>
void WALE<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = Ck_*this->delta()*sqrt(this->k(fvc::grad(this->U_)));
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
WALE<BasicMomentumTransportModel>::WALE
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    LESeddyViscosity<BasicMomentumTransportModel>
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
    Cw_("Cw", this->coeffDict(), 0.325)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool WALE<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());
        Cw_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> WALE<BasicMomentumTransportModel>::epsilon() const
{
    volScalarField k(this->k(fvc::grad(this->U_)));

    return volScalarField::New
    (
        this->groupName("epsilon"),
        this->Ce_*k*sqrt(k)/this->delta()
    );
}


template<class BasicMomentumTransportModel>
void WALE<BasicMomentumTransportModel>::correct()
{
    LESeddyViscosity<BasicMomentumTransportModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
