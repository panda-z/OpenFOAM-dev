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

\*---------------------------------------------------------------------------*/

#include "SpalartAllmarasIDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::IDDESalpha() const
{
    return volScalarField::Internal::New
    (
        typedName("alpha"),
        max(0.25 - this->y()()/IDDESDelta_.hmax(), scalar(-5))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::ft
(
    const volScalarField::Internal& magGradU
) const
{
    return volScalarField::Internal::New
    (
        typedName("ft"),
        tanh(pow3(sqr(ct_)*rd(this->nut_, magGradU)))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::fl
(
    const volScalarField::Internal& magGradU
) const
{
    return volScalarField::Internal::New
    (
        typedName("fl"),
        tanh(pow(sqr(cl_)*rd(this->nu(), magGradU), 10))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::rd
(
    const volScalarField::Internal& nur,
    const volScalarField::Internal& magGradU
) const
{
    return volScalarField::Internal::New
    (
        typedName("rd"),
        min
        (
            nur
           /(
               max
               (
                   magGradU,
                   dimensionedScalar(magGradU.dimensions(), small)
               )*sqr(this->kappa_*this->y()())
            ),
            scalar(10)
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::fd
(
    const volScalarField::Internal& magGradU
) const
{
    return volScalarField::Internal::New
    (
        typedName("fd"),
        1 - tanh(pow3(8*rd(this->nuEff(), magGradU)))
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::dTilda
(
    const volScalarField::Internal& chi,
    const volScalarField::Internal& fv1,
    const volTensorField::Internal& gradU
) const
{
    const volScalarField::Internal alpha(IDDESalpha());

    const volScalarField::Internal expTerm
    (
        typedName("expTerm"),
        exp(sqr(alpha))
    );

    const volScalarField::Internal magGradU(typedName("magGradU"), mag(gradU));

    tmp<volScalarField::Internal> fHill
    (
        volScalarField::Internal::New
        (
            typedName("fHill"),
            2*(pos0(alpha)*pow(expTerm, -11.09) + neg(alpha)*pow(expTerm, -9.0))
        )
    );

    tmp<volScalarField::Internal> fStep
    (
        volScalarField::Internal::New
        (
            typedName("fStep"),
            min(2*pow(expTerm, -9.0), scalar(1))
        )
    );

    const volScalarField::Internal fHyb
    (
        typedName("fHyb"),
        max(1 - fd(magGradU), fStep)
    );

    tmp<volScalarField::Internal> fAmp
    (
        volScalarField::Internal::New
        (
            typedName("fAmp"),
            1 - max(ft(magGradU), fl(magGradU))
        )
    );

    tmp<volScalarField::Internal> fRestore
    (
        volScalarField::Internal::New
        (
            typedName("fRestore"),
            max(fHill - 1, scalar(0))*fAmp
        )
    );

    // IGNORING ft2 terms
    const volScalarField::Internal Psi
    (
        typedName("Psi"),
        sqrt
        (
            min
            (
                scalar(100),
                (
                    1
                  - this->Cb1_*this->fv2(chi, fv1)
                   /(this->Cw1_*sqr(this->kappa_)*fwStar_)
                 )/max(small, fv1)
            )
        )
    );

    return volScalarField::Internal::New
    (
        typedName("dTilda"),
        max
        (
            dimensionedScalar(dimLength, small),
            fHyb*(1 + fRestore*Psi)*this->y()
          + (1 - fHyb)*this->CDES_*Psi*this->delta()
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
SpalartAllmarasIDDES<BasicMomentumTransportModel>::SpalartAllmarasIDDES
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
    SpalartAllmarasDES<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    fwStar_("fwStar", this->coeffDict(), 0.424),
    cl_("cl", this->coeffDict(), 3.55),
    ct_("ct", this->coeffDict(), 1.63),
    IDDESDelta_(refCast<IDDESDelta>(this->delta_()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool SpalartAllmarasIDDES<BasicMomentumTransportModel>::read()
{
    if (SpalartAllmarasDES<BasicMomentumTransportModel>::read())
    {
        fwStar_.readIfPresent(this->coeffDict());
        cl_.readIfPresent(this->coeffDict());
        ct_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
