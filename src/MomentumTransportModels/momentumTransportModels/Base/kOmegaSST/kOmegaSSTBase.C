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

#include "kOmegaSSTBase.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
void
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::boundOmega()
{
    omega_ = max(omega_, k_/(this->nutMaxCoeff_*this->nu()));
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*this->y()),
                scalar(500)*this->nu()/(sqr(this->y())*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(this->y()))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*this->y()),
            scalar(500)*this->nu()/(sqr(this->y())*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*this->nu()/(omega_*sqr(this->y())),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& F2
)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F2*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::
correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))), F23());
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return betaStar_*omega_();
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::
omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::kOmegaSST
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
    MomentumTransportModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    alphaK1_("alphaK1", this->coeffDict(), 0.85),
    alphaK2_("alphaK2", this->coeffDict(), 1.0),
    alphaOmega1_("alphaOmega1", this->coeffDict(), 0.5),
    alphaOmega2_("alphaOmega2", this->coeffDict(), 0.856),
    gamma1_("gamma1", this->coeffDict(), 5.0/9.0),
    gamma2_("gamma2", this->coeffDict(), 0.44),
    beta1_("beta1", this->coeffDict(), 0.075),
    beta2_("beta2", this->coeffDict(), 0.0828),
    betaStar_("betaStar", this->coeffDict(), 0.09),
    a1_("a1", this->coeffDict(), 0.31),
    b1_("b1", this->coeffDict(), 1.0),
    c1_("c1", this->coeffDict(), 10.0),
    F3_(this->coeffDict().template lookupOrDefault<Switch>("F3", false)),

    k_
    (
        IOobject
        (
            this->groupName("k"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            this->groupName("omega"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    boundOmega();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MomentumTransportModel, class BasicMomentumTransportModel>
bool kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::read()
{
    if (MomentumTransportModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class MomentumTransportModel, class BasicMomentumTransportModel>
void kOmegaSST<MomentumTransportModel, BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    MomentumTransportModel::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu(dev(twoSymm(tgradU()())) && tgradU()());
    volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma
           *min
            (
                GbyNu,
                (c1_/a1_)*betaStar_*omega_()
               *max(a1_*omega_(), b1_*F23()*sqrt(S2()))
            )
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvModels.source(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvConstraints.constrain(omega_);
        boundOmega();
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, F23), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);
    boundOmega();

    correctNut(S2, F23);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
