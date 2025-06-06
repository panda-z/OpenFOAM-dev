/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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

Global
    MULES

Description
    CMULES: Multidimensional universal limiter for
    explicit corrected implicit solution.

    Solve a convective-only transport equation using an explicit universal
    multi-dimensional limiter to correct an implicit conservative bounded
    obtained using rigorously bounded schemes such as Euler-implicit in time
    upwind in space.

    Parameters are the variable to solve, the normal convective flux and the
    actual explicit flux of the variable which is also used to return limited
    flux used in the bounded-solution.

SourceFiles
    CMULES.C
    CMULESTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef CMULES_H
#define CMULES_H

#include "MULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "gaussConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace MULES
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RdeltaTType, class RhoType, class SpType>
void correct
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp
);

template<class RhoType>
void correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiCorr
);

template<class RhoType, class SpType>
void correct
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp
);

template<class RhoType, class PsiMaxType, class PsiMinType>
void correct
(
    const control& controls,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RhoType,
    class SpType,
    class PsiMaxType,
    class PsiMinType
>
void correct
(
    const control& controls,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class PsiMaxType,
    class PsiMinType
>
void limitCorr
(
    const control& controls,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);

template
<
    class RhoType,
    class SpType,
    class PsiMaxType,
    class PsiMinType
>
void limitCorr
(
    const control& controls,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phiBD,
    surfaceScalarField& phiCorr,
    const SpType& Sp,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace MULES
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CMULESTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
