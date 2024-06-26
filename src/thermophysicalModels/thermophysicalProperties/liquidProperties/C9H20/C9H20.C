/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "C9H20.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C9H20, 0);
    addToRunTimeSelectionTable(liquidProperties, C9H20,);
    addToRunTimeSelectionTable(liquidProperties, C9H20, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C9H20::C9H20()
:
    liquidProperties
    (
        typeName,
        128.258,
        594.60,
        2.29e+6,
        0.544,
        0.252,
        219.66,
        4.3058e-1,
        423.97,
        0.0,
        0.4435,
        1.56e+4
    ),
    rho_("rho", 62.06019846, 0.26147, 594.6, 0.28281),
    pv_("pv", 109.35, -90304.0, -12.882, 7.8544e-06, 2.0),
    hl_("hl", 594.60, 470691.886665939, 0.38522, 0.0, 0.0, 0.0),
    Cp_
    (
        "Cp",
        2986.79224687739,
       -8.88677509395125,
        0.0211300659607978,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -2825628.50868792,
        2986.79224687739,
       -4.44338754697563,
        0.00704335532026592,
        0.0,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        1183.16206396482,
        3832.11963386299,
        1644.8,
        2705.48425829188,
        749.6
    ),
    B_
    (
        "B",
        0.00304542406711472,
       -3.65357326638494,
       -520825.211682702,
       -6.15400208953827e+18,
        1.41901479829718e+21
    ),
    mu_("mu", -21.149, 1658, 1.454, 0.0, 0.0),
    mug_("mug", 1.0344e-07, 0.77301, 220.47, 0.0),
    kappa_("kappa", 0.209, -0.000264, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", -0.065771, 0.27198, -3482.3, -1580300.0),
    sigma_("sigma", 594.60, 0.054975, 1.2897, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 128.258, 28.0), // note: Same as nHeptane
    hf_(h_.value(Tstd))
{}


Foam::C9H20::C9H20
(
    const liquidProperties& l,
    const Function1s::NSRDS5& density,
    const Function1s::NSRDS1& vapourPressure,
    const Function1s::NSRDS6& heatOfVapourisation,
    const Function1s::NSRDS0& heatCapacity,
    const Function1s::NSRDS0& enthalpy,
    const Function1s::NSRDS7& idealGasHeatCapacity,
    const Function1s::NSRDS4& secondVirialCoeff,
    const Function1s::NSRDS1& dynamicViscosity,
    const Function1s::NSRDS2& vapourDynamicViscosity,
    const Function1s::NSRDS0& thermalConductivity,
    const Function1s::NSRDS2& vapourThermalConductivity,
    const Function1s::NSRDS6& surfaceTension,
    const Function2s::APIdiffCoef& vapourDiffusivity
)
:
    liquidProperties(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    Cp_(heatCapacity),
    h_(enthalpy),
    Cpg_(idealGasHeatCapacity),
    B_(secondVirialCoeff),
    mu_(dynamicViscosity),
    mug_(vapourDynamicViscosity),
    kappa_(thermalConductivity),
    kappag_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffusivity),
    hf_(h_.value(Tstd))
{}


Foam::C9H20::C9H20(const dictionary& dict)
:
    C9H20()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C9H20::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
