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

#include "leanInhomogeneousMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::leanInhomogeneousMixture<ThermoType>::leanInhomogeneousMixture
(
    const dictionary& dict
)
:
    stoicRatio_("stoichiometricAirFuelMassRatio", dimless, dict),
    fuel_("fuel", dict.subDict("fuel")),
    oxidant_("oxidant", dict.subDict("oxidant")),
    products_("burntProducts", dict.subDict("burntProducts")),
    mixture_("mixture", fuel_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::leanInhomogeneousMixture<ThermoType>::fres
(
    const scalar ft
) const
{
     return max(ft - (scalar(1) - ft)/stoicRatio_.value(), scalar(0));
}


template<class ThermoType>
Foam::scalar Foam::leanInhomogeneousMixture<ThermoType>::fres
(
    const scalarFieldListSlice& Y
) const
{
    return fres(Y[FT]);
}


template<class ThermoType>
const ThermoType& Foam::leanInhomogeneousMixture<ThermoType>::mixture
(
    const scalar ft,
    const scalar b
) const
{
    if (ft < 0.0001)
    {
        return oxidant_;
    }
    else
    {
        const scalar fu = b*ft + (1 - b)*fres(ft);
        const scalar ox = 1 - ft - (ft - fu)*stoicRatio_.value();
        const scalar pr = 1 - fu - ox;

        mixture_ = fu*fuel_;
        mixture_ += ox*oxidant_;
        mixture_ += pr*products_;

        return mixture_;
    }
}


template<class ThermoType>
const typename Foam::leanInhomogeneousMixture<ThermoType>::thermoMixtureType&
Foam::leanInhomogeneousMixture<ThermoType>::thermoMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT], Y[B]);
}


template<class ThermoType>
const typename Foam::leanInhomogeneousMixture<ThermoType>::transportMixtureType&
Foam::leanInhomogeneousMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT], Y[B]);
}


template<class ThermoType>
const typename Foam::leanInhomogeneousMixture<ThermoType>::transportMixtureType&
Foam::leanInhomogeneousMixture<ThermoType>::transportMixture
(
    const scalarFieldListSlice&,
    const thermoMixtureType& mixture
) const
{
    return mixture;
}


template<class ThermoType>
const typename Foam::leanInhomogeneousMixture<ThermoType>::thermoType&
Foam::leanInhomogeneousMixture<ThermoType>::reactants
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT], 1);
}


template<class ThermoType>
const typename Foam::leanInhomogeneousMixture<ThermoType>::thermoType&
Foam::leanInhomogeneousMixture<ThermoType>::products
(
    const scalarFieldListSlice& Y
) const
{
    return mixture(Y[FT], 0);
}


template<class ThermoType>
void Foam::leanInhomogeneousMixture<ThermoType>::read(const dictionary& dict)
{
    stoicRatio_ =
        dimensionedScalar("stoichiometricAirFuelMassRatio", dimless, dict);

    fuel_ = ThermoType("fuel", dict.subDict("fuel"));
    oxidant_ = ThermoType("oxidant", dict.subDict("oxidant"));
    products_ = ThermoType("burntProducts", dict.subDict("burntProducts"));
}


// ************************************************************************* //
