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

#include "uniform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(uniform, 0);
    addToRunTimeSelectionTable(distribution, uniform, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
auto Foam::distributions::uniform::Phi(const Type& x, const label q)
{
    if (q == -1)
    {
        return log(x);
    }
    else
    {
        return integerPow(x, 1 + q)/(1 + q);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::uniform::uniform
(
    const unitConversion& units,
    const dictionary& dict,
    const label sampleQ,
    randomGenerator&& rndGen
)
:
    FieldDistribution<distribution, uniform>
    (
        typeName,
        units,
        dict,
        sampleQ,
        std::move(rndGen)
    ),
    min_(dict.lookupBackwardsCompatible<scalar>({"min", "minValue"}, units)),
    max_(dict.lookupBackwardsCompatible<scalar>({"max", "maxValue"}, units)),
    Phi0_(Phi(min_, q())),
    Phi1_(Phi(max_, q()))
{
    validateBounds(dict);
    if (q() != 0) validatePositive(dict);
}


Foam::distributions::uniform::uniform(const uniform& d, const label sampleQ)
:
    FieldDistribution<distribution, uniform>(d, sampleQ),
    min_(d.min_),
    max_(d.max_),
    Phi0_(Phi(min_, q())),
    Phi1_(Phi(max_, q()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::uniform::~uniform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::uniform::sample() const
{
    const scalar s = rndGen_.sample01<scalar>();

    if (q() == -1)
    {
        return min_*pow(max_/min_, s);
    }
    else
    {
        const scalar PhiS = (1 - s)*Phi0_ + s*Phi1_;
        return integerRoot((1 + q())*PhiS, 1 + q());
    }
}


Foam::scalar Foam::distributions::uniform::min() const
{
    return min_;
}


Foam::scalar Foam::distributions::uniform::max() const
{
    return max_;
}


Foam::scalar Foam::distributions::uniform::mean() const
{
    return (Phi(max_, q() + 1) - Phi(min_, q() + 1))/(Phi1_ - Phi0_);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::uniform::integralPDFxPow
(
    const scalarField& x,
    const label e,
    const bool
) const
{
    const scalarField xClip(Foam::min(Foam::max(x, min()), max()));
    return (Phi(xClip, q() + e) - Phi(min_, q() + e))/(Phi1_ - Phi0_);
}


void Foam::distributions::uniform::write
(
    Ostream& os,
    const unitConversion& units
) const
{
    distribution::write(os, units);

    writeEntry(os, "min", units, min_);
    writeEntry(os, "max", units, max_);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::uniform::plotPDF(const scalarField& x) const
{
    if (q() == -1)
    {
        return clipPDF(x, 1/x/log(max_/min_));
    }
    else
    {
        return clipPDF(x, integerPow(x, q())/(Phi1_ - Phi0_));
    }
}


// ************************************************************************* //
