/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::logPolynomialTransport<Thermo, PolySize>::logPolynomialTransport
(
    const Thermo& t,
    const Polynomial<PolySize>& muLogCoeffs,
    const Polynomial<PolySize>& kappaLogCoeffs
)
:
    Thermo(t),
    muLogCoeffs_(muLogCoeffs),
    kappaLogCoeffs_(kappaLogCoeffs)
{}


template<class Thermo, int PolySize>
inline Foam::logPolynomialTransport<Thermo, PolySize>::logPolynomialTransport
(
    const word& name,
    const logPolynomialTransport& pt
)
:
    Thermo(name, pt),
    muLogCoeffs_(pt.muLogCoeffs_),
    kappaLogCoeffs_(pt.kappaLogCoeffs_)
{}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::logPolynomialTransport<Thermo, PolySize>>
Foam::logPolynomialTransport<Thermo, PolySize>::clone() const
{
    return autoPtr<logPolynomialTransport<Thermo, PolySize>>
    (
        new logPolynomialTransport<Thermo, PolySize>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::scalar Foam::logPolynomialTransport<Thermo, PolySize>::mu
(
    const scalar p,
    const scalar T
) const
{
    return exp(muLogCoeffs_.value(log(T)));
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::logPolynomialTransport<Thermo, PolySize>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return exp(kappaLogCoeffs_.value(log(T)));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline void Foam::logPolynomialTransport<Thermo, PolySize>::operator+=
(
    const logPolynomialTransport<Thermo, PolySize>& pt
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(pt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = pt.Y()/this->Y();

        muLogCoeffs_ = Y1*muLogCoeffs_ + Y2*pt.muLogCoeffs_;
        kappaLogCoeffs_ = Y1*kappaLogCoeffs_ + Y2*pt.kappaLogCoeffs_;
    }
}


template<class Thermo, int PolySize>
inline void Foam::logPolynomialTransport<Thermo, PolySize>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::logPolynomialTransport<Thermo, PolySize> Foam::operator+
(
    const logPolynomialTransport<Thermo, PolySize>& pt1,
    const logPolynomialTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
    );

    if (mag(t.Y()) < small)
    {
        return logPolynomialTransport<Thermo>
        (
            t,
            0,
            pt1.muLogCoeffs_,
            pt1.kappaLogCoeffs_
        );
    }
    else
    {
        scalar Y1 = pt1.Y()/t.Y();
        scalar Y2 = pt2.Y()/t.Y();

        return logPolynomialTransport<Thermo, PolySize>
        (
            t,
            Y1*pt1.muLogCoeffs_ + Y2*pt2.muLogCoeffs_,
            Y1*pt1.kappaLogCoeffs_ + Y2*pt2.kappaLogCoeffs_
        );
    }
}


template<class Thermo, int PolySize>
inline Foam::logPolynomialTransport<Thermo, PolySize> Foam::operator*
(
    const scalar s,
    const logPolynomialTransport<Thermo, PolySize>& pt
)
{
    return logPolynomialTransport<Thermo, PolySize>
    (
        s*static_cast<const Thermo&>(pt),
        pt.muLogCoeffs_,
        pt.kappaLogCoeffs_
    );
}


// ************************************************************************* //
