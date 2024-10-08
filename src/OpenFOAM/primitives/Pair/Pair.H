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

Class
    Foam::Pair

Description
    An ordered pair of two objects of type \<T\> with first() and second()
    elements.

See also
    Foam::Tuple2 for storing two objects of dissimilar types.

\*---------------------------------------------------------------------------*/

#ifndef Pair_H
#define Pair_H

#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class Pair;

template<class Type>
void writeEntry(Ostream& os, const Pair<Type>&);

template<class Type>
inline Istream& operator>>(Istream&, Pair<Type>&);

template<class Type>
inline Ostream& operator<<(Ostream&, const Pair<Type>&);


/*---------------------------------------------------------------------------*\
                            Class Pair Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Pair
:
    public FixedList<Type, 2>
{

public:

    // Constructors

        //- Null constructor
        inline Pair()
        {}

        //- Construct from components
        inline Pair(const Type& f, const Type& s)
        {
            first() = f;
            second() = s;
        }

        //- Construct from FixedList
        inline Pair(const FixedList<Type, 2>& fl)
        :
            FixedList<Type, 2>(fl)
        {}

        //- Construct from Istream
        inline Pair(Istream& is)
        {
            is >> *this;
        }


    // Member Functions

        //- Return first
        inline const Type& first() const
        {
            return this->operator[](0);
        }

        //- Return first
        inline Type& first()
        {
            return this->operator[](0);
        }

        //- Return second
        inline const Type& second() const
        {
            return this->operator[](1);
        }

        //- Return second
        inline Type& second()
        {
            return this->operator[](1);
        }

        //- Return other
        inline const Type& other(const Type& a) const
        {
            if (first() == second())
            {
                FatalErrorInFunction
                    << "Call to other only valid for Pair with differing"
                    << " elements:" << *this << abort(FatalError);
                return first();
            }
            else if (first() == a)
            {
                return second();
            }
            else
            {
                if (second() != a)
                {
                    FatalErrorInFunction
                        << "Pair " << *this
                        << " does not contain " << a << abort(FatalError);
                }
                return first();
            }
        }


        //- Compare Pairs
        //  Returning:
        //  -  0: different
        //  - +1: identical
        //  - -1: same pair, but reversed order
        static inline int compare(const Pair<Type>& a, const Pair<Type>& b)
        {
            if (a == b)
            {
                return 1;
            }
            else if (a == reverse(b))
            {
                return -1;
            }
            else
            {
                return 0;
            }
        }
};


template<class Type>
Pair<Type> reverse(const Pair<Type>& p)
{
    return Pair<Type>(p.second(), p.first());
}


template<class Type>
bool operator==(const Pair<Type>& a, const Pair<Type>& b)
{
    return (a.first() == b.first() && a.second() == b.second());
}


template<class Type>
bool operator!=(const Pair<Type>& a, const Pair<Type>& b)
{
    return !(a == b);
}


template<class Type>
bool operator<(const Pair<Type>& a, const Pair<Type>& b)
{
    return
    (
        a.first() < b.first()
     ||
        (
            !(b.first() < a.first())
         && a.second() < b.second()
        )
    );
}


template<class Type>
bool operator<=(const Pair<Type>& a, const Pair<Type>& b)
{
    return !(b < a);
}


template<class Type>
bool operator>(const Pair<Type>& a, const Pair<Type>& b)
{
    return (b < a);
}


template<class Type>
bool operator>=(const Pair<Type>& a, const Pair<Type>& b)
{
    return !(a < b);
}


template<class Type>
inline Istream& operator>>(Istream& is, Pair<Type>& p)
{
    is.readBegin("Pair");
    is >> p.first() >> p.second();
    is.readEnd("Pair");

    // Check state of Istream
    is.check("operator>>(Istream&, Pair<Type>&)");

    return is;
}


template<class Type>
inline Ostream& operator<<(Ostream& os, const Pair<Type>& p)
{
    os  << token::BEGIN_LIST
        << p.first() << token::SPACE << p.second()
        << token::END_LIST;

    return os;
}


template<class Type>
void writeEntry(Ostream& os, const Pair<Type>& p)
{
    os << p;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
