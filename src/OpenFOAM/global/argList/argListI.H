/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "argList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::string& Foam::argList::commandLine() const
{
    return argListStr_;
}


inline const Foam::word& Foam::argList::executable() const
{
    return executable_;
}


inline const Foam::fileName& Foam::argList::rootPath() const
{
    return rootPath_;
}


inline const Foam::fileName& Foam::argList::caseName() const
{
    return case_;
}


inline const Foam::fileName& Foam::argList::globalCaseName() const
{
    return globalCase_;
}


inline const Foam::ParRunControl& Foam::argList::parRunControl() const
{
    return parRunControl_;
}


inline Foam::fileName Foam::argList::path() const
{
    return rootPath()/caseName();
}


inline const Foam::stringList& Foam::argList::args() const
{
    return args_;
}


inline Foam::stringList& Foam::argList::args()
{
    return args_;
}


inline const Foam::string& Foam::argList::arg(const label index) const
{
    return args_[index];
}


inline Foam::label Foam::argList::size() const
{
    return args_.size();
}


inline const Foam::HashTable<Foam::string>& Foam::argList::options() const
{
    return options_;
}


inline Foam::HashTable<Foam::string>& Foam::argList::options()
{
    return options_;
}


inline const Foam::string& Foam::argList::option(const word& opt) const
{
    return options_[opt];
}


inline bool Foam::argList::optionFound(const word& opt) const
{
    return options_.found(opt);
}


inline Foam::IStringStream Foam::argList::optionLookup(const word& opt) const
{
    return IStringStream(options_[opt]);
}


// * * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * //

namespace Foam
{
    // Template specialisation for string
    template<>
    inline Foam::string
    Foam::argList::argRead<Foam::string>(const label index) const
    {
        return args_[index];
    }

    // Template specialisation for word
    template<>
    inline Foam::word
    Foam::argList::argRead<Foam::word>(const label index) const
    {
        return args_[index];
    }

    // Template specialisation for fileName
    template<>
    inline Foam::fileName
    Foam::argList::argRead<Foam::fileName>(const label index) const
    {
        return args_[index];
    }

    // Template specialisation for string
    template<>
    inline Foam::string
    Foam::argList::optionRead<Foam::string>(const word& opt) const
    {
        return options_[opt];
    }

    // Template specialisation for word
    template<>
    inline Foam::word
    Foam::argList::optionRead<Foam::word>(const word& opt) const
    {
        return options_[opt];
    }

    // Template specialisation for fileName
    template<>
    inline Foam::fileName
    Foam::argList::optionRead<Foam::fileName>(const word& opt) const
    {
        return options_[opt];
    }

    inline Foam::label
    Foam::argList::bracketLevel(const Foam::string& str) const
    {
        label level = 0;
        label maxLevel = 0;

        for
        (
            string::const_iterator iter = str.begin();
            iter != str.end();
            ++iter
        )
        {
            const char c = *iter;

            if (c == '(')
            {
                level++;
            }
            else if (c == ')')
            {
                level--;
            }

            maxLevel = max(maxLevel, level);
        }

        return maxLevel;
    }

    template<class T>
    inline Foam::List<T>
    Foam::argList::optionReadList(const word& opt) const
    {
        string argstr(optionLookup(opt).str());

        OStringStream os;
        const T x;
        os << x;

        label level = bracketLevel(os.str()) + 1;

        while (bracketLevel(argstr) < level)
        {
            argstr.insert(0, "(");
            argstr.insert(argstr.size(), ")");
        }

        return readList<T>(IStringStream(argstr)());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline T Foam::argList::argRead(const label index) const
{
    T val;

    IStringStream(args_[index])() >> val;
    return val;
}


template<class T>
inline T Foam::argList::optionRead(const word& opt) const
{
    T val;

    optionLookup(opt)() >> val;
    return val;
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& opt,
    T& val
) const
{
    if (optionFound(opt))
    {
        val = optionRead<T>(opt);
        return true;
    }
    else
    {
        return false;
    }
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& opt,
    T& val,
    const T& deflt
) const
{
    if (optionReadIfPresent<T>(opt, val))
    {
        return true;
    }
    else
    {
        val = deflt;
        return false;
    }
}


template<class T>
inline T Foam::argList::optionLookupOrDefault
(
    const word& opt,
    const T& deflt
) const
{
    if (optionFound(opt))
    {
        return optionRead<T>(opt);
    }
    else
    {
        return deflt;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::string& Foam::argList::operator[](const label index) const
{
    return args_[index];
}


inline const Foam::string& Foam::argList::operator[](const word& opt) const
{
    return options_[opt];
}


// ************************************************************************* //
