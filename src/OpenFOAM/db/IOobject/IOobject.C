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

#include "IOobject.H"
#include "Time.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IOobject, 0);
}

const Foam::NamedEnum<Foam::IOobject::fileCheckTypes, 4>
Foam::IOobject::fileCheckTypesNames
{
    "timeStamp",
    "timeStampMaster",
    "inotify",
    "inotifyMaster"
};

// Default fileCheck type
Foam::IOobject::fileCheckTypes Foam::IOobject::fileModificationChecking
(
    Foam::debug::namedEnumOptimisationSwitch
    (
        "fileModificationChecking",
        fileCheckTypesNames,
        fileModificationChecking
    )
);

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

bool Foam::IOobject::fileNameComponents
(
    const fileName& path,
    fileName& instance,
    fileName& local,
    word& name
)
{
    instance.clear();
    local.clear();
    name.clear();

    // called with directory
    if (isDir(path))
    {
        WarningInFunction
            << " called with directory: " << path << endl;

        return false;
    }

    if (path.isAbsolute())
    {
        string::size_type last = path.rfind('/');
        instance = path.substr(0, last);

        // Check afterwards
        name.string::operator=(path.substr(last+1));
    }
    else
    {
        string::size_type first = path.find('/');

        if (first == string::npos)
        {
            // no '/' found - no instance or local

            // check afterwards
            name.string::operator=(path);
        }
        else
        {
            instance = path.substr(0, first);

            string::size_type last = path.rfind('/');
            if (last > first)
            {
                // with local
                local = path.substr(first+1, last-first-1);
            }

            // check afterwards
            name.string::operator=(path.substr(last+1));
        }
    }


    // Check for valid (and stripped) name, regardless of the debug level
    if (name.empty() || string::stripInvalid<word>(name))
    {
        WarningInFunction
            << "has invalid word for name: \"" << name
            << "\"\nwhile processing path: " << path << endl;

        return false;
    }

    return true;
}


Foam::word Foam::IOobject::group(const word& name)
{
    word::size_type i = name.find_last_of('.');

    if (i == word::npos || i == 0)
    {
        return word::null;
    }
    else
    {
        return name.substr(i+1, word::npos);
    }
}


Foam::word Foam::IOobject::member(const word& name)
{
    word::size_type i = name.find_last_of('.');

    if (i == word::npos || i == 0)
    {
        return name;
    }
    else
    {
        return name.substr(0, i);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject
)
:
    name_(name),
    headerClassName_(typeName),
    note_(),
    instance_(instance),
    local_(),
    db_(registry),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    objState_(GOOD)
{
    if (objectRegistry::debug)
    {
        InfoInFunction
            << "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const fileName& local,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject
)
:
    name_(name),
    headerClassName_(typeName),
    note_(),
    instance_(instance),
    local_(local),
    db_(registry),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    objState_(GOOD)
{
    if (objectRegistry::debug)
    {
        InfoInFunction
            << "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const fileName& path,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject
)
:
    name_(),
    headerClassName_(typeName),
    note_(),
    instance_(),
    local_(),
    db_(registry),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    objState_(GOOD)
{
    if (!fileNameComponents(path, instance_, local_, name_))
    {
        FatalErrorInFunction
            << " invalid path specification"
            << exit(FatalError);
    }

    if (objectRegistry::debug)
    {
        InfoInFunction
            << "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const IOobject& io,
    const objectRegistry& registry
)
:
    name_(io.name_),
    headerClassName_(io.headerClassName_),
    note_(io.note_),
    instance_(io.instance_),
    local_(io.local_),
    db_(registry),
    rOpt_(io.rOpt_),
    wOpt_(io.wOpt_),
    registerObject_(io.registerObject_),
    objState_(io.objState_)
{}


Foam::IOobject::IOobject
(
    const IOobject& io,
    const word& name
)
:
    name_(name),
    headerClassName_(io.headerClassName_),
    note_(io.note_),
    instance_(io.instance_),
    local_(io.local_),
    db_(io.db_),
    rOpt_(io.rOpt_),
    wOpt_(io.wOpt_),
    registerObject_(io.registerObject_),
    objState_(io.objState_)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::IOobject::~IOobject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::IOobject::db() const
{
    return db_;
}


const Foam::Time& Foam::IOobject::time() const
{
    return db_.time();
}


Foam::word Foam::IOobject::group() const
{
    return group(name_);
}


Foam::word Foam::IOobject::member() const
{
    return member(name_);
}


const Foam::fileName& Foam::IOobject::rootPath() const
{
    return time().rootPath();
}


const Foam::fileName& Foam::IOobject::caseName(const bool global) const
{
    if (global)
    {
        return time().globalCaseName();
    }
    else
    {
        return time().caseName();
    }
}


Foam::fileName& Foam::IOobject::instance() const
{
    return instance_;
}


void Foam::IOobject::updateInstance() const
{
    if
    (
        !instance_.isAbsolute()
     && instance_ != time().system()
     && instance_ != time().constant()
     && instance_ != time().name()
    )
    {
        scalar timeValue;
        if (readScalar(instance_.c_str(), timeValue))
        {
            instance_ = time().name();
        }
    }
}


void Foam::IOobject::updateTimeInstance() const
{
    instance_ = time().name();
}


Foam::fileName Foam::IOobject::path(const bool global) const
{
    if (instance_.isAbsolute())
    {
        return instance_;
    }
    else
    {
        return rootPath()/caseName(global)/instance()/db_.dbDir()/local();
    }
}


Foam::fileName Foam::IOobject::relativePath() const
{
    if (instance().isAbsolute())
    {
        return instance();
    }
    else
    {
        return instance()/db_.dbDir()/local();
    }
}


Foam::fileName Foam::IOobject::filePath(const bool global) const
{
    return fileHandler().filePath(global, *this);
}


void Foam::IOobject::operator=(const IOobject& io)
{
    name_ = io.name_;
    headerClassName_ = io.headerClassName_;
    note_ = io.note_;
    instance_ = io.instance_;
    local_ = io.local_;
    rOpt_ = io.rOpt_;
    wOpt_ = io.wOpt_;
    objState_ = io.objState_;
}


// ************************************************************************* //
