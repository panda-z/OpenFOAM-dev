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

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "vtkPVFoamReader.h"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "Cloud.H"
#include "LagrangianMesh.H"
#include "surfaceFields.H"

// Local includes
#include "vtkPVFoamAddToSelection.H"

// VTK includes
#include "vtkDataArraySelection.h"


// * * * * * * * * * * * * * * * Private Classes * * * * * * * * * * * * * * //

namespace Foam
{

//- A class for reading zone information without requiring a mesh
template<class ZonesType>
class zonesEntries
:
    public regIOobject,
    public PtrList<entry>
{

public:

    // Constructors

        explicit zonesEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream(ZonesType::typeName))
        {
            close();
        }

   // Member functions

        bool writeData(Ostream&) const
        {
            NotImplemented;
            return false;
        }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZonesType>
Foam::wordList Foam::vtkPVFoam::getZoneNames(const ZonesType& zones) const
{
    wordList names(zones.size());
    label nZone = 0;

    forAll(zones, zoneI)
    {
        if (zones[zoneI].size())
        {
            names[nZone++] = zones[zoneI].name();
        }
    }
    names.setSize(nZone);

    return names;
}


template<class ZonesType>
Foam::wordList Foam::vtkPVFoam::getZoneNames(const word& zonesName) const
{
    wordList names;

    // Mesh not loaded - read from file
    typeIOobject<ZonesType> ioObj
    (
        zonesName,
        dbPtr_().findInstance
        (
            meshDir_,
            zonesName,
            IOobject::READ_IF_PRESENT
        ),
        meshDir_,
        dbPtr_(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (ioObj.headerOk())
    {
        const zonesEntries<ZonesType> zones(ioObj);

        names.setSize(zones.size());
        forAll(zones, zoneI)
        {
            names[zoneI] = zones[zoneI].keyword();
        }
    }

    return names;
}


void Foam::vtkPVFoam::updateInfoInternalMesh
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // Determine mesh parts (internalMesh, patches...)
    //- Add internal mesh as first entry
    arrayRangeVolume_.reset(arraySelection->GetNumberOfArrays());
    arraySelection->AddArray
    (
        "internalMesh"
    );
    arrayRangeVolume_ += 1;

    if (debug)
    {
        // Just for debug info
        getSelectedArrayEntries(arraySelection);
    }
}


void Foam::vtkPVFoam::updateInfolagrangian
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        InfoInFunction << nl
            << "    " << dbPtr_->timePath()/lagrangian::cloud::prefix << endl;
    }

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(lagrangian::cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/lagrangian::cloud::prefix;
    }

    arrayRangelagrangian_.reset(arraySelection->GetNumberOfArrays());

    // Generate a list of lagrangian clouds across all times
    HashSet<fileName> cloudDirs;

    // Get times list. Flush first to force refresh.
    fileHandler().flush();
    instantList times = dbPtr_().times();
    forAll(times, timei)
    {
        cloudDirs +=
            fileHandler().readDir
            (
                dbPtr_->path()/times[timei].name()/lagrangianPrefix,
                fileType::directory
            );
    }

    forAllConstIter(HashSet<fileName>, cloudDirs, cloudIter)
    {
        arraySelection->AddArray
        (
            (cloudIter.key() + " - lagrangian").c_str()
        );
    }

    arrayRangelagrangian_ += cloudDirs.size();

    if (debug)
    {
        // Just for debug info
        getSelectedArrayEntries(arraySelection);
    }
}


void Foam::vtkPVFoam::updateInfoLagrangian
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        InfoInFunction << nl
            << "    " << dbPtr_->timePath()/LagrangianMesh::prefix << endl;
    }

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    const fileName LagrangianPrefix =
        meshRegion_ == polyMesh::defaultRegion
      ? fileName(LagrangianMesh::prefix)
      : meshRegion_/LagrangianMesh::prefix;

    arrayRangeLagrangian_.reset(arraySelection->GetNumberOfArrays());

    // Generate a list of Lagrangian meshes across all times
    HashSet<fileName> LagrangianDirs;

    // Get times list. Flush first to force refresh.
    fileHandler().flush();
    instantList times = dbPtr_().times();
    forAll(times, timei)
    {
        LagrangianDirs +=
            fileHandler().readDir
            (
                dbPtr_->path()/times[timei].name()/LagrangianPrefix,
                fileType::directory
            );
    }

    forAllConstIter(HashSet<fileName>, LagrangianDirs, LagrangianIter)
    {
        arraySelection->AddArray
        (
            (LagrangianIter.key() + " - Lagrangian").c_str()
        );
    }

    arrayRangeLagrangian_ += LagrangianDirs.size();

    if (debug)
    {
        // Just for debug info
        getSelectedArrayEntries(arraySelection);
    }
}


void Foam::vtkPVFoam::updateInfoPatches
(
    vtkDataArraySelection* arraySelection,
    stringList& enabledEntries,
    const bool first
)
{
    if (debug)
    {
        InfoInFunction
            << " [meshPtr=" << (meshPtr_ ? "set" : "nullptr") << "]" << endl;
    }


    HashSet<string> enabledEntriesSet(enabledEntries);

    arrayRangePatches_.reset(arraySelection->GetNumberOfArrays());

    int nPatches = 0;
    if (meshPtr_)
    {
        const polyBoundaryMesh& patches = meshPtr_->boundaryMesh();
        const HashTable<labelList, word>& groups = patches.groupPatchIndices();
        const wordList allPatchNames = patches.names();

        // Add patch groups

        for
        (
            HashTable<labelList, word>::const_iterator iter = groups.begin();
            iter != groups.end();
            ++iter
        )
        {
            const word& groupName = iter.key();
            const labelList& patchIDs = iter();

            label nFaces = 0;
            forAll(patchIDs, i)
            {
                nFaces += patches[patchIDs[i]].size();
            }

            // Valid patch if nFace > 0 - add patch to GUI list
            if (nFaces)
            {
                string vtkGrpName = groupName + " - group";
                arraySelection->AddArray(vtkGrpName.c_str());

                ++nPatches;

                if (enabledEntriesSet.found(vtkGrpName))
                {
                    if (!reader_->GetShowGroupsOnly())
                    {
                        enabledEntriesSet.erase(vtkGrpName);
                        forAll(patchIDs, i)
                        {
                            const polyPatch& pp = patches[patchIDs[i]];
                            if (pp.size())
                            {
                                string vtkPatchName
                                (
                                    pp.name() + " - " + pp.type()
                                );

                                enabledEntriesSet.insert(vtkPatchName);
                            }
                        }
                    }
                }
            }
        }


        // Add patches

        if (!reader_->GetShowGroupsOnly())
        {
            forAll(patches, patchi)
            {
                const polyPatch& pp = patches[patchi];

                if (pp.size())
                {
                    const string vtkPatchName = pp.name() + " - " + pp.type();

                    // Add patch to GUI list
                    arraySelection->AddArray(vtkPatchName.c_str());

                    ++nPatches;
                }
            }
        }
    }
    else
    {
        // Mesh not loaded - read from file
        // but this could fail if we've supplied a bad region name
        typeIOobject<polyBoundaryMesh> ioObj
        (
            "boundary",
            dbPtr_().findInstance
            (
                meshDir_,
                "boundary",
                IOobject::READ_IF_PRESENT
            ),
            meshDir_,
            dbPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        // This should only ever fail if the mesh region doesn't exist
        if (ioObj.headerOk())
        {
            polyBoundaryMeshEntries patchEntries(ioObj);


            // Read patches and determine sizes

            wordList names(patchEntries.size());
            labelList sizes(patchEntries.size());

            forAll(patchEntries, patchi)
            {
                const dictionary& patchDict = patchEntries[patchi].dict();

                sizes[patchi] = patchDict.lookup<label>("nFaces");
                names[patchi] = patchEntries[patchi].keyword();
            }


            // Add (non-zero) patch groups to the list of mesh parts

            HashTable<labelList, word> groups(patchEntries.size());

            forAll(patchEntries, patchi)
            {
                const dictionary& patchDict = patchEntries[patchi].dict();

                wordList groupNames;
                patchDict.readIfPresent("inGroups", groupNames);

                forAll(groupNames, groupI)
                {
                    HashTable<labelList, word>::iterator iter = groups.find
                    (
                        groupNames[groupI]
                    );
                    if (iter != groups.end())
                    {
                        iter().append(patchi);
                    }
                    else
                    {
                        groups.insert(groupNames[groupI], labelList(1, patchi));
                    }
                }
            }

            for
            (
                HashTable<labelList, word>::const_iterator iter =
                    groups.begin();
                iter != groups.end();
                ++iter
            )
            {
                const word& groupName = iter.key();
                const labelList& patchIDs = iter();

                label nFaces = 0;
                forAll(patchIDs, i)
                {
                    nFaces += sizes[patchIDs[i]];
                }

                // Valid patch if nFace > 0 - add patch to GUI list
                if (nFaces)
                {
                    string vtkGrpName = groupName + " - group";
                    arraySelection->AddArray(vtkGrpName.c_str());

                    ++nPatches;

                    if (enabledEntriesSet.found(vtkGrpName))
                    {
                        if (!reader_->GetShowGroupsOnly())
                        {
                            enabledEntriesSet.erase(vtkGrpName);
                            forAll(patchIDs, i)
                            {
                                if (sizes[patchIDs[i]])
                                {
                                    const word patchType
                                    (
                                        patchEntries[patchIDs[i]].dict().lookup
                                        (
                                            "type"
                                        )
                                    );

                                    string vtkPatchName
                                    (
                                        names[patchIDs[i]] + " - " + patchType
                                    );

                                    enabledEntriesSet.insert(vtkPatchName);
                                }
                            }
                        }
                    }
                }
            }


            // Add (non-zero) patches to the list of mesh parts

            if (!reader_->GetShowGroupsOnly())
            {
                const wordReList defaultPatchTypes
                (
                    configDict_.lookupOrDefault
                    (
                        "defaultPatchTypes",
                        wordReList(wordList({"patch", "wall"}))
                    )
                );

                forAll(names, patchi)
                {
                    // Valid patch if nFace > 0 - add patch to GUI list
                    if (sizes[patchi])
                    {
                        const word patchType
                        (
                            patchEntries[patchi].dict().lookup("type")
                        );

                        const string vtkPatchName
                        (
                            names[patchi] + " - " + patchType
                        );

                        arraySelection->AddArray(vtkPatchName.c_str());

                        if (first)
                        {
                            if (findStrings(defaultPatchTypes, patchType))
                            {
                                enabledEntriesSet.insert(vtkPatchName);
                            }
                        }

                        ++nPatches;
                    }
                }
            }
        }
    }
    arrayRangePatches_ += nPatches;

    // Update enabled entries in case of group selection
    enabledEntries = enabledEntriesSet.toc();

    if (debug)
    {
        // Just for debug info
        getSelectedArrayEntries(arraySelection);
    }
}


void Foam::vtkPVFoam::updateInfoZones
(
    vtkDataArraySelection* arraySelection
)
{
    if (!reader_->GetIncludeZones())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction
            << " [meshPtr=" << (meshPtr_ ? "set" : "nullptr") << "]" << endl;
    }

    wordList namesLst;

    // cellZones information
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->cellZones());
    }
    else
    {
        namesLst = getZoneNames<cellZoneList>("cellZones");
    }

    arrayRangeCellZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            (namesLst[elemI] + " - cellZone").c_str()
        );
    }
    arrayRangeCellZones_ += namesLst.size();


    // faceZones information
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->faceZones());
    }
    else
    {
        namesLst = getZoneNames<faceZoneList>("faceZones");
    }

    arrayRangeFaceZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            (namesLst[elemI] + " - faceZone").c_str()
        );
    }
    arrayRangeFaceZones_ += namesLst.size();


    // pointZones information
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->pointZones());
    }
    else
    {
        namesLst = getZoneNames<pointZoneList>("pointZones");
    }

    arrayRangePointZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            (namesLst[elemI] + " - pointZone").c_str()
        );
    }
    arrayRangePointZones_ += namesLst.size();

    if (debug)
    {
        // Just for debug info
        getSelectedArrayEntries(arraySelection);
    }
}


void Foam::vtkPVFoam::updateInfoSets
(
    vtkDataArraySelection* arraySelection
)
{
    if (!reader_->GetIncludeSets())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction << endl;
    }

    // Add names of sets. Search for last time directory with a sets
    // subdirectory. Take care not to search beyond the last mesh.

    word facesInstance = dbPtr_().findInstance
    (
        meshDir_,
        "faces",
        IOobject::READ_IF_PRESENT
    );

    word setsInstance = dbPtr_().findInstance
    (
        meshDir_/"sets",
        word::null,
        IOobject::READ_IF_PRESENT,
        facesInstance
    );

    IOobjectList objects(dbPtr_(), setsInstance, meshDir_/"sets");

    if (debug)
    {
        Info<< "     Foam::vtkPVFoam::updateInfoSets read "
            << objects.names() << " from " << setsInstance << endl;
    }


    arrayRangeCellSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangeCellSets_ += addToSelection<cellSet>
    (
        arraySelection,
        objects,
        " - cellSet"
    );

    arrayRangeFaceSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangeFaceSets_ += addToSelection<faceSet>
    (
        arraySelection,
        objects,
        " - faceSet"
    );

    arrayRangePointSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangePointSets_ += addToSelection<pointSet>
    (
        arraySelection,
        objects,
        " - pointSet"
    );

    if (debug)
    {
        // Just for debug info
        getSelectedArrayEntries(arraySelection);
    }
}


void Foam::vtkPVFoam::updateInfoFields()
{
    if (debug)
    {
        InfoInFunction
            << " [meshPtr=" << (meshPtr_ ? "set" : "nullptr") << "]"
            << endl;
    }

    vtkDataArraySelection* fieldSelection = reader_->GetFieldSelection();

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    word regionPrefix;
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        regionPrefix = meshRegion_;
    }

    const Time& runTime = dbPtr_();
    const instantList times = runTime.times();

    stringList enabledEntries;

    if (fieldSelection->GetNumberOfArrays() == 0 && !meshPtr_)
    {
        const wordReList defaultFieldRes
        (
            configDict_.lookupOrDefault
            (
                "defaultFields",
                wordReList(wordList{"p", "U"})
            )
        );

        wordHashSet objectNameSet;
        forAll(times, timei)
        {
            // Search for list of objects for this time and mesh region
            IOobjectList objects(runTime, times[timei].name(), regionPrefix);

            forAllConstIter(IOobjectList, objects, iter)
            {
                objectNameSet.insert(iter.key());
            }
        }

        const wordList objectNames(objectNameSet.toc());
        const labelList defaultFields
        (
            findStrings(defaultFieldRes, objectNames)
        );

        enabledEntries.setSize(defaultFields.size());
        forAll(defaultFields, i)
        {
            enabledEntries[i] = objectNames[defaultFields[i]];
        }
    }
    else
    {
        // Preserve the enabled selections
        enabledEntries = getSelectedArrayEntries(fieldSelection);
    }

    fieldSelection->RemoveAllArrays();

    forAll(times, timei)
    {
        // Search for list of objects for this time and mesh region
        IOobjectList objects(runTime, times[timei].name(), regionPrefix);

        addFieldsToSelection<volMesh>(fieldSelection, objects);
        addInternalFieldsToSelection<volMesh>(fieldSelection, objects);
        addFieldsToSelection<surfaceMesh>(fieldSelection, objects);
        addFieldsToSelection<pointMesh>(fieldSelection, objects);
    }

    // Restore the enabled selections
    setSelectedArrayEntries(fieldSelection, enabledEntries);
}


void Foam::vtkPVFoam::updateInfolagrangianFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    vtkDataArraySelection* fieldSelection =
        reader_->GetlagrangianFieldSelection();

    // Preserve the enabled selections
    stringList enabledEntries = getSelectedArrayEntries(fieldSelection);
    fieldSelection->RemoveAllArrays();

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(lagrangian::cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/lagrangian::cloud::prefix;
    }

    // Add the available fields from all clouds and all time directories.
    // Differing sets of fields from multiple clouds get combined into a single
    // set. ParaView will display "(partial)" after field names that only apply
    // to some of the clouds.
    const arrayRange& range = arrayRangelagrangian_;

    fileHandler().flush();
    for (label partId = range.start(); partId < range.end(); ++ partId)
    {
        const instantList times = dbPtr_().times();
        forAll(times, timei)
        {
            IOobjectList objects
            (
                dbPtr_(),
                times[timei].name(),
                lagrangianPrefix/getPartName(partId)
            );

            addToSelection<IOField<label>>(fieldSelection, objects);
            addToSelection<IOField<scalar>>(fieldSelection, objects);
            addToSelection<IOField<vector>>(fieldSelection, objects);
            addToSelection<IOField<sphericalTensor>>(fieldSelection, objects);
            addToSelection<IOField<symmTensor>>(fieldSelection, objects);
            addToSelection<IOField<tensor>>(fieldSelection, objects);
        }
    }

    // Restore the enabled selections
    setSelectedArrayEntries(fieldSelection, enabledEntries);
}


void Foam::vtkPVFoam::updateInfoLagrangianFields()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    vtkDataArraySelection* fieldSelection =
        reader_->GetLagrangianFieldSelection();

    // Preserve the enabled selections
    stringList enabledEntries = getSelectedArrayEntries(fieldSelection);
    fieldSelection->RemoveAllArrays();

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    const fileName LagrangianPrefix =
        meshRegion_ == polyMesh::defaultRegion
      ? fileName(LagrangianMesh::prefix)
      : meshRegion_/LagrangianMesh::prefix;

    // Add the available fields from all clouds and all time directories.
    // Differing sets of fields from multiple clouds get combined into a single
    // set. ParaView will display "(partial)" after field names that only apply
    // to some of the clouds.
    const arrayRange& range = arrayRangeLagrangian_;

    fileHandler().flush();

    const instantList times = dbPtr_().times();

    for (label partId = range.start(); partId < range.end(); ++ partId)
    {
        forAll(times, timei)
        {
            IOobjectList objects
            (
                dbPtr_(),
                times[timei].name(),
                LagrangianPrefix/getPartName(partId)
            );

            #define ADD_TO_SELECTION(Type, GeoField) \
                addToSelection<GeoField<Type>>(fieldSelection, objects);
            ADD_TO_SELECTION(label, LagrangianField);
            FOR_ALL_FIELD_TYPES(ADD_TO_SELECTION, LagrangianField);
            ADD_TO_SELECTION(label, LagrangianInternalField);
            FOR_ALL_FIELD_TYPES(ADD_TO_SELECTION, LagrangianInternalField);
            #undef ADD_TO_SELECTION
        }
    }

    // Restore the enabled selections
    setSelectedArrayEntries(fieldSelection, enabledEntries);
}


// ************************************************************************* //
