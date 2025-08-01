/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

Application
    snappyHexMeshConfig

Description
    Preconfigures blockMeshDict, surfaceFeaturesDict and snappyHexMeshDict
    files based on the case surface geometry files.

    Starting from a standard OpenFOAM case, this utility locates surface
    geometry files, e.g. OBJ, STL format, in the constant/geometry directory.
    It writes out the configuration files for mesh generation with
    snappyHexMesh based on assumptions which can be overridden by options on
    the command line.

    The utility processes the surface geometry files, attempting to anticipate
    their intended purpose, trying in particular to recognise whether the
    domain represents an external or internal flow. If there is a surface
    which is closed, and is either single or surrounds all other surfaces,
    then it is assumed that it forms the external boundary of an internal
    flow. This assumption is overridden if the bounds of the background mesh
    are specified using the '-bounds' option and they are more than 50% larger
    than the surface bounds.

    Surfaces which form boundaries of the domain may contain named regions
    that are intended to become patches in the final mesh. Any surface region
    whose name begins with 'inlet' or 'outlet' will become a patch of the same
    name in the final mesh. On an external surface (for an internal flow),
    regions can be identified as inlets and outlets using the '-inletRegions'
    and '-outletRegions' options, respectively. When either option specifies a
    single region, the resulting patch name will be specifically 'inlet' or
    'outlet', respectively. Surfaces which are contained within the domain,
    which do not surround or intersect other surfaces, are assumed by default
    to be wall patches. Any closed surface which surrounds another (but not an
    external surface) is used to form a cellZone within the mesh. Any surface
    can be specifically identified as a cellZone with the '-cellZones' option,
    with the additional '-baffles' and '-rotatingZones' options available to
    assign a surface to a more specific use.

    The background mesh for snappyHexMesh is a single block generated by
    blockMesh, configured using a blockMeshDict file. The block bounds are
    automatically calculated, but can be overridden by the '-bounds'
    option. The number of cells is calculated to produce a fairly small
    prototype mesh. The cell density can be overridden by the '-nCells' option
    or can be scaled up by an integer factor using the '-refineBackground'
    option. When the background mesh is required to form patches in the final
    mesh, e.g. for an external flow, the user can specify the names and types
    of the patches corresponding to the six block faces using options such as
    '-xMinPatch', '-xMaxPatch', etc. The name and type of the default patch,
    formed from block faces which are not configured, can also be specified
    with the '-defaultPatch' option. The utility provides placeholder entries
    for all block faces unless the '-clearBoundary' option is used. A special
    '-cylindricalBackground' option generates a cylindrical background mesh,
    oriented along the z-axis along x = y = 0.

    The snappyHexMesh configuration is generated automatically, applying a set
    of defaults to the main configuration parameters. By default, implicit
    feature capturing is configured.  Explicit feature capturing can
    alternatively be selected with the '-explicitFeatures' option, when an
    additional surfaceFeaturesDict file is written for the user to generate the
    features files with the surfaceFeatures utility.  Refinement levels can be
    controlled with a range of options including: '-refinementLevel' for the
    baseline refinement level; '-refinementSurfaces' for levels on specific
    surfaces; '-refinementRegions' for levels inside specific surfaces;
    '-refinementBoxes' for quick, box-shaped refinement regions specified by min
    and max bounds; '-refinementDists' for distance-based refinement; and
    '-nCellsBetweenLevels' to control the transition between refinement
    levels. A '-layers' option controls additional layers of cells at specified
    surfaces. The insidePoint parameter is set to '(0 0 0)' by default but can
    be overridden using the '-insidePoint' option.  There is an alternative
    '-insidePoints' option to specify multiple insidePoints to mesh multiple
    disconnected mesh regions.

Usage
    \b snappyHexMeshConfig [OPTIONS]

    Options:

      - \par -baffles \<list\>
        Surfaces that form baffles, e.g. 'helical'

      - \par -bounds \<box\>
        Bounding box of the mesh, e.g. '(-10 -5 0) (10 5 10)'

      - \par -cellZones \<list\>
        Surfaces that form cellZones, e.g. 'porousZone heatSource'

      - \par -clearBoundary,
        Do not set default patch entries, i.e. xMin, xMax, yMin, etc...

      - \par -closedDomain
        Domain does not contain inlets or outlets

      - \par -cylindricalBackground
        Generate a cylindrical background mesh aligned with the z-axis

      - \par -defaultPatch \<entry\>
        Name and type of default patch, '(\<name\> \<type\>)'

      - \par -explicitFeatures,
        Use explicit feature capturing, default is implicit

      - \par -firstLayerThickness \<value\>
        Specify the thickness of the near wall cells for layer addition

      - \par -inletRegions \<list\>
        Inlet regions on an external surface, e.g. 'inletA inletB'

      - \par -insidePoint \<point\>
        Point location inside the region of geometry to be meshed

      - \par -insidePoints \<list\>
        Point locations inside geometry to be meshed, e.g. '(0 0 0) (0 1 0)'

      - \par -layerExpansionRatio \<value\>
        Specify the expansion ratio between layers, default 1.2

      - \par -layers \<entry\>
        Number of layers on specified surfaces, e.g. '(car 3) (ground 4)'

      - \par -minDimCells \<cells\>
        Number of cells in the shortest direction, e.g. 10

      - \par -nCells \<cells\>
        Number of cells in each direction, e.g. '(10 20 30)'

      - \par -nCellsBetweenLevels \<int\>
        Number of cells at successive refinement levels, default 3

      - \par -noBackground
        Do not write a blockMeshDict file

      - \par -outletRegions \<list\>
        Outlet regions on an external surface, e.g. 'outletA outletB'

      - \par -refineBackground \<int\>
        Integer multiplier for the number of cells (>= 1)

      - \par -refinementBoxes \<entry\>
        Refinement boxes specified by '(\<min\> \<max\> \<level\>) (...) '

      - \par -refinementDists \<entry\>
        Refinement distance specified by
        '( (\<surface\> \<dist\> \<level\>) (...) )'

      - \par -refinementLevel \<int\>
        Refinement level used by snappyHexMesh, default 2

      - \par -refinementRegions \<entry\>
        Refinement regions specified by '(\<surface\> \<level\>) (...)'

      - \par -region \<name\>
        Specify alternative mesh region

      - \par -rotatingZones \<list\>
        Surfaces that form rotatingZones, e.g. 'rotatingZone'

      - \par -surface \<file\>
        Single surface geometry file for meshing

      - \par -surfaceLevels \<entry\>
        Refinement level at specified surfaces, e.g. '(pipe 2) (baffles 1)'

      - \par -xMinPatch (-xMaxPatch, -yMinPatch, etc...) \<entry\>
        Name and type of the xMin (xMax, yMin, etc...) patch,
        '(\<name\> \<type\>)'

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "meshingSurface.H"
#include "blockMeshCartesianConfiguration.H"
#include "blockMeshCylindricalConfiguration.H"
#include "snappyHexMeshConfiguration.H"
#include "meshQualityConfiguration.H"
#include "surfaceFeaturesConfiguration.H"
#include "boundBox.H"
#include "searchableSurface.H"
#include "Tuple3.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void readPatchOption
(
    const argList& args,
    HashTable<Pair<word>>& opts,
    const word& name
)
{
    if (args.optionFound(name))
    {
        List<word> patchOption(args.optionReadList<word>(name));

        if (patchOption.size() == 2)
        {
            opts.insert(name, Pair<word>(patchOption[0], patchOption[1]));
        }
        else
        {
            FatalErrorInFunction
                << "Argument to option '-" << name
                << "' is " << args.option(name) << nl
                << "It should be of the form \"<name> <type>\""
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::usageMin = 32;
    argList::usageMax = 105;

    argList::addNote
    (
        "Writes blockMeshDict, snappyHexMeshDict, surfaceFeaturesDict "
        "and meshQualityDict files.\n"
        "Requires surface geometry files as inputs.\n"
        "For more information, see 'Description' in snappyHexMeshConfig.C "
        "or run\n\n    foamInfo snappyHexMeshConfig"
    );

    #include "removeCaseOptions.H"
    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "rm",
        "delete mesh configuration files"
    );

    argList::addOption
    (
        "surface",
        "file",
        "single surface geometry file for meshing"
    );

    argList::addOption
    (
        "nCells",
        "cells",
        "number of cells in each direction, e.g. '(10 20 30)'"
    );

    argList::addOption
    (
        "minDimCells",
        "int",
        "number of cells in the shortest direction, e.g. 10"
    );

    argList::addOption
    (
        "bounds",
        "box",
        "bounding box of the mesh, e.g. '(-10 -5 0) (10 5 10)'"
    );

    argList::addBoolOption
    (
        "cylindricalBackground",
        "generate a cylindrical background mesh aligned with the z-axis"
    );

    argList::addBoolOption
    (
        "noBackground",
        "do not write a blockMeshDict file"
    );

    argList::addOption
    (
        "refineBackground",
        "int",
        "integer multiplier for the number of cells (>= 1)"
    );

    argList::addOption
    (
        "refinementLevel",
        "int",
        "refinement level used by snappyHexMesh, default 2"
    );

    argList::addOption
    (
        "surfaceLevels",
        "entry",
        "refinement level at specified surfaces, e.g. '(pipe 2) (baffles 1)'"
    );

    argList::addOption
    (
        "refinementRegions",
        "entry",
        "refinement regions specified by '(<surface> <level>) (...)'"
    );

    argList::addOption
    (
        "refinementBoxes",
        "entry",
        "refinement boxes specified by '(<min> <max> <level>) (...)'"
    );

    argList::addOption
    (
        "refinementDists",
        "entry",
        "refinement distance specified by "
        "'(<surface> <dist> <level>) (...)'"
    );

    argList::addOption
    (
        "defaultPatch",
        "entry",
        "name and type of default patch, '<name> <type>'"
    );

    List<word> patches(blockMeshCartesianConfiguration::patches);

    forAll(patches, i)
    {
        argList::addOption
        (
            patches[i] + "Patch",
            "entry",
            "patch in the "
          + patches[i]
          + " direction, format '<name> <type>'"
        );
    }

    argList::addBoolOption
    (
        "clearBoundary",
        "do not set default patch entries, i.e. xMin, xMax, etc"
    );

    argList::addBoolOption
    (
        "explicitFeatures",
        "use explicit feature capturing"
    );

    argList::addOption
    (
        "layers",
        "entry",
        "number of layers on specified surfaces, e.g. '(car 3) (ground 4)'"
    );

    argList::addOption
    (
        "firstLayerThickness",
        "value",
        "specify the thickness of the near wall cells for layer addition"
    );

    argList::addOption
    (
        "layerExpansionRatio",
        "value",
        "specify the expansion ratio between layers, default 1.2"
    );

    argList::addOption
    (
        "cellZones",
        "list",
        "surfaces that form cellZones, e.g. 'porousZone heatSource'"
    );

    argList::addOption
    (
        "rotatingZones",
        "list",
        "surfaces that form rotatingZones, e.g. 'rotatingZone'"
    );

    argList::addOption
    (
        "baffles",
        "list",
        "surfaces that form baffles, e.g. 'helical'"
    );

    argList::addOption
    (
        "insidePoint",
        "point",
        "point location inside the region of geometry to be meshed"
    );

    argList::addOption
    (
        "insidePoints",
        "list",
        "point locations inside the geometry, e.g. '(0 1 0) (1 1 1)'"
    );

    argList::addOption
    (
        "nCellsBetweenLevels",
        "int",
        "number of cells at successive refinement levels, default 3"
    );

    argList::addOption
    (
        "inletRegions",
        "list",
        "inlet regions on an external surface, e.g. 'inletA inletB'"
    );

    argList::addOption
    (
        "outletRegions",
        "list",
        "outlet regions on an external surface, e.g. 'outletA outletB'"
    );

    argList::addBoolOption
    (
        "closedDomain",
        "domain does not contain inlets or outlets"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    word regionName;
    word dir(runTime.system());

    if (args.optionReadIfPresent("region", regionName))
    {
        dir = runTime.system()/regionName;

        if (!isDir(dir))
        {
            mkDir(dir);
        }
    }

    if (args.optionFound("rm"))
    {
        wordList dicts
        {
            "snappyHexMeshDict",
            "blockMeshDict",
            "meshQualityDict",
            "surfaceFeaturesDict"
        };

        Info<< "Deleting mesh configuration files in '"
            << dir << "'" << endl;
        label count = 0;

        forAll(dicts, i)
        {
            if (rm(dir/dicts[i]))
            {
                Info<< "+ " << dicts[i] << endl;
                ++count;
            }
        }

        if (count == 0)
        {
            Info<< "+ No files to delete" << endl;
        }

        Info<< "\nEnd\n" << endl;

        return 0;
    }

    Info<< "Writing mesh configuration files to '"
        << dir << "'" << nl << endl;

    fileNameList surfaceNames;

    if (args.optionFound("surface"))
    {
        surfaceNames.append(args.optionRead<fileName>("surface"));
    }
    else
    {
        const fileName surfDir
        (
            runTime.constant()/searchableSurface::geometryDir(runTime)
        );

        // Reads files, removing "gz" extensions
        fileNameList files(readDir(surfDir));

        // Check valid extensions and add the path to the names
        forAll(files, i)
        {
            if (!meshingSurface::isSurfaceExt(files[i]))
            {
                continue;
            }

            surfaceNames.append(surfDir/files[i]);
        }

        // Need to exit if no surface geometry files found
        if (surfaceNames.empty())
        {
            FatalErrorInFunction
                << "No surface geometry files found in "
                << surfDir << nl
                << "or provided using the '-surface' option"
                << exit(FatalError);
        }
    }

    wordList cellZoneNames;
    if (args.optionFound("cellZones"))
    {
        cellZoneNames.append(args.optionReadList<word>("cellZones"));
    }

    wordList rotatingZoneNames;
    if (args.optionFound("rotatingZones"))
    {
        rotatingZoneNames.append(args.optionReadList<word>("rotatingZones"));
    }

    wordList baffleNames;
    if (args.optionFound("baffles"))
    {
        baffleNames.append(args.optionReadList<word>("baffles"));
    }

    boundBox bb;
    if (args.optionFound("bounds"))
    {
        List<vector> bounds(args.optionReadList<vector>("bounds"));

        if (bounds.size() != 2)
        {
            FatalErrorInFunction
                << "Argument to '-bounds'"
                << " should be of the form '(<min> <max>)'" << nl
                << "with the <min> and <max> bounds of a bounding box"
                << "\n\nFound instead the argument: "
                << bounds
                << exit(FatalError);
        }

        bb = boundBox(bounds[0], bounds[1]);
        Info<< "Bounding box specified by '-bounds' option: "
            << bb << endl;
    }

    wordList inletRegions;
    if (args.optionFound("inletRegions"))
    {
        inletRegions.append(args.optionReadList<word>("inletRegions"));
    }

    wordList outletRegions;
    if (args.optionFound("outletRegions"))
    {
        outletRegions.append(args.optionReadList<word>("outletRegions"));
    }

    const bool closedDomain(args.optionFound("closedDomain"));

    meshingSurfaceList surfaces
    (
        runTime,
        surfaceNames,
        cellZoneNames,
        rotatingZoneNames,
        baffleNames,
        bb,
        inletRegions,
        outletRegions,
        closedDomain
    );

    const Vector<label> nCells
    (
        args.optionLookupOrDefault("nCells", Vector<label>::zero)
    );

    const label minDimCells
    (
        args.optionLookupOrDefault("minDimCells", 0)
    );

    const label refineFactor
    (
        args.optionLookupOrDefault("refineBackground", 1)
    );

    HashTable<Pair<word>> patchOpts(7);
    patches.append("default");
    forAll(patches, i)
    {
        readPatchOption(args, patchOpts, patches[i] + "Patch");
    }

    const bool clearBoundary(args.optionFound("clearBoundary"));

    if (!args.optionFound("noBackground"))
    {
        if (args.optionFound("cylindricalBackground"))
        {
            blockMeshCylindricalConfiguration blockMeshConfig
            (
                "blockMeshDict",
                dir,
                runTime,
                surfaces,
                args.optionFound("bounds"),
                nCells,
                refineFactor,
                patchOpts,
                clearBoundary
            );

            blockMeshConfig.write();
        }
        else
        {
            blockMeshCartesianConfiguration blockMeshConfig
            (
                "blockMeshDict",
                dir,
                runTime,
                surfaces,
                args.optionFound("bounds"),
                nCells,
                minDimCells,
                refineFactor,
                patchOpts,
                clearBoundary
            );

            blockMeshConfig.write();
        }
    }

    // snappyHexMeshDict options
    const label refinementLevel
    (
        args.optionLookupOrDefault<label>("refinementLevel", 2)
    );

    List<Tuple2<word, label>> surfaceLevels;
    if (args.optionFound("surfaceLevels"))
    {
        surfaceLevels.append
        (
            args.optionReadList<Tuple2<word, label>>("surfaceLevels")
        );
    }

    List<Tuple2<word, label>> refinementRegions;
    if (args.optionFound("refinementRegions"))
    {
        refinementRegions.append
        (
            args.optionReadList<Tuple2<word, label>>("refinementRegions")
        );
    }

    List<Tuple3<vector, vector, label>> refinementBoxes;
    if (args.optionFound("refinementBoxes"))
    {
        refinementBoxes.append
        (
            args.optionReadList<Tuple3<vector, vector, label>>
            (
                "refinementBoxes"
            )
        );
    }

    List<Tuple3<word, scalar, label>> refinementDists;
    if (args.optionFound("refinementDists"))
    {
        refinementDists.append
        (
            args.optionReadList<Tuple3<word, scalar, label>>("refinementDists")
        );
    }

    const bool explicitFeatures(args.optionFound("explicitFeatures"));

    List<Tuple2<word, label>> layers;
    if (args.optionFound("layers"))
    {
        layers.append
        (
            args.optionReadList<Tuple2<word, label>>("layers")
        );
    }

    const scalar firstLayerThickness
    (
        args.optionLookupOrDefault<scalar>("firstLayerThickness", 0)
    );

    const scalar layerExpansionRatio
    (
        args.optionLookupOrDefault<scalar>("layerExpansionRatio", 1.2)
    );

    if (args.optionFound("insidePoint") && args.optionFound("insidePoints"))
    {
        FatalErrorInFunction
            << "Options '-insidePoint' and '-insidePoints' "
            << "cannot both be selected"
            << exit(FatalError);
    }

    List<point> insidePoints;
    bool insidePointsOpt(false);

    if (args.optionFound("insidePoints"))
    {
        insidePoints.append
        (
            args.optionReadList<point>("insidePoints")
        );

        insidePointsOpt = true;
    }
    else
    {
        insidePoints.append
        (
            args.optionLookupOrDefault<point>("insidePoint", point::zero)
        );
    }

    const label nCellsBetweenLevels
    (
        args.optionLookupOrDefault<label>("nCellsBetweenLevels", 3)
    );

    if (explicitFeatures)
    {
        surfaceFeaturesConfiguration surfaceFeaturesConfig
        (
            "surfaceFeaturesDict",
            dir,
            runTime,
            surfaces
        );

        surfaceFeaturesConfig.write();
    }

    snappyHexMeshConfiguration snappyConfig
    (
        "snappyHexMeshDict",
        dir,
        runTime,
        surfaces,
        refinementLevel,
        surfaceLevels,
        refinementRegions,
        refinementBoxes,
        refinementDists,
        explicitFeatures,
        layers,
        firstLayerThickness,
        layerExpansionRatio,
        insidePointsOpt,
        insidePoints,
        nCellsBetweenLevels
    );

    snappyConfig.write();

    meshQualityConfiguration meshQualityConfig
    (
        "meshQualityDict",
        dir,
        runTime
    );

    meshQualityConfig.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
