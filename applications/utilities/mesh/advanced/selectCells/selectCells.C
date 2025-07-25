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

Application
    selectCells

Description
    Select cells in relation to surface.

    Divides cells into three sets:
    - cutCells : cells cut by surface or close to surface.
    - outside  : cells not in cutCells and reachable from set of
      user-defined points (outsidePoints)
    - inside   : same but not reachable.

    Finally the wanted sets are combined into a cellSet 'selected'. Apart
    from straightforward adding the contents there are a few extra rules to
    make sure that the surface of the 'outside' of the mesh is singly
    connected.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IOdictionary.H"
#include "twoDPointCorrector.H"
#include "OFstream.H"
#include "meshTools.H"

#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"
#include "cellClassification.H"
#include "cellSet.H"
#include "cellInfo.H"
#include "FaceCellWave.H"
#include "edgeStats.H"
#include "treeDataTriSurface.H"
#include "indexedOctree.H"
#include "globalMeshData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// cellType for cells included/not included in mesh.
static const label MESH = cellClassification::INSIDE;
static const label NONMESH = cellClassification::OUTSIDE;


void writeSet(const cellSet& cells, const string& msg)
{
    Info<< "Writing " << msg << " (" << cells.size() << ") to cellSet "
        << cells.instance()/cells.local()/cells.name()
        << endl << endl;
    cells.write();
}


void getType(const labelList& elems, const label type, labelHashSet& set)
{
    forAll(elems, i)
    {
        if (elems[i] == type)
        {
            set.insert(i);
        }
    }
}


void cutBySurface
(
    const polyMesh& mesh,
    const triSurfaceSearch& querySurf,

    const pointField& outsidePts,
    const bool selectCut,
    const bool selectInside,
    const bool selectOutside,
    const scalar nearDist,

    cellClassification& cellType
)
{
    // Cut with surface and classify as inside/outside/cut
    cellType = cellClassification(mesh, querySurf, outsidePts);

    // Get inside/outside/cutCells cellSets.
    cellSet inside(mesh, "inside", mesh.nCells()/10);
    getType(cellType, cellClassification::INSIDE, inside);
    writeSet(inside, "inside cells");

    cellSet outside(mesh, "outside", mesh.nCells()/10);
    getType(cellType, cellClassification::OUTSIDE, outside);
    writeSet(outside, "outside cells");

    cellSet cutCells(mesh, "cutCells", mesh.nCells()/10);
    getType(cellType, cellClassification::CUT, cutCells);
    writeSet(cutCells, "cells cut by surface");


    // Change cellType to reflect selected part of mesh. Use
    // MESH to denote selected part, NONMESH for all
    // other cells.
    // Is a bit of a hack but allows us to reuse all the functionality
    // in cellClassification.

    forAll(cellType, celli)
    {
        label cType = cellType[celli];

        if (cType == cellClassification::CUT)
        {
            if (selectCut)
            {
                cellType[celli] = MESH;
            }
            else
            {
                cellType[celli] = NONMESH;
            }
        }
        else if (cType == cellClassification::INSIDE)
        {
            if (selectInside)
            {
                cellType[celli] = MESH;
            }
            else
            {
                cellType[celli] = NONMESH;
            }
        }
        else if (cType == cellClassification::OUTSIDE)
        {
            if (selectOutside)
            {
                cellType[celli] = MESH;
            }
            else
            {
                cellType[celli] = NONMESH;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Multiple mesh regions in original mesh" << endl
                << "Please use splitMeshRegions to separate these"
                << exit(FatalError);
        }
    }


    if (nearDist > 0)
    {
        Info<< "Removing cells with points closer than " << nearDist
            << " to the surface ..." << nl << endl;

        const pointField& pts = mesh.points();
        const indexedOctree<treeDataTriSurface>& tree = querySurf.tree();

        label nRemoved = 0;

        forAll(pts, pointi)
        {
            const point& pt = pts[pointi];

            pointIndexHit hitInfo = tree.findNearest(pt, sqr(nearDist));

            if (hitInfo.hit())
            {
                const labelList& pCells = mesh.pointCells()[pointi];

                forAll(pCells, i)
                {
                    if (cellType[pCells[i]] != NONMESH)
                    {
                        cellType[pCells[i]] = NONMESH;
                        nRemoved++;
                    }
                }
            }
        }

//        tmp<pointField> tnearest = querySurf.calcNearest(pts);
//        const pointField& nearest = tnearest();
//
//        label nRemoved = 0;
//
//        forAll(nearest, pointi)
//        {
//            if (mag(nearest[pointi] - pts[pointi]) < nearDist)
//            {
//                const labelList& pCells = mesh.pointCells()[pointi];
//
//                forAll(pCells, i)
//                {
//                    if (cellType[pCells[i]] != NONMESH)
//                    {
//                        cellType[pCells[i]] = NONMESH;
//                        nRemoved++;
//                    }
//                }
//            }
//        }

        Info<< "Removed " << nRemoved << " cells since too close to surface"
            << nl << endl;
    }
}



// We're meshing the outside. Subset the currently selected mesh cells with the
// ones reachable from the outsidepoints.
label selectOutsideCells
(
    const polyMesh& mesh,
    const pointField& outsidePts,
    cellClassification& cellType
)
{
    //
    // Check all outsidePts and for all of them inside a mesh cell
    // collect the faces to start walking from
    //

    // Outside faces
    labelHashSet outsideFacesMap(outsidePts.size() * 6 * 2);
    DynamicList<label> outsideFaces(outsideFacesMap.size());
    // CellInfo on outside faces
    DynamicList<cellInfo> outsideFacesInfo(outsideFacesMap.size());

    // cellInfo for mesh cell
    const cellInfo meshInfo(MESH);

    forAll(outsidePts, outsidePtI)
    {
        // Find cell containing point. Linear search.
        label celli = meshSearch::findCellNoTree(mesh, outsidePts[outsidePtI]);

        if (celli != -1 && cellType[celli] == MESH)
        {
            Info<< "Marking cell " << celli << " containing outside point "
                << outsidePts[outsidePtI] << " with type " << cellType[celli]
                << " ..." << endl;

            //
            // Mark this cell and its faces to start walking from
            //

            // Mark faces of celli
            const labelList& cFaces = mesh.cells()[celli];
            forAll(cFaces, i)
            {
                label facei = cFaces[i];

                if (outsideFacesMap.insert(facei))
                {
                    outsideFaces.append(facei);
                    outsideFacesInfo.append(meshInfo);
                }
            }
        }
    }

    // Floodfill starting from outsideFaces (of type meshInfo)
    List<cellInfo> faceInfoList(mesh.nFaces()), cellInfoList(mesh.nCells());
    FaceCellWave<cellInfo> regionCalc
    (
        mesh,
        outsideFaces.shrink(),
        outsideFacesInfo.shrink(),
        faceInfoList,
        cellInfoList,
        mesh.globalData().nTotalCells() + 1 // max iterations
    );

    // Now cellInfoList should hold info on cells that are reachable from
    // changedFaces. Use these to subset cellType
    label nChanged = 0;
    forAll(cellInfoList, celli)
    {
        if (cellType[celli] == MESH)
        {
            // Original cell was selected for meshing. Check if cell was
            // reached from outsidePoints
            if (cellInfoList[celli].type() != MESH)
            {
                cellType[celli] = NONMESH;
                nChanged++;
            }
        }
    }

    return nChanged;
}



int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    // Mesh edge statistics calculator
    edgeStats edgeCalc(mesh);


    IOdictionary refineDict
    (
        IOobject
        (
            "selectCellsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    fileName surfName(refineDict.lookup("surface"));
    pointField outsidePts(refineDict.lookup("outsidePoints"));
    bool useSurface(readBool(refineDict.lookup("useSurface")));
    bool selectCut(readBool(refineDict.lookup("selectCut")));
    bool selectInside(readBool(refineDict.lookup("selectInside")));
    bool selectOutside(readBool(refineDict.lookup("selectOutside")));
    scalar nearDist(refineDict.lookup<scalar>("nearDistance"));


    if (useSurface)
    {
        Info<< "Cells to be used for meshing (0=false, 1=true):" << nl
            << "    cells cut by surface            : " << selectCut << nl
            << "    cells inside of surface         : " << selectInside << nl
            << "    cells outside of surface        : " << selectOutside << nl
            << "    cells with points further than  : " << nearDist << nl
            << endl;
    }
    else
    {
        Info<< "Cells to be used for meshing (0=false, 1=true):" << nl
            << "    cells reachable from outsidePoints:" << selectOutside << nl
            << endl;
    }

    // Print edge stats on original mesh.
    (void)edgeCalc.minLen(Info);

    // Check all 'outside' points
    forAll(outsidePts, outsideI)
    {
        const point& outsidePoint = outsidePts[outsideI];

        label celli = meshSearch::findCellNoTree(mesh, outsidePoint);

        if (returnReduce(celli, maxOp<label>()) == -1)
        {
            FatalErrorInFunction
                << "outsidePoint " << outsidePoint
                << " is not inside any cell"
                << exit(FatalError);
        }
    }

    // Cell status (compared to surface if provided): inside/outside/cut.
    // Start off from everything selected and cut later.
    cellClassification cellType
    (
        mesh,
        labelList
        (
            mesh.nCells(),
            cellClassification::MESH
        )
    );


    // Surface
    autoPtr<triSurface> surf(nullptr);
    // Search engine on surface.
    autoPtr<triSurfaceSearch> querySurf(nullptr);

    if (useSurface)
    {
        surf.reset(new triSurface(surfName));

        // Dump some stats
        surf().writeStats(Info);

        // Search engine on surface.
        querySurf.reset(new triSurfaceSearch(surf));

        // Set cellType[celli] according to relation to surface
        cutBySurface
        (
            mesh,
            querySurf,

            outsidePts,
            selectCut,
            selectInside,
            selectOutside,
            nearDist,

            cellType
        );
    }


    // Now 'trim' all the corners from the mesh so meshing/surface extraction
    // becomes easier.

    label nHanging, nRegionEdges, nRegionPoints, nOutside;

    do
    {
        Info<< "Removing cells which after subsetting would have all points"
            << " on outside ..." << nl << endl;

        nHanging = cellType.fillHangingCells
        (
            MESH,       // meshType
            NONMESH,    // fill type
            mesh.nCells()
        );


        Info<< "Removing edges connecting cells unconnected by faces ..."
            << nl << endl;

        nRegionEdges = cellType.fillRegionEdges
        (
            MESH,       // meshType
            NONMESH,    // fill type
            mesh.nCells()
        );


        Info<< "Removing points connecting cells unconnected by faces ..."
            << nl << endl;

        nRegionPoints = cellType.fillRegionPoints
        (
            MESH,       // meshType
            NONMESH,    // fill type
            mesh.nCells()
        );

        nOutside = 0;
        if (selectOutside)
        {
            // Since we're selecting the cells reachable from outsidePoints
            // and the set might have changed, redo the outsideCells
            // calculation
            nOutside = selectOutsideCells(mesh, outsidePts, cellType);
        }
    } while
    (
        nHanging != 0
     || nRegionEdges != 0
     || nRegionPoints != 0
     || nOutside != 0
    );

    cellSet selectedCells(mesh, "selected", mesh.nCells()/10);
    getType(cellType, MESH, selectedCells);

    writeSet(selectedCells, "cells selected for meshing");


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
