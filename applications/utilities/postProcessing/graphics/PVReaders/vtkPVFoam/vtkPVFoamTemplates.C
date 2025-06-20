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
#include "vtkOpenFOAMPoints.H"

// OpenFOAM includes
#include "polyPatch.H"
#include "primitivePatch.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PatchType>
vtkPolyData* Foam::vtkPVFoam::patchVTKMesh(const PatchType& p)
{
    vtkPolyData* vtkmesh = vtkPolyData::New();

    // Convert OpenFOAM mesh vertices to VTK
    const Foam::pointField& points = p.localPoints();

    vtkPoints* vtkpoints = vtkPoints::New();
    vtkpoints->Allocate(points.size());
    forAll(points, i)
    {
        vtkInsertNextOpenFOAMPoint(vtkpoints, points[i]);
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    // Add faces as polygons
    const faceList& faces = p.localFaces();

    vtkCellArray* vtkcells = vtkCellArray::New();
    vtkcells->Allocate(faces.size());

    DynamicList<vtkIdType> nodeIds(4);

    forAll(faces, facei)
    {
        const face& f = faces[facei];

        nodeIds.resize(f.size());
        forAll(f, fp)
        {
            nodeIds[fp] = f[fp];
        }

        vtkcells->InsertNextCell(f.size(), nodeIds.data());
    }

    vtkmesh->SetPolys(vtkcells);
    vtkcells->Delete();

    return vtkmesh;
}


// ************************************************************************* //
