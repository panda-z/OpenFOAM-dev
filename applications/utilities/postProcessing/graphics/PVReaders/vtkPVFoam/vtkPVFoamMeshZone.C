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
#include "fvMesh.H"

// VTK includes
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkPolyData* Foam::vtkPVFoam::pointZoneVTKMesh
(
    const fvMesh& mesh,
    const labelList& pointLabels
)
{
    vtkPolyData* vtkmesh = vtkPolyData::New();

    const pointField& meshPoints = mesh.points();

    vtkPoints* vtkpoints = vtkPoints::New();
    vtkpoints->Allocate(pointLabels.size());

    forAll(pointLabels, pointi)
    {
        vtkInsertNextOpenFOAMPoint(vtkpoints, meshPoints[pointLabels[pointi]]);
    }

    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    return vtkmesh;
}


// ************************************************************************* //
