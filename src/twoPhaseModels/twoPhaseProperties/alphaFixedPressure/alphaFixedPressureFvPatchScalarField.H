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

Class
    Foam::alphaFixedPressureFvPatchScalarField

Description
    A fixed-pressure alphaContactAngle boundary

SourceFiles
    alphaFixedPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaFixedPressureFvPatchScalarField_H
#define alphaFixedPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class alphaFixedPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class alphaFixedPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Fixed pressure
        scalarField p_;


public:

    //- Runtime type information
    TypeName("alphaFixedPressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        alphaFixedPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphaFixedPressureFvPatchScalarField
        //  onto a new patch
        alphaFixedPressureFvPatchScalarField
        (
            const alphaFixedPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        alphaFixedPressureFvPatchScalarField
        (
            const alphaFixedPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        alphaFixedPressureFvPatchScalarField
        (
            const alphaFixedPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new alphaFixedPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the alphaFixed pressure
            const scalarField& p() const
            {
                return p_;
            }

            //- Return reference to the alphaFixed pressure to allow adjustment
            scalarField& p()
            {
                return p_;
            }


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
