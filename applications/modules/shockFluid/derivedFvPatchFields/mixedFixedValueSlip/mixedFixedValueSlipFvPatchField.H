/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::mixedFixedValueSlipFvPatchField

Description
    A mixed boundary type that blends between fixedValue and slip, as opposed
    to the standard mixed condition that blends between fixedValue and
    fixedGradient; required to implement maxwellSlipU condition.

SourceFiles
    mixedFixedValueSlipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mixedFixedValueSlipFvPatchField_H
#define mixedFixedValueSlipFvPatchField_H

#include "transformFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class mixedFixedValueSlipFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixedFixedValueSlipFvPatchField
:
    public transformFvPatchField<Type>
{
    // Private Data

        //- Value field used for boundary condition
        Field<Type> refValue_;

        //- Fraction (0-1) of value used for boundary condition
        scalarField valueFraction_;


protected:

    // Protected functions

        //- Construct from patch and internal field
        mixedFixedValueSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );


public:

    //- Runtime type information
    TypeName("mixedFixedValueSlip");


    // Constructors

        //- Construct from patch, internal field and dictionary
        mixedFixedValueSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedFixedValueSlipFvPatchField
        //  onto a new patch
        mixedFixedValueSlipFvPatchField
        (
            const mixedFixedValueSlipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mixedFixedValueSlipFvPatchField
        (
            const mixedFixedValueSlipFvPatchField<Type>&
        ) = delete;

        //- Copy constructor setting internal field reference
        mixedFixedValueSlipFvPatchField
        (
            const mixedFixedValueSlipFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new mixedFixedValueSlipFvPatchField<Type>(*this, iF)
            );
        }

    // Member Functions

        // Access

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchField<Type>&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchField<Type>&);


        // Return defining fields

            virtual Field<Type>& refValue()
            {
                return refValue_;
            }

            virtual const Field<Type>& refValue() const
            {
                return refValue_;
            }

            virtual scalarField& valueFraction()
            {
                return valueFraction_;
            }

            virtual const scalarField& valueFraction() const
            {
                return valueFraction_;
            }

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=
                    Pstream::commsTypes::blocking
            );

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixedFixedValueSlipFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
