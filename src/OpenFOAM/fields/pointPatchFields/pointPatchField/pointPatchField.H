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
    Foam::pointPatchField

Description
    Abstract base class for point-mesh patch fields.

    The base-field does not store values as they are part of the
    "internal field".  There are derived classes to store constraint values
    e.g. fixedValuePointPatchField derived from the generic
    valuePointPatchField which ensures the values in the "internal field"
    are reset to the fixed-values by applying the stored values.

SourceFiles
    pointPatchField.C
    newPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatchField_H
#define pointPatchField_H

#include "pointPatch.H"
#include "DimensionedField.H"
#include "pointPatchFieldMapperFwd.H"
#include "fieldTypes.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class objectRegistry;
class dictionary;
class fieldMapper;
class pointMesh;

// Forward declaration of friend functions and operators

template<class Type>
class pointPatchField;

template<class Type>
class calculatedPointPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const pointPatchField<Type>&);


/*---------------------------------------------------------------------------*\
                       Class pointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointPatchField
{
    // Private Data

        //- Reference to patch
        const pointPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, pointMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;


public:

    typedef Type value_type;
    typedef pointPatch Patch;
    typedef calculatedPointPatchField<Type> Calculated;


    //- Runtime type information
    TypeName("pointPatchField");

    //- Debug switch to disallow the use of genericPointPatchField
    static int disallowGenericPointPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchField,
            pointPatch,
            (
                const pointPatch& p,
                const DimensionedField<Type, pointMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchField,
            patchMapper,
            (
                const pointPatchField<Type>& ptf,
                const pointPatch& p,
                const DimensionedField<Type, pointMesh>& iF,
                const fieldMapper& m
            ),
            (dynamic_cast<const pointPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchField,
            dictionary,
            (
                const pointPatch& p,
                const DimensionedField<Type, pointMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        pointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pointPatchField<Type> onto a new patch
        pointPatchField
        (
            const pointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        pointPatchField(const pointPatchField<Type>&) = delete;

        //- Disallow clone without setting internal field reference
        autoPtr<pointPatchField<Type>> clone() const
        {
            NotImplemented;
            return autoPtr<pointPatchField<Type>>(nullptr);
        }

        //- Construct as copy setting internal field reference
        pointPatchField
        (
            const pointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const = 0;

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const pointPatch& patch,
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            NotImplemented;
            return autoPtr<pointPatchField<Type>>(nullptr);
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static autoPtr<pointPatchField<Type>> New
        (
            const word&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static autoPtr<pointPatchField<Type>> New
        (
            const word&,
            const word& actualPatchType,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given pointPatchField mapped onto a new patch
        static autoPtr<pointPatchField<Type>> New
        (
            const pointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const fieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static autoPtr<pointPatchField<Type>> New
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedPointPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static autoPtr<pointPatchField<Type>>
        NewCalculatedType
        (
            const pointPatchField<Type2>&
        );


    //- Destructor
    virtual ~pointPatchField()
    {}


    // Member Functions

        // Attributes

            //- Return the type of the calculated form of pointPatchField
            static const word& calculatedType();

            //- Return size
            label size() const
            {
                return patch().size();
            }

            //- Return true if this patch field fixes a value
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- Return true if this overrides the underlying constraint type
            bool overridesConstraint() const
            {
                if (type() == patch_.type())
                {
                    return false;
                }

                typename pointPatchConstructorTable::iterator patchTypeCstrIter
                    = pointPatchConstructorTablePtr_->find(patch_.type());

                return
                    patchTypeCstrIter
                 != pointPatchConstructorTablePtr_->end();
            }


        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const pointPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, pointMesh>& internalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& primitiveField() const
            {
                return internalField_;
            }

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }

            //- Return field created from appropriate internal field values
            tmp<Field<Type>> patchInternalField() const;

            //- Return field created from appropriate internal field values
            //  given internal field reference
            template<class Type1>
            tmp<Field<Type1>> patchInternalField
            (
                const Field<Type1>& iF
            ) const;

            //- Return field created from selected internal field values
            //  given internal field reference
            template<class Type1>
            tmp<Field<Type1>> patchInternalField
            (
                const Field<Type1>& iF,
                const labelList& meshPoints
            ) const;

            //- Given the internal field and a patch field,
            //  add the patch field to the internal field
            template<class Type1>
            void addToInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;

            //- Given the internal field and a patch field,
            //  add selected elements of the patch field to the internal field
            template<class Type1>
            void addToInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF,
                const labelList& points
            ) const;

            //- Given the internal field and a patch field,
            //  set the patch field in the internal field
            template<class Type1>
            void setInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF,
                const labelList& meshPoints
            ) const;

            //- Given the internal field and a patch field,
            //  set the patch field in the internal field
            template<class Type1>
            void setInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;


        // Mapping functions

            //- Map the given pointPatchField onto this pointPatchField
            virtual void map(const pointPatchField<Type>&, const fieldMapper&)
            {}

            //- Reset the pointPatchField to the given pointPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const pointPatchField<Type>&)
            {}


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs()
            {
                updated_ = true;
            }

            //- Initialise evaluation of the patch field (do nothing)
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {}

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const pointPatchField<Type>&) {}
        virtual void operator=(pointPatchField<Type>&&) {}
        virtual void operator+=(const pointPatchField<Type>&) {}
        virtual void operator-=(const pointPatchField<Type>&) {}
        virtual void operator*=(const pointPatchField<scalar>&) {}
        virtual void operator/=(const pointPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}


        // Force an assignment irrespective of form of patch. These do nothing
        // unless the patch actually has boundary values.

        virtual void operator==(const pointPatchField<Type>&) {}
        virtual void operator==(const Field<Type>&) {}
        virtual void operator==(const Type&) {}


    // Ostream operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const pointPatchField<Type>&
        );
};


// This function is added as a hack to enable simple backward compatibility
// with versions using referenceLevel in GeometicField
template<class Type>
const pointPatchField<Type>& operator+
(
    const pointPatchField<Type>& ppf,
    const Type&
)
{
    return ppf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointPatchFieldFunctions.H"

#ifdef NoRepository
    #include "pointPatchField.C"
    #include "calculatedPointPatchField.H"
#endif


#define makePointPatchField(pointPatchTypeField)                               \
    defineNamedTemplateTypeNameAndDebug(pointPatchTypeField, 0);               \
    template<>                                                                 \
    int pointPatchTypeField::disallowGenericPointPatchField                    \
    (                                                                          \
        debug::debugSwitch("disallowGenericPointPatchField", 0)                \
    );                                                                         \
    defineTemplateRunTimeSelectionTable(pointPatchTypeField, pointPatch);      \
    defineTemplateRunTimeSelectionTable(pointPatchTypeField, patchMapper);     \
    defineTemplateRunTimeSelectionTable(pointPatchTypeField, dictionary);


#define addToPointPatchFieldRunTimeSelection(PatchTypeField,typePatchTypeField)\
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


#define addToNullConstructablePointPatchFieldRunTimeSelection(                 \
    PatchTypeField, typePatchTypeField)                                        \
                                                                               \
    addToPointPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)   \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        pointPatch                                                             \
    );


// For non-templated patch fields
#define makePointPatchTypeField(PatchTypeField, typePatchTypeField)            \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToPointPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


// For non-templated patch fields
#define makeNullConstructablePointPatchTypeField(                              \
    PatchTypeField, typePatchTypeField)                                        \
                                                                               \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToNullConstructablePointPatchFieldRunTimeSelection                      \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField                                                     \
    )


// For templated patch fields
#define makeTemplatePointPatchTypeField(fieldType, type)                       \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, PointPatch, CAPITALIZE(fieldType), Field),                  \
        0                                                                      \
    );                                                                         \
    addToPointPatchFieldRunTimeSelection                                       \
    (                                                                          \
        CAT3(pointPatch, CAPITALIZE(fieldType), Field),                        \
        CAT4(type, PointPatch, CAPITALIZE(fieldType), Field)                   \
    )


// For templated patch fields
#define makeNullConstructableTemplatePointPatchTypeField(fieldType, type)      \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, PointPatch, CAPITALIZE(fieldType), Field),                  \
        0                                                                      \
    );                                                                         \
    addToNullConstructablePointPatchFieldRunTimeSelection                      \
    (                                                                          \
        CAT3(pointPatch, CAPITALIZE(fieldType), Field),                        \
        CAT4(type, PointPatch, CAPITALIZE(fieldType), Field)                   \
    )


#define makePointPatchFields(type)                                             \
    FOR_ALL_FIELD_TYPES(makeTemplatePointPatchTypeField, type)


#define makeNullConstructablePointPatchFields(type)                            \
    FOR_ALL_FIELD_TYPES(makeNullConstructableTemplatePointPatchTypeField, type)


#define makePointPatchFieldTypeName(fieldType, type)                           \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, PointPatch, CAPITALIZE(fieldType), Field),                  \
        0                                                                      \
    );


#define makePointPatchFieldTypeNames(type)                                     \
    FOR_ALL_FIELD_TYPES(makePointPatchFieldTypeName, type)


#define makePointPatchFieldTypedef(fieldType, type)                            \
    typedef type##PointPatchField<fieldType>                                   \
        CAT4(type, PointPatch, CAPITALIZE(fieldType), Field);


#define makePointPatchFieldTypedefs(type)                                      \
    FOR_ALL_FIELD_TYPES(makePointPatchFieldTypedef, type)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
