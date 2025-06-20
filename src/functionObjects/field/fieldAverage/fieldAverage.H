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
    Foam::functionObjects::fieldAverage

Description
    Calculates average quantities for a user-specified selection of volumetric
    and surface fields.

    Fields are entered as a list of sub-dictionaries, which indicate the type of
    averages to perform, and can be updated during the calculation.  The current
    options include:
    - \c mean: arithmetic mean
    - \c prime2Mean: prime-squared mean
    - \c base: average over 'time', or 'iteration'
    - \c window: optional averaging window, specified in 'base' units

    Average field names are constructed by concatenating the base field with
    the averaging type, e.g. when averaging field 'U', the resultant fields
    are:
    - arithmetic mean field, \c UMean
    - prime-squared field, \c UPrime2Mean

    Information regarding the number of averaging steps, and total averaging
    time are written on a per-field basis to the \c "<functionObject
    name>Properties" dictionary, located in \<time\>/uniform

    When restarting form a previous calculation, the averaging is continuous or
    may be restarted using the \c restartOnRestart option.

    The averaging process may be restarted after each calculation output time
    using the \c restartOnOutput option or restarted periodically using the \c
    periodicRestart option and setting \c restartPeriod to the required
    averaging period.

    Example of function object specification:
    \verbatim
    fieldAverage1
    {
        type                fieldAverage;
        libs                ("libfieldFunctionObjects.so");

        writeControl        writeTime;

        restartOnRestart    false;
        restartOnOutput     false;
        periodicRestart     false;
        restartPeriod       0.002;

        base                time;
        window              10.0;
        windowName          w1;

        mean                yes;
        prime2Mean          yes;

        fields              (U p);
    }
    \endverbatim

Usage
    \table
        Property          | Description               | Required | Default
        type              | type name: fieldAverage              | yes |
        restartOnRestart  | Restart the averaging on restart     | no  | no
        restartOnOutput   | Restart the averaging on output      | no  | no
        periodicRestart   | Periodically restart the averaging   | no  | no
        restartPeriod     | Periodic restart period              | conditional |
        fields            | list of fields and averaging options | yes |
    \endtable

    Note:
        To employ the \c prime2Mean option, the \c mean option must be selected.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    fieldAverage.C
    fieldAverageTemplates.C
    fieldAverageItem.C

\*---------------------------------------------------------------------------*/

#ifndef fieldAverage_functionObject_H
#define fieldAverage_functionObject_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

// Forward declaration of classes
class fieldAverageItem;

/*---------------------------------------------------------------------------*\
                        Class fieldAverage Declaration
\*---------------------------------------------------------------------------*/

class fieldAverage
:
    public fvMeshFunctionObject
{
public:

        //- Enumeration defining the averaging base type
        enum class baseType
        {
            iter,
            time
        };


protected:

    // Protected data

        //- Time at last call, prevents repeated averaging
        label prevTimeIndex_;

        //- Restart the averaging process on restart
        Switch restartOnRestart_;

        //- Restart the averaging process on output
        Switch restartOnOutput_;

        //- Periodically restart the averaging process
        Switch periodicRestart_;

        //- Restart period
        scalar restartPeriod_;

        //- Averaging base type names
        static const NamedEnum<baseType, 2> baseTypeNames_;

        //- Averaging base type
        baseType base_;

        //- Averaging window - defaults to -1 for 'all iters/time'
        scalar window_;

        //- Averaging window name - defaults to 'window'
        word windowName_;

        //- Compute mean flag
        Switch mean_;

        //- Compute prime-squared mean flag
        Switch prime2Mean_;

        //- List of field average items, describing what averages to be
        //  calculated and output
        PtrList<fieldAverageItem> faItems_;

        //- Iteration steps counter
        List<label> totalIter_;

        //- Total time counter
        List<scalar> totalTime_;


    // Protected Member Functions

        // Initialisation routines

            //- Reset lists (clear existing values) and initialise averaging.
            //  Check requested field averages are valid, populate field lists
            void initialise();

            //- Restart averaging for restartOnOutput
            void restart();

            //- Read the mean average field
            template<class Type>
            void readMeanFieldType(const label fieldi);

            //- Read the mean average field
            template<class Type>
            void readMeanField(const label fieldi);

            //- Initialise the mean average field
            template<class Type>
            void initialiseMeanFieldType(const label fieldi);

            //- Initialise the mean average field
            template<class Type>
            void initialiseMeanField(const label fieldi);

            //- Read the prime-squared average field
            template<class Type1, class Type2>
            void readPrime2MeanFieldType(const label fieldi);

            //- Read the prime-squared average field
            template<class Type1, class Type2>
            void readPrime2MeanField(const label fieldi);

            //- Initialise the prime-squared average field
            template<class Type1, class Type2>
            void initialisePrime2MeanFieldType(const label fieldi);

            //- Initialise the prime-squared average field
            template<class Type1, class Type2>
            void initialisePrime2MeanField(const label fieldi);


        // Calculation functions

            //- Main calculation routine
            virtual void calcAverages();

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFieldType(const label fieldi) const;

            //- Calculate mean average fields
            template<class Type>
            void calculateMeanFields() const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFieldType(const label fieldi) const;

            //- Calculate prime-squared average fields
            template<class Type1, class Type2>
            void calculatePrime2MeanFields() const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2MeanType(const label fieldi) const;

            //- Add mean-squared field value to prime-squared mean field
            template<class Type1, class Type2>
            void addMeanSqrToPrime2Mean() const;


        // I-O

            //- Write averages
            virtual void writeAverages() const;

            //- Write fields
            template<class Type>
            void writeFieldType(const word& fieldName) const;

            //- Write fields
            template<class Type>
            void writeFields() const;

            //- Read
            void read(const dictionary& dict, const bool construct);


public:

    friend class fieldAverageItem;


    //- Runtime type information
    TypeName("fieldAverage");


    // Constructors

        //- Construct from Time and dictionary
        fieldAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        fieldAverage(const fieldAverage&) = delete;


    //- Destructor
    virtual ~fieldAverage();


    // Member Functions

        // Access

            //- Return averaging base type name
            const word base() const
            {
                return baseTypeNames_[base_];
            }

            //- Return true if base is iter
            Switch iterBase() const
            {
                return base_ == baseType::iter;
            }

            //- Return true if base is time
            Switch timeBase() const
            {
                return base_ == baseType::time;
            }

            scalar window() const
            {
                return window_;
            }

            const word& windowName() const
            {
                return windowName_;
            }


        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Do not average at the start of the run
        virtual bool executeAtStart() const
        {
            return false;
        }

        //- Calculate the field averages
        virtual bool execute();

        //- Write the field averages
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldAverage&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
