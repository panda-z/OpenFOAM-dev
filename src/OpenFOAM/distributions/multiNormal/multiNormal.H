/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::distributions::multiNormal

Description
    Multiple superimposed normal distributions

    \f[
        PDF(x) = \sum_i S_i \frac{1}{\sigma_i \sqrt{2 \pi}} \\
        \exp \left( \frac{1}{2} \left( \frac{x - \mu_i}{\sigma_i} \\
        \right)^2 \right)
    \f]

Usage
    Example usage:
    \verbatim
    {
        type    multiNormal;
        Q       0;
        min     0.001;
        max     0.019;
        mu      (0.005 0.011 0.015);
        sigma   (0.002 0.001 0.0015);
        strength (2 1 3);
    }
    \endverbatim

SourceFiles
    multiNormal.C

See also
    Foam::distribution
    Foam::distributions::normal

\*---------------------------------------------------------------------------*/

#ifndef multiNormal_H
#define multiNormal_H

#include "normal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{

/*---------------------------------------------------------------------------*\
                        Class multiNormal Declaration
\*---------------------------------------------------------------------------*/

class multiNormal
:
    public FieldDistribution<unintegrableForNonZeroQ, multiNormal>
{
    // Private Data

        //- Relative strengths of the different distributions
        const scalarList cumulativeStrengths_;

        //- Normal distributions
        PtrList<normal> distributions_;


    // Private Member Functions

        //- Normalise the strengths
        static scalarList readCumulativeStrengths(const dictionary& dict);

        //- Return values of the un-normalised PDF for the given size exponent
        //  and x-coordinates.
        virtual tmp<scalarField> phi
        (
            const label q,
            const scalarField& x
        ) const;

        //- Return values of the un-normalised CDF for zero effective size
        //  exponent and given x-coordinates
        virtual tmp<scalarField> PhiForZeroQ(const scalarField& x) const;


public:

    //- Runtime type information
    TypeName("multiNormal");


    // Constructors

        //- Construct from a dictionary
        multiNormal
        (
            const unitConversion& units,
            const dictionary& dict,
            const label sampleQ,
            randomGenerator&& rndGen
        );

        //- Construct copy
        multiNormal(const multiNormal& d, const label sampleQ);

        //- Construct and return a clone
        virtual autoPtr<distribution> clone(const label sampleQ) const
        {
            return autoPtr<distribution>(new multiNormal(*this, sampleQ));
        }


    //- Destructor
    virtual ~multiNormal();


    // Member Functions

        //- Sample the distribution for zero effective size exponent
        virtual scalar sampleForZeroQ() const;

        //- Sample the distribution
        using FieldDistribution<unintegrableForNonZeroQ, multiNormal>::sample;

        //- Return the minimum value
        virtual scalar min() const;

        //- Return the maximum value
        virtual scalar max() const;

        //- Write stream
        virtual void write(Ostream& os, const unitConversion& units) const;

        //- Return coordinates to plot across the range of the distribution
        virtual tmp<scalarField> plotX(const label n) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
