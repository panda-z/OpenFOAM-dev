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

#include "forcesBase.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "incompressibleMomentumTransportModel.H"
#include "compressibleMomentumTransportModel.H"
#include "phaseIncompressibleMomentumTransportModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "fluidThermo.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forcesBase, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::forcesBase::createFileNames
(
    const dictionary& dict
) const
{
    DynamicList<word> names(1);

    const word forceType(dict.lookup("type"));

    // Name for file(fileID::mainFile=0)
    names.append(forceType);

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));
        const label nb = binDict.lookup<label>("nBin");
        if (nb > 0)
        {
            // Name for file(fileID::binsFile=1)
            names.append(forceType + "_bins");
        }
    }

    return names;
}


void Foam::functionObjects::forcesBase::writeFileHeader(const label i)
{
    const word forceTypes
    (
        porosity_
      ? "(pressure viscous porous)"
      : "(pressure viscous)"
    );

    switch (fileID(i))
    {
        case fileID::mainFile:
        {
            // force data

            writeHeader(file(i), "Forces");
            writeCoRValueHeader(file(i));
            writeCommented(file(i), "Time");

            writeCoRHeader(file(i));

            file(i)
                << "forces" << forceTypes << tab
                << "moments" << forceTypes;

            break;
        }
        case fileID::binsFile:
        {
            // bin data

            writeHeader(file(i), "Force bins");
            writeHeaderValue(file(i), "bins", nBin_);
            writeHeaderValue(file(i), "start", binMin_);
            writeHeaderValue(file(i), "delta", binDx_);
            writeHeaderValue(file(i), "direction", binDir_);

            vectorField binPoints(nBin_);
            writeCommented(file(i), "x co-ords  :");
            forAll(binPoints, pointi)
            {
                binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
                file(i) << tab << binPoints[pointi].x();
            }
            file(i) << nl;

            writeCommented(file(i), "y co-ords  :");
            forAll(binPoints, pointi)
            {
                file(i) << tab << binPoints[pointi].y();
            }
            file(i) << nl;

            writeCommented(file(i), "z co-ords  :");
            forAll(binPoints, pointi)
            {
                file(i) << tab << binPoints[pointi].z();
            }
            file(i) << nl;

            writeCommented(file(i), "Time");

            for (label j = 0; j < nBin_; j++)
            {
                const word jn('(' + Foam::name(j) + ')');
                const word f("forces" + jn + forceTypes);
                const word m("moments" + jn + forceTypes);

                file(i)<< tab << f << tab << m;
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }

    file(i)<< endl;
}


void Foam::functionObjects::forcesBase::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!obr_.foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database."
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !obr_.foundObject<volVectorField>(UName_)
         || !obr_.foundObject<volScalarField>(pName_)

        )
        {
            FatalErrorInFunction
                << "Could not find " << UName_ << ", " << pName_
                << exit(FatalError);
        }

        if
        (
            rhoName_ != "rhoInf"
         && !obr_.foundObject<volScalarField>(rhoName_)
        )
        {
            FatalErrorInFunction
                << "Could not find " << rhoName_
                << exit(FatalError);
        }
    }

    initialised_ = true;
}




Foam::scalar Foam::functionObjects::forcesBase::rho
(
    const volScalarField& p
) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1.0;
    }
    else
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorInFunction
                << "Dynamic pressure is expected but kinematic is provided."
                << exit(FatalError);
        }

        return rhoRef_;
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forcesBase::alpha() const
{
    if (phaseName_ == word::null)
    {
        return volScalarField::New
        (
            "alpha",
            mesh_,
            dimensionedScalar(dimless, 1)
        );
    }
    else
    {
        return obr_.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phaseName_)
        );
    }
}


Foam::tmp<Foam::scalarField> Foam::functionObjects::forcesBase::alpha
(
    const label patchi
) const
{
    if (phaseName_ == word::null)
    {
        return tmp<scalarField>
        (
            new scalarField(mesh_.boundary()[patchi].size(), 1)
        );
    }
    else
    {
        return obr_.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phaseName_)
        ).boundaryField()[patchi];
    }
}


void Foam::functionObjects::forcesBase::applyBins
(
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP,
    const vectorField& d
)
{
    if (nBin_ == 1)
    {
        force_[0][0] += sum(fN);
        force_[1][0] += sum(fT);
        force_[2][0] += sum(fP);
        moment_[0][0] += sum(Md^fN);
        moment_[1][0] += sum(Md^fT);
        moment_[2][0] += sum(Md^fP);
    }
    else
    {
        scalarField dd((d & binDir_) - binMin_);

        forAll(dd, i)
        {
            label bini = min(max(floor(dd[i]/binDx_), 0), force_[0].size() - 1);

            force_[0][bini] += fN[i];
            force_[1][bini] += fT[i];
            force_[2][bini] += fP[i];
            moment_[0][bini] += Md[i]^fN[i];
            moment_[1][bini] += Md[i]^fT[i];
            moment_[2][bini] += Md[i]^fP[i];
        }
    }
}


void Foam::functionObjects::forcesBase::writeCoRValueHeader(Ostream& file)
{}


void Foam::functionObjects::forcesBase::writeCoRHeader(Ostream& file)
{
    file << "CofR" << tab;
}


void Foam::functionObjects::forcesBase::writeCofR(Ostream& file)
{
    file << CofR();
}


void Foam::functionObjects::forcesBase::writeForces()
{
    Log << type() << " " << name() << " write:" << nl
        << "    sum of forces:" << nl
        << "        pressure : " << sum(force_[0]) << nl
        << "        viscous  : " << sum(force_[1]) << nl
        << "        porous   : " << sum(force_[2]) << nl
        << "    sum of moments:" << nl
        << "        pressure : " << sum(moment_[0]) << nl
        << "        viscous  : " << sum(moment_[1]) << nl
        << "        porous   : " << sum(moment_[2])
        << endl;

    writeTime(file(fileID::mainFile));
    writeCofR(file(fileID::mainFile));

    if (porosity_)
    {
        file(fileID::mainFile) << tab << setw(1) << '('
            << sum(force_[0]) << setw(1) << ' '
            << sum(force_[1]) << setw(1) << ' '
            << sum(force_[2]) << setw(3) << ") ("
            << sum(moment_[0]) << setw(1) << ' '
            << sum(moment_[1]) << setw(1) << ' '
            << sum(moment_[2]) << setw(1) << ')';
    }
    else
    {
        file(fileID::mainFile) << tab << setw(1) << '('
            << sum(force_[0]) << setw(1) << ' '
            << sum(force_[1]) << setw(3) << ") ("
            << sum(moment_[0]) << setw(1) << ' '
            << sum(moment_[1]) << setw(1) << ')';
    }

    file(fileID::mainFile) << endl;
}


void Foam::functionObjects::forcesBase::writeBins()
{
    if (nBin_ == 1)
    {
        return;
    }

    List<Field<vector>> f(force_);
    List<Field<vector>> m(moment_);

    if (binCumulative_)
    {
        for (label i = 1; i < f[0].size(); i++)
        {
            f[0][i] += f[0][i-1];
            f[1][i] += f[1][i-1];
            f[2][i] += f[2][i-1];

            m[0][i] += m[0][i-1];
            m[1][i] += m[1][i-1];
            m[2][i] += m[2][i-1];
        }
    }

    writeTime(file(fileID::binsFile));


    forAll(f[0], i)
    {
        if (porosity_)
        {
            file(fileID::binsFile)
                << tab << setw(1) << '('
                << f[0][i] << setw(1) << ' '
                << f[1][i] << setw(1) << ' '
                << f[2][i] << setw(3) << ") ("
                << m[0][i] << setw(1) << ' '
                << m[1][i] << setw(1) << ' '
                << m[2][i] << setw(1) << ')';
        }
        else
        {
            file(fileID::binsFile)
                << tab << setw(1) << '('
                << f[0][i] << setw(1) << ' '
                << f[1][i] << setw(3) << ") ("
                << m[0][i] << setw(1) << ' '
                << m[1][i] << setw(1) << ')';
        }
    }

    file(fileID::binsFile) << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forcesBase::forcesBase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    force_(3),
    moment_(3),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    phaseName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(vGreat),
    pRef_(0),
    porosity_(false),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(great),
    binPoints_(),
    binCumulative_(true),
    initialised_(false)
{
    read(dict);
}


Foam::functionObjects::forcesBase::forcesBase
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, obr),
    logFiles(obr_, name),
    force_(3),
    moment_(3),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    phaseName_(word::null),
    directForceDensity_(false),
    fDName_(""),
    rhoRef_(vGreat),
    pRef_(0),
    porosity_(false),
    nBin_(1),
    binDir_(Zero),
    binDx_(0.0),
    binMin_(great),
    binPoints_(),
    binCumulative_(true),
    initialised_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forcesBase::~forcesBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forcesBase::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    initialised_ = false;

    Log << type() << " " << name() << ":" << nl;

    directForceDensity_ = dict.lookupOrDefault("directForceDensity", false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ = mesh_.boundaryMesh().patchSet(dict);

    if (directForceDensity_)
    {
        // Optional entry for fDName
        fDName_ = dict.lookupOrDefault<word>("fD", "fD");
    }
    else
    {
        // Optional phase entry
        phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

        // Optional U, p and rho entries
        pName_ =
            dict.lookupOrDefault<word>
            (
                "p",
                IOobject::groupName("p", phaseName_)
            );
        UName_ =
            dict.lookupOrDefault<word>
            (
                "U",
                IOobject::groupName("U", phaseName_)
            );
        rhoName_ =
            dict.lookupOrDefault<word>
            (
                "rho",
                IOobject::groupName("rho", phaseName_)
            );

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            dict.lookup("rhoInf") >> rhoRef_;
        }

        // Reference pressure, 0 by default
        pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Log << "    Including porosity effects" << endl;
    }
    else
    {
        Log << "    Not including porosity effects" << endl;
    }

    if (dict.found("binData"))
    {
        const dictionary& binDict(dict.subDict("binData"));
        binDict.lookup("nBin") >> nBin_;

        if (nBin_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Number of bins (nBin) must be zero or greater"
                << exit(FatalIOError);
        }
        else if ((nBin_ == 0) || (nBin_ == 1))
        {
            nBin_ = 1;
            forAll(force_, i)
            {
                force_[i].setSize(1);
                moment_[i].setSize(1);
            }
        }

        if (nBin_ > 1)
        {
            binDict.lookup("direction") >> binDir_;
            binDir_ /= mag(binDir_);

            binMin_ = great;
            scalar binMax = -great;
            forAllConstIter(labelHashSet, patchSet_, iter)
            {
                const label patchi = iter.key();
                const polyPatch& pp = pbm[patchi];
                const scalarField d(pp.faceCentres() & binDir_);
                binMin_ = min(min(d), binMin_);
                binMax = max(max(d), binMax);
            }
            reduce(binMin_, minOp<scalar>());
            reduce(binMax, maxOp<scalar>());

            // slightly boost binMax so that region of interest is fully
            // within bounds
            binMax = 1.0001*(binMax - binMin_) + binMin_;

            binDx_ = (binMax - binMin_)/scalar(nBin_);

            // create the bin points used for writing
            binPoints_.setSize(nBin_);
            forAll(binPoints_, i)
            {
                binPoints_[i] = (i + 0.5)*binDir_*binDx_;
            }

            binDict.lookup("cumulative") >> binCumulative_;

            // allocate storage for forces and moments
            forAll(force_, i)
            {
                force_[i].setSize(nBin_);
                moment_[i].setSize(nBin_);
            }
        }
    }

    if (nBin_ == 1)
    {
        // allocate storage for forces and moments
        force_[0].setSize(1);
        force_[1].setSize(1);
        force_[2].setSize(1);
        moment_[0].setSize(1);
        moment_[1].setSize(1);
        moment_[2].setSize(1);
    }

    resetNames(createFileNames(dict));

    return true;
}

Foam::tmp<Foam::surfaceVectorField>
Foam::functionObjects::forcesBase::devTau() const
{
    typedef incompressible::momentumTransportModel icoModel;
    typedef compressible::momentumTransportModel cmpModel;
    typedef phaseIncompressible::momentumTransportModel phaseIcoModel;
    typedef phaseCompressible::momentumTransportModel phaseCmpModel;

    const word& modelName = momentumTransportModel::typeName;
    const word phaseModelName =
        phaseName_ == word::null
      ? word::null
      : IOobject::groupName(momentumTransportModel::typeName, phaseName_);

    if (obr_.foundObject<icoModel>(modelName))
    {
        const incompressible::momentumTransportModel& model =
            obr_.lookupObject<icoModel>(modelName);

        return fvc::interpolate(alpha()*rho())*model.devSigma();
    }
    else if (obr_.foundObject<cmpModel>(modelName))
    {
        const cmpModel& model =
            obr_.lookupObject<cmpModel>(modelName);

        return fvc::interpolate(alpha())*model.devTau();
    }
    else if (obr_.foundObject<phaseIcoModel>(phaseModelName))
    {
        const phaseIcoModel& model =
            obr_.lookupObject<phaseIcoModel>(phaseModelName);

        return fvc::interpolate(rho())*model.devSigma();
    }
    else if (obr_.foundObject<phaseCmpModel>(phaseModelName))
    {
        const phaseCmpModel& model =
            obr_.lookupObject<phaseCmpModel>(phaseModelName);

        return model.devTau();
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return surfaceVectorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forcesBase::mu() const
{
    typedef incompressible::momentumTransportModel icoModel;
    typedef compressible::momentumTransportModel cmpModel;
    typedef phaseIncompressible::momentumTransportModel phaseIcoModel;
    typedef phaseCompressible::momentumTransportModel phaseCmpModel;

    const word& modelName = momentumTransportModel::typeName;
    const word phaseModelName =
        phaseName_ == word::null
      ? word::null
      : IOobject::groupName(momentumTransportModel::typeName, phaseName_);

    if (obr_.foundObject<icoModel>(modelName))
    {
        const incompressible::momentumTransportModel& model =
            obr_.lookupObject<icoModel>(modelName);

        return rho()*model.nu();
    }
    else if (obr_.foundObject<cmpModel>(modelName))
    {
        const cmpModel& model =
            obr_.lookupObject<cmpModel>(modelName);

        return model.rho()*model.nu();
    }
    else if (obr_.foundObject<phaseIcoModel>(phaseModelName))
    {
        const phaseIcoModel& model =
            obr_.lookupObject<phaseIcoModel>(phaseModelName);

        return rho()*model.nu();
    }
    else if (obr_.foundObject<phaseCmpModel>(phaseModelName))
    {
        const phaseCmpModel& model =
            obr_.lookupObject<phaseCmpModel>(phaseModelName);

        return model.rho()*model.nu();
    }
    else if (obr_.foundObject<dictionary>("physicalProperties"))
    {
        // Legacy support for icoFoam

        const dictionary& physicalProperties =
             obr_.lookupObject<dictionary>("physicalProperties");

        const dimensionedScalar nu
        (
            "nu",
            dimKinematicViscosity,
            physicalProperties.lookup("nu")
        );

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forcesBase::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return volScalarField::New
        (
            "rho",
            mesh_,
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }
    else
    {
        return(obr_.lookupObject<volScalarField>(rhoName_));
    }
}


void Foam::functionObjects::forcesBase::calcForcesMoments(const vector& CofR)
{
    initialise();

    force_[0] = Zero;
    force_[1] = Zero;
    force_[2] = Zero;

    moment_[0] = Zero;
    moment_[1] = Zero;
    moment_[2] = Zero;

    if (directForceDensity_)
    {
        const volVectorField& fD = obr_.lookupObject<volVectorField>(fDName_);

        const surfaceVectorField::Boundary& Sfb =
            mesh_.Sf().boundaryField();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            const label patchi = iter.key();

            const vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - CofR
            );

            const scalarField sA(mag(Sfb[patchi]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            const vectorField fN
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );

            // Tangential force (total force minus normal fN)
            const vectorField fT(sA*fD.boundaryField()[patchi] - fN);

            //- Porous force
            const vectorField fP(Md.size(), Zero);

            applyBins(Md, fN, fT, fP, mesh_.C().boundaryField()[patchi]);
        }
    }
    else
    {
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        const surfaceVectorField::Boundary& Sfb =
            mesh_.Sf().boundaryField();

        const surfaceScalarField::Boundary& magSfb =
            mesh_.magSf().boundaryField();

        tmp<surfaceVectorField> tdevTau = devTau();
        const surfaceVectorField::Boundary& devTaub = tdevTau().boundaryField();

        // Scale pRef by density for incompressible simulations
        const scalar pRef = pRef_/rho(p);

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            const label patchi = iter.key();

            const vectorField Md
            (
                mesh_.C().boundaryField()[patchi] - CofR
            );

            const vectorField fN
            (
                alpha(patchi)
               *rho(p)
               *Sfb[patchi]
               *(p.boundaryField()[patchi] - pRef)
            );

            const vectorField fT(magSfb[patchi] * devTaub[patchi]);

            const vectorField fP(Md.size(), Zero);

            applyBins(Md, fN, fT, fP, mesh_.C().boundaryField()[patchi]);
        }
    }

    if (porosity_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel*>, models, iter)
        {
            const porosityModel& pm = *iter();

            const vectorField fPTot(pm.force(U, rho, mu));

            const cellZone& cZone = mesh_.cellZones()[pm.zoneName()];
            const vectorField d(mesh_.C(), cZone);
            const vectorField fP(fPTot, cZone);
            const vectorField Md(d - CofR);

            const vectorField fDummy(Md.size(), Zero);

            applyBins(Md, fDummy, fDummy, fP, d);
        }
    }

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineScatter(moment_);
}


void Foam::functionObjects::forcesBase::calcForcesMoments()
{
    forcesBase::calcForcesMoments(CofR());
}


Foam::vector Foam::functionObjects::forcesBase::forceEff() const
{
    return sum(force_[0]) + sum(force_[1]) + sum(force_[2]);
}


Foam::vector Foam::functionObjects::forcesBase::momentEff() const
{
    return sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2]);
}


bool Foam::functionObjects::forcesBase::execute()
{
    return true;
}


bool Foam::functionObjects::forcesBase::write()
{
    calcForcesMoments(CofR());

    if (Pstream::master())
    {
        logFiles::write();

        writeForces();

        writeBins();

        Log << endl;
    }

    return true;
}


// ************************************************************************* //
