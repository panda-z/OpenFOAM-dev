/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "HeatTransferPhaseSystem.H"
#include "fvmSup.H"
#include "rhoFluidMulticomponentThermo.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::HeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefs
(
    const phaseSystem::dmdtfTable& dmdtfs,
    phaseSystem::heatTransferTable& eqns
) const
{
    // Loop the pairs
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfIter)
    {
        const phaseInterface interface(*this, dmdtfIter.key());

        const volScalarField& dmdtf = *dmdtfIter();
        const volScalarField dmdtf21(posPart(dmdtf));
        const volScalarField dmdtf12(negPart(dmdtf));

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();
        const rhoFluidThermo& thermo1 = phase1.fluidThermo();
        const rhoFluidThermo& thermo2 = phase2.fluidThermo();
        const volScalarField& he1 = thermo1.he();
        const volScalarField& he2 = thermo2.he();
        const volScalarField hs1(thermo1.hs());
        const volScalarField hs2(thermo2.hs());
        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        // Transfer of sensible enthalpy within the phases
        *eqns[phase1.name()] +=
            dmdtf*hs1 + fvm::Sp(dmdtf12, he1) - dmdtf12*he1;
        *eqns[phase2.name()] -=
            dmdtf*hs2 + fvm::Sp(dmdtf21, he2) - dmdtf21*he2;

        // Transfer of sensible enthalpy between the phases
        *eqns[phase1.name()] += dmdtf21*(hs2 - hs1);
        *eqns[phase2.name()] -= dmdtf12*(hs1 - hs2);

        // Transfer of kinetic energy
        *eqns[phase1.name()] += dmdtf21*K2 + dmdtf12*K1;
        *eqns[phase2.name()] -= dmdtf12*K1 + dmdtf21*K2;
    }
}


template<class BasePhaseSystem>
void Foam::HeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHefs
(
    const phaseSystem::dmidtfTable& dmidtfs,
    phaseSystem::heatTransferTable& eqns
) const
{
    static const dimensionedScalar one(dimless, 1);

    // Loop the pairs
    forAllConstIter(phaseSystem::dmidtfTable, dmidtfs, dmidtfIter)
    {
        const phaseInterface interface(*this, dmidtfIter.key());

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();
        const rhoFluidThermo& thermo1 = phase1.fluidThermo();
        const rhoFluidThermo& thermo2 = phase2.fluidThermo();
        const rhoFluidMulticomponentThermo* mcThermoPtr1 =
            isA<rhoFluidMulticomponentThermo>(thermo1)
          ? &refCast<const rhoFluidMulticomponentThermo>(thermo1)
          : static_cast<const rhoFluidMulticomponentThermo*>(nullptr);
        const rhoFluidMulticomponentThermo* mcThermoPtr2 =
            isA<rhoFluidMulticomponentThermo>(thermo2)
          ? &refCast<const rhoFluidMulticomponentThermo>(thermo2)
          : static_cast<const rhoFluidMulticomponentThermo*>(nullptr);
        const volScalarField& he1 = thermo1.he();
        const volScalarField& he2 = thermo2.he();
        const volScalarField hs1(thermo1.hs());
        const volScalarField hs2(thermo2.hs());
        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        // Loop the species
        forAllConstIter(HashPtrTable<volScalarField>, *dmidtfIter(), dmidtfJter)
        {
            const word& specie = dmidtfJter.key();

            // Mass transfer rates
            const volScalarField& dmidtf = *dmidtfJter();
            const volScalarField dmidtf21(posPart(dmidtf));
            const volScalarField dmidtf12(negPart(dmidtf));

            // Specie indices
            const label speciei1 =
                mcThermoPtr1 ? mcThermoPtr1->species()[specie] : -1;
            const label speciei2 =
                mcThermoPtr2 ? mcThermoPtr2->species()[specie] : -1;

            // Enthalpies
            const volScalarField hsi1
            (
                mcThermoPtr1
              ? mcThermoPtr1->hsi(speciei1, thermo1.p(), thermo1.T())
              : tmp<volScalarField>(hs1)
            );
            const volScalarField hsi2
            (
                mcThermoPtr2
              ? mcThermoPtr2->hsi(speciei2, thermo2.p(), thermo2.T())
              : tmp<volScalarField>(hs2)
            );

            // Limited mass fractions
            tmp<volScalarField> tYi1, tYi2;
            if (residualY_ > 0)
            {
                tYi1 =
                    mcThermoPtr1
                  ? max(mcThermoPtr1->Y(speciei1), residualY_)
                  : volScalarField::New("Yi1", this->mesh(), one);
                tYi2 =
                    mcThermoPtr2
                  ? max(mcThermoPtr2->Y(speciei2), residualY_)
                  : volScalarField::New("Yi2", this->mesh(), one);
            }

            // Transfer of sensible enthalpy within the phases
            *eqns[phase1.name()] += dmidtf*hsi1;
            *eqns[phase2.name()] -= dmidtf*hsi2;
            if (residualY_ > 0)
            {
                *eqns[phase1.name()] +=
                    fvm::Sp(dmidtf12/tYi1(), he1) - dmidtf12/tYi1()*he1;
                *eqns[phase2.name()] -=
                    fvm::Sp(dmidtf21/tYi2(), he2) - dmidtf21/tYi2()*he2;
            }

            // Transfer of sensible enthalpy between the phases
            *eqns[phase1.name()] += dmidtf21*(hsi2 - hsi1);
            *eqns[phase2.name()] -= dmidtf12*(hsi1 - hsi2);

            // Transfer of kinetic energy
            *eqns[phase1.name()] += dmidtf21*K2 + dmidtf12*K1;
            *eqns[phase2.name()] -= dmidtf12*K1 + dmidtf21*K2;
        }
    }
}


template<class BasePhaseSystem>
void Foam::HeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefsWithoutL
(
    const phaseSystem::dmdtfTable& dmdtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const latentHeatScheme scheme,
    phaseSystem::heatTransferTable& eqns
) const
{
    // Loop the pairs
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfIter)
    {
        const phaseInterface interface(*this, dmdtfIter.key());

        const volScalarField& dmdtf = *dmdtfIter();
        const volScalarField dmdtf21(posPart(dmdtf));
        const volScalarField dmdtf12(negPart(dmdtf));

        const volScalarField& Tf = *Tfs[dmdtfIter.key()];

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();
        const rhoFluidThermo& thermo1 = phase1.fluidThermo();
        const rhoFluidThermo& thermo2 = phase2.fluidThermo();
        const volScalarField& he1 = thermo1.he();
        const volScalarField& he2 = thermo2.he();
        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        // Interface enthalpies
        const volScalarField hsf1(thermo1.hs(thermo1.p(), Tf));
        const volScalarField hsf2(thermo2.hs(thermo1.p(), Tf));

        // Transfer of energy from the interface into the bulk
        switch (scheme)
        {
            case latentHeatScheme::symmetric:
            {
                *eqns[phase1.name()] += dmdtf*hsf1;
                *eqns[phase2.name()] -= dmdtf*hsf2;

                break;
            }
            case latentHeatScheme::upwind:
            {
                // Bulk enthalpies
                const volScalarField hs1(thermo1.hs());
                const volScalarField hs2(thermo2.hs());

                *eqns[phase1.name()] += dmdtf21*hsf1 + dmdtf12*hs1;
                *eqns[phase2.name()] -= dmdtf12*hsf2 + dmdtf21*hs2;

                break;
            }
        }
        *eqns[phase1.name()] += fvm::Sp(dmdtf12, he1) - dmdtf12*he1;
        *eqns[phase2.name()] -= fvm::Sp(dmdtf21, he2) - dmdtf21*he2;

        // Transfer of kinetic energy
        *eqns[phase1.name()] += dmdtf21*K2 + dmdtf12*K1;
        *eqns[phase2.name()] -= dmdtf12*K1 + dmdtf21*K2;
    }
}


template<class BasePhaseSystem>
void Foam::HeatTransferPhaseSystem<BasePhaseSystem>::addDmdtL
(
    const phaseSystem::dmdtfTable& dmdtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const scalar weight,
    const latentHeatScheme scheme,
    phaseSystem::heatTransferTable& eqns
) const
{
    // Loop the pairs
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfIter)
    {
        const phaseInterface interface(*this, dmdtfIter.key());

        const volScalarField& dmdtf = *dmdtfIter();
        const volScalarField dmdtf21(posPart(dmdtf));
        const volScalarField dmdtf12(negPart(dmdtf));

        const volScalarField& Tf = *Tfs[dmdtfIter.key()];

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();

        // Latent heat contribution
        const volScalarField L(this->L(interface, dmdtf, Tf, scheme));
        *eqns[phase1.name()] += ((1 - weight)*dmdtf12 + weight*dmdtf21)*L;
        *eqns[phase2.name()] += ((1 - weight)*dmdtf21 + weight*dmdtf12)*L;
    }
}


template<class BasePhaseSystem>
void Foam::HeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefs
(
    const phaseSystem::dmdtfTable& dmdtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const scalar weight,
    const latentHeatScheme scheme,
    phaseSystem::heatTransferTable& eqns
) const
{
    addDmdtHefsWithoutL(dmdtfs, Tfs, scheme, eqns);
    addDmdtL(dmdtfs, Tfs, weight, scheme, eqns);
}


template<class BasePhaseSystem>
void Foam::HeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHefsWithoutL
(
    const phaseSystem::dmidtfTable& dmidtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const latentHeatScheme scheme,
    phaseSystem::heatTransferTable& eqns
) const
{
    static const dimensionedScalar one(dimless, 1);

    // Loop the pairs
    forAllConstIter(phaseSystem::dmidtfTable, dmidtfs, dmidtfIter)
    {
        const phaseInterface interface(*this, dmidtfIter.key());

        const volScalarField& Tf = *Tfs[dmidtfIter.key()];

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();
        const rhoFluidThermo& thermo1 = phase1.fluidThermo();
        const rhoFluidThermo& thermo2 = phase2.fluidThermo();
        const rhoFluidMulticomponentThermo* mcThermoPtr1 =
            isA<rhoFluidMulticomponentThermo>(thermo1)
          ? &refCast<const rhoFluidMulticomponentThermo>(thermo1)
          : static_cast<const rhoFluidMulticomponentThermo*>(nullptr);
        const rhoFluidMulticomponentThermo* mcThermoPtr2 =
            isA<rhoFluidMulticomponentThermo>(thermo2)
          ? &refCast<const rhoFluidMulticomponentThermo>(thermo2)
          : static_cast<const rhoFluidMulticomponentThermo*>(nullptr);
        const volScalarField& he1 = thermo1.he();
        const volScalarField& he2 = thermo2.he();
        const volScalarField K1(phase1.K());
        const volScalarField K2(phase2.K());

        // Interface enthalpies
        const volScalarField hsf1(thermo1.hs(thermo1.p(), Tf));
        const volScalarField hsf2(thermo2.hs(thermo2.p(), Tf));

        // Loop the species
        forAllConstIter(HashPtrTable<volScalarField>, *dmidtfIter(), dmidtfJter)
        {
            const word& specie = dmidtfJter.key();

            // Mass transfer rates
            const volScalarField& dmidtf = *dmidtfJter();
            const volScalarField dmidtf21(posPart(dmidtf));
            const volScalarField dmidtf12(negPart(dmidtf));

            // Specie indices
            const label speciei1 =
                mcThermoPtr1 ? mcThermoPtr1->species()[specie] : -1;
            const label speciei2 =
                mcThermoPtr2 ? mcThermoPtr2->species()[specie] : -1;

            // Interface enthalpies
            const volScalarField hsfi1
            (
                mcThermoPtr1
              ? mcThermoPtr1->hsi(speciei1, thermo1.p(), Tf)
              : tmp<volScalarField>(hsf1)
            );
            const volScalarField hsfi2
            (
                mcThermoPtr2
              ? mcThermoPtr2->hsi(speciei2, thermo2.p(), Tf)
              : tmp<volScalarField>(hsf2)
            );

            // Limited mass fractions
            tmp<volScalarField> tYi1, tYi2;
            if (this->residualY_ > 0)
            {
                tYi1 =
                    mcThermoPtr1
                  ? max(mcThermoPtr1->Y(speciei1), this->residualY_)
                  : volScalarField::New("Yi1", this->mesh(), one);
                tYi2 =
                    mcThermoPtr2
                  ? max(mcThermoPtr2->Y(speciei2), this->residualY_)
                  : volScalarField::New("Yi2", this->mesh(), one);
            }

            // Transfer of energy from the interface into the bulk
            switch (scheme)
            {
                case latentHeatScheme::symmetric:
                {
                    *eqns[phase1.name()] += dmidtf*hsfi1;
                    *eqns[phase2.name()] -= dmidtf*hsfi2;

                    break;
                }
                case latentHeatScheme::upwind:
                {
                    // Bulk enthalpies
                    const volScalarField hsi1
                    (
                        mcThermoPtr1
                      ? mcThermoPtr1->hsi(speciei1, thermo1.p(), thermo1.T())
                      : thermo1.hs()
                    );
                    const volScalarField hsi2
                    (
                        mcThermoPtr2
                      ? mcThermoPtr2->hsi(speciei2, thermo2.p(), thermo2.T())
                      : thermo2.hs()
                    );

                    *eqns[phase1.name()] += dmidtf21*hsfi1 + dmidtf12*hsi1;
                    *eqns[phase2.name()] -= dmidtf12*hsfi2 + dmidtf21*hsi2;

                    break;
                }
            }
            if (this->residualY_ > 0)
            {
                *eqns[phase1.name()] +=
                    fvm::Sp(dmidtf12/tYi1(), he1) - dmidtf12/tYi1()*he1;
            }
            if (this->residualY_ > 0)
            {
                *eqns[phase2.name()] -=
                    fvm::Sp(dmidtf21/tYi2(), he2) - dmidtf21/tYi2()*he2;
            }


            // Transfer of kinetic energy
            *eqns[phase1.name()] += dmidtf21*K2 + dmidtf12*K1;
            *eqns[phase2.name()] -= dmidtf12*K1 + dmidtf21*K2;
        }
    }
}


template<class BasePhaseSystem>
void Foam::HeatTransferPhaseSystem<BasePhaseSystem>::addDmidtL
(
    const phaseSystem::dmidtfTable& dmidtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const scalar weight,
    const latentHeatScheme scheme,
    phaseSystem::heatTransferTable& eqns
) const
{
    // Loop the pairs
    forAllConstIter(phaseSystem::dmidtfTable, dmidtfs, dmidtfIter)
    {
        const phaseInterface interface(*this, dmidtfIter.key());

        const volScalarField& Tf = *Tfs[dmidtfIter.key()];

        const phaseModel& phase1 = interface.phase1();
        const phaseModel& phase2 = interface.phase2();

        // Loop the species
        forAllConstIter(HashPtrTable<volScalarField>, *dmidtfIter(), dmidtfJter)
        {
            const word& specie = dmidtfJter.key();

            // Mass transfer rates
            const volScalarField& dmidtf = *dmidtfJter();
            const volScalarField dmidtf21(posPart(dmidtf));
            const volScalarField dmidtf12(negPart(dmidtf));

            // Latent heat contribution
            const volScalarField Li
            (
                this->Li(interface, specie, dmidtf, Tf, scheme)
            );
            *eqns[phase1.name()] +=
                ((1 - weight)*dmidtf12 + weight*dmidtf21)*Li;
            *eqns[phase2.name()] +=
                ((1 - weight)*dmidtf21 + weight*dmidtf12)*Li;
        }
    }
}


template<class BasePhaseSystem>
void Foam::HeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHefs
(
    const phaseSystem::dmidtfTable& dmidtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const scalar weight,
    const latentHeatScheme scheme,
    phaseSystem::heatTransferTable& eqns
) const
{
    addDmidtHefsWithoutL(dmidtfs, Tfs, scheme, eqns);
    addDmidtL(dmidtfs, Tfs, weight, scheme, eqns);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::HeatTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    heatTransferPhaseSystem(),
    BasePhaseSystem(mesh),
    residualY_(this->template lookupOrDefault<scalar>("residualY", -1))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::~HeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::L
(
    const phaseInterface& interface,
    const volScalarField& dmdtf,
    const volScalarField& Tf,
    const latentHeatScheme scheme
) const
{
    const rhoFluidThermo& thermo1 = interface.phase1().fluidThermo();
    const rhoFluidThermo& thermo2 = interface.phase2().fluidThermo();

    // Interface enthalpies
    const volScalarField haf1(thermo1.ha(thermo1.p(), Tf));
    const volScalarField haf2(thermo2.ha(thermo2.p(), Tf));

    switch (scheme)
    {
        case latentHeatScheme::symmetric:
        {
            return haf2 - haf1;
        }
        case latentHeatScheme::upwind:
        {
            // Bulk enthalpies
            const volScalarField ha1(thermo1.ha());
            const volScalarField ha2(thermo2.ha());

            return
                neg0(dmdtf)*haf2 + pos(dmdtf)*ha2
              - pos0(dmdtf)*haf1 - neg(dmdtf)*ha1;
        }
    }

    return tmp<volScalarField>(nullptr);
}


template<class BasePhaseSystem>
Foam::tmp<Foam::scalarField>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::L
(
    const phaseInterface& interface,
    const scalarField& dmdtf,
    const scalarField& Tf,
    const labelUList& cells,
    const latentHeatScheme scheme
) const
{
    const rhoFluidThermo& thermo1 = interface.phase1().fluidThermo();
    const rhoFluidThermo& thermo2 = interface.phase2().fluidThermo();

    // Interface enthalpies
    const scalarField haf1(thermo1.ha(Tf, cells));
    const scalarField haf2(thermo2.ha(Tf, cells));

    switch (scheme)
    {
        case latentHeatScheme::symmetric:
        {
            return haf2 - haf1;
        }
        case latentHeatScheme::upwind:
        {
            const scalarField T1(UIndirectList<scalar>(thermo1.T(), cells));
            const scalarField T2(UIndirectList<scalar>(thermo2.T(), cells));

            // Bulk enthalpies
            const scalarField ha1(thermo1.ha(T1, cells));
            const scalarField ha2(thermo2.ha(T2, cells));

            return
                neg0(dmdtf)*haf2 + pos(dmdtf)*ha2
              - pos0(dmdtf)*haf1 - neg(dmdtf)*ha1;
        }
    }

    return tmp<scalarField>(nullptr);
}


template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::Li
(
    const phaseInterface& interface,
    const word& specie,
    const volScalarField& dmdtf,
    const volScalarField& Tf,
    const latentHeatScheme scheme
) const
{
    const rhoFluidThermo& thermo1 = interface.phase1().fluidThermo();
    const rhoFluidThermo& thermo2 = interface.phase2().fluidThermo();
    const rhoFluidMulticomponentThermo* mcThermoPtr1 =
        isA<rhoFluidMulticomponentThermo>(thermo1)
      ? &refCast<const rhoFluidMulticomponentThermo>(thermo1)
      : static_cast<const rhoFluidMulticomponentThermo*>(nullptr);
    const rhoFluidMulticomponentThermo* mcThermoPtr2 =
        isA<rhoFluidMulticomponentThermo>(thermo2)
      ? &refCast<const rhoFluidMulticomponentThermo>(thermo2)
      : static_cast<const rhoFluidMulticomponentThermo*>(nullptr);
    const label speciei1 =
        mcThermoPtr1 ? mcThermoPtr1->species()[specie] : -1;
    const label speciei2 =
        mcThermoPtr2 ? mcThermoPtr2->species()[specie] : -1;

    // Interface enthalpies
    const volScalarField hafi1
    (
        mcThermoPtr1
      ? mcThermoPtr1->hai(speciei1, thermo1.p(), Tf)
      : thermo1.ha(thermo1.p(), Tf)
    );
    const volScalarField hafi2
    (
        mcThermoPtr2
      ? mcThermoPtr2->hai(speciei2, thermo2.p(), Tf)
      : thermo2.ha(thermo1.p(), Tf)
    );

    switch (scheme)
    {
        case latentHeatScheme::symmetric:
        {
            return hafi2 - hafi1;
        }
        case latentHeatScheme::upwind:
        {
            // Bulk enthalpies
            const volScalarField hai1
            (
                mcThermoPtr1
              ? mcThermoPtr1->hai(speciei1, thermo1.p(), thermo1.T())
              : thermo1.ha()
            );
            const volScalarField hai2
            (
                mcThermoPtr2
              ? mcThermoPtr2->hai(speciei2, thermo2.p(), thermo2.T())
              : thermo2.ha()
            );

            return
                neg0(dmdtf)*hafi2 + pos(dmdtf)*hai2
              - pos0(dmdtf)*hafi1 - neg(dmdtf)*hai1;
        }
    }

    return tmp<volScalarField>(nullptr);
}


template<class BasePhaseSystem>
Foam::tmp<Foam::scalarField>
Foam::HeatTransferPhaseSystem<BasePhaseSystem>::Li
(
    const phaseInterface& interface,
    const word& specie,
    const scalarField& dmdtf,
    const scalarField& Tf,
    const labelUList& cells,
    const latentHeatScheme scheme
) const
{
    const rhoFluidThermo& thermo1 = interface.phase1().fluidThermo();
    const rhoFluidThermo& thermo2 = interface.phase2().fluidThermo();
    const rhoFluidMulticomponentThermo* mcThermoPtr1 =
        isA<rhoFluidMulticomponentThermo>(thermo1)
      ? &refCast<const rhoFluidMulticomponentThermo>(thermo1)
      : static_cast<const rhoFluidMulticomponentThermo*>(nullptr);
    const rhoFluidMulticomponentThermo* mcThermoPtr2 =
        isA<rhoFluidMulticomponentThermo>(thermo2)
      ? &refCast<const rhoFluidMulticomponentThermo>(thermo2)
      : static_cast<const rhoFluidMulticomponentThermo*>(nullptr);
    const label speciei1 =
        mcThermoPtr1 ? mcThermoPtr1->species()[specie] : -1;
    const label speciei2 =
        mcThermoPtr2 ? mcThermoPtr2->species()[specie] : -1;

    const scalarField p1(UIndirectList<scalar>(thermo1.p(), cells));
    const scalarField p2(UIndirectList<scalar>(thermo2.p(), cells));

    // Interface enthalpies
    const scalarField hafi1
    (
        mcThermoPtr1
      ? mcThermoPtr1->hai(speciei1, p1, Tf)
      : thermo1.ha(Tf, cells)
    );
    const scalarField hafi2
    (
        mcThermoPtr2
      ? mcThermoPtr2->hai(speciei2, p2, Tf)
      : thermo2.ha(Tf, cells)
    );

    switch (scheme)
    {
        case latentHeatScheme::symmetric:
        {
            return hafi2 - hafi1;
        }
        case latentHeatScheme::upwind:
        {
            const scalarField T1(UIndirectList<scalar>(thermo1.T(), cells));
            const scalarField T2(UIndirectList<scalar>(thermo2.T(), cells));

            // Bulk enthalpies
            const scalarField hai1
            (
                mcThermoPtr1
              ? mcThermoPtr1->hai(speciei1, p1, T1)
              : thermo1.ha(T1, cells)
            );
            const scalarField hai2
            (
                mcThermoPtr2
              ? mcThermoPtr2->hai(speciei2, p2, T2)
              : thermo2.ha(T2, cells)
            );

            return
                neg0(dmdtf)*hafi2 + pos(dmdtf)*hai2
              - pos0(dmdtf)*hafi1 - neg(dmdtf)*hai1;
        }
    }

    return tmp<scalarField>(nullptr);
}


template<class BasePhaseSystem>
bool Foam::HeatTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
