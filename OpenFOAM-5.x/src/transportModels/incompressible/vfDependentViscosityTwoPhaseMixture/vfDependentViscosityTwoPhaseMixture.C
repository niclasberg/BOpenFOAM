/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "vfDependentViscosityTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vfDependentViscosityTwoPhaseMixture, 0);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::vfDependentViscosityTwoPhaseMixture::calcNu()
{
    nuModel_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vfDependentViscosityTwoPhaseMixture::vfDependentViscosityTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),

    nuModel_
    (
        viscosityModelC::New
        (
            "nu",
            *this,
            U,
            phi,
            alpha1_
        )
    ),

    rho1_("rho", dimDensity, subDict(phase1Name_)),
    rho2_("rho", dimDensity, subDict(phase2Name_)),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::vfDependentViscosityTwoPhaseMixture::mu() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            (limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_)*nuModel_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::vfDependentViscosityTwoPhaseMixture::muf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
            (alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)*fvc::interpolate(nuModel_->nu())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::vfDependentViscosityTwoPhaseMixture::nuf() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            fvc::interpolate(nuModel_->nu())
        )
    );
}


bool Foam::vfDependentViscosityTwoPhaseMixture::read()
{
    if (regIOobject::read())
    {
        if(nuModel_().read(*this))
        {
        	subDict(phase1Name_).lookup("rho") >> rho1_;
        	subDict(phase2Name_).lookup("rho") >> rho2_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
