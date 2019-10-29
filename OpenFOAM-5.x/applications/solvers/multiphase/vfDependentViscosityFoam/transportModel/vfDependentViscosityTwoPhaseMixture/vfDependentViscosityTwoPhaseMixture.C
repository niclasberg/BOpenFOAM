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
    muModel_->correct();

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

    muModel_
    (
        viscosityModelC::New
        (
            "mu1",
            *this,
            U,
            phi,
            alpha1_
        )
    ),

	diffusionModel_
	(
		diffusionModel::New
		(
			*this,
			U, 
			phi,
            muModel_
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
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::surfaceScalarField>
Foam::vfDependentViscosityTwoPhaseMixture::muf() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
            fvc::interpolate(mu())
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
            fvc::interpolate(nu())
        )
    );
}


bool Foam::vfDependentViscosityTwoPhaseMixture::read()
{
    if (regIOobject::read())
    {
        if(muModel_().read(*this))
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
