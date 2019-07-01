/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "viscosityModelC.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscosityModelC, 0);
    defineRunTimeSelectionTable(viscosityModelC, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModelC::viscosityModelC
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& alpha1
)
:
    name_(name),
    viscosityProperties_(viscosityProperties),
    U_(U),
    phi_(phi),
    alpha1_(alpha1),
    mu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::viscosityModelC::strainRate(const volVectorField & U) const
{
    return sqrt(2.0)*mag(symm(fvc::grad(U)));
}

void Foam::viscosityModelC::correct()
{
    this->mu_ = this->calcMu(alpha1_, U_);
}

bool Foam::viscosityModelC::read(const dictionary& viscosityProperties)
{
    viscosityProperties_ = viscosityProperties;

    return true;
}


// ************************************************************************* //
