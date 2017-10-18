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

#include "QuemadaC.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(QuemadaC, 0);

    addToRunTimeSelectionTable
    (
        viscosityModelC,
        QuemadaC,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::QuemadaC::QuemadaC
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& alpha1
)
:

    viscosityModelC(name, viscosityProperties, U, phi, alpha1),
    QuemadaCCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    mup_(QuemadaCCoeffs_.lookup("mup")),
    muMax_(QuemadaCCoeffs_.lookup("muMax"))
{
	correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::viscosityModels::QuemadaC::correct()
{
	const volScalarField limitedVF( min(max(VF(), scalar(0)), scalar(1)) );
	const volScalarField k0 = exp(3.874 - 10.41*limitedVF + 13.8*pow(limitedVF, 2) - 6.738*pow(limitedVF,3));
	const volScalarField kinf = exp(1.3435 - 2.803*limitedVF + 2.711*pow(limitedVF, 2) - 0.6479*pow(limitedVF, 3));
	const volScalarField sqrtGammaOverGammaC =
			sqrt(max(dimensionedScalar("one", dimTime, 1.0)*strainRate(), dimensionedScalar("VSMALL", dimless, VSMALL))
			/ exp(-6.1508 + 27.923*limitedVF - 25.6 * pow(limitedVF, 2) + 3.697 * pow(limitedVF, 3)));

    this->mu_ = min(muMax_, mup_ * pow(1.0 - 0.5 * limitedVF * (k0 + kinf*sqrtGammaOverGammaC) / (1.0 + sqrtGammaOverGammaC), -2.0));
}

bool Foam::viscosityModels::QuemadaC::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModelC::read(viscosityProperties);

    QuemadaCCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    QuemadaCCoeffs_.lookup("mup") >> mup_;
    QuemadaCCoeffs_.lookup("muMax") >> muMax_;

    return true;
}


// ************************************************************************* //
