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

#include "WalburnSchneckC.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(WalburnSchneckC, 0);

    addToRunTimeSelectionTable
    (
        viscosityModelC,
        WalburnSchneckC,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::viscosityModels::WalburnSchneckC::calcMu(const volScalarField & alpha, const volVectorField & U) const
{
	volScalarField alpha1LimitedPercent(
			scalar(100.0)*min(max(alpha, scalar(0.25)), scalar(1))
	);

    return max(
		muMin_,
		min(
			muMax_,
			C1_*exp(0.0608*alpha1LimitedPercent)
			*exp(14.585 * TPMA_/ sqr(alpha1LimitedPercent))
			*pow(
				max
				(
					dimensionedScalar("one", dimTime, 1.0)*strainRate(U),
					dimensionedScalar("VSMALL", dimless, VSMALL)
				),
				(-scalar(0.00499) * alpha1LimitedPercent)
			)
		)
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::WalburnSchneckC::WalburnSchneckC
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& alpha1
)
:
    viscosityModelC(name, viscosityProperties, U, phi, alpha1),
    WalburnSchneckCCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    C1_(WalburnSchneckCCoeffs_.lookup("C1")),
    TPMA_(WalburnSchneckCCoeffs_.lookup("TPMA")),
    muMin_(WalburnSchneckCCoeffs_.lookup("muMin")),
    muMax_(WalburnSchneckCCoeffs_.lookup("muMax"))
{
	this->correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::WalburnSchneckC::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModelC::read(viscosityProperties);

    WalburnSchneckCCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    WalburnSchneckCCoeffs_.lookup("C1") >> C1_;
    WalburnSchneckCCoeffs_.lookup("TPMA") >> TPMA_;
    WalburnSchneckCCoeffs_.lookup("muMin") >> muMin_;
    WalburnSchneckCCoeffs_.lookup("muMax") >> muMax_;

    return true;
}


// ************************************************************************* //
