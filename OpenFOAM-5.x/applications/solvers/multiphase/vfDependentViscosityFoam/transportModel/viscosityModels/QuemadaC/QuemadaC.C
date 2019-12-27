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

// Helper functions
namespace Foam {
/*void quemadaCoefficients(
    const volScalarField & alpha, 
    volScalarField & k0, 
    volScalarField & kinf, 
    volScalarField & gammaC
)
{
    k0 = exp(3.874 - 10.41*alpha + 13.8*pow(alpha, 2) - 6.738*pow(alpha,3));
	kinf = exp(1.3435 - 2.803*alpha + 2.711*pow(alpha, 2) - 0.6479*pow(alpha, 3));
    gammaC = dimensionedScalar("oneOverTime", dimless/dimTime, 1.0) * exp(-6.1508 + 27.923*alpha - 25.6 * pow(alpha, 2) + 3.697 * pow(alpha, 3));
}*/

/*void quemadaCoefficientsAndDerivs(
    const volScalarField & alpha, 
    volScalarField & k0,
    volScalarField & dk0_dalpha,
    volScalarField & kinf,
    volScalarField & dkinf_dalpha,
    volScalarField & gammaC,
    volScalarField & dgammaC_dalpha
)
{
    quemadaCoefficients(alpha, k0, kinf, gammaC);
    dk0_dalpha = k0 * (-10.41 + 2*13.8*alpha - 3*6.738*pow(alpha, 2));
    dkinf_dalpha = kinf * (-2.803 + 2*2.711*alpha - 3*0.6479*pow(alpha, 2));
    dgammaC_dalpha = gammaC * (27.923 - 2*25.6 * alpha + 3*3.697 * pow(alpha, 2));
}*/

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::viscosityModels::QuemadaC::calcMu(const volScalarField & alpha, const volVectorField & U) const
{
	const volScalarField alphaLim( min(max(alpha, scalar(0)), scalar(1)) );
    volScalarField k0 = exp(3.874 - 10.41*alphaLim + 13.8*pow(alphaLim, 2) - 6.738*pow(alphaLim,3));
	volScalarField kinf = exp(1.3435 - 2.803*alphaLim + 2.711*pow(alphaLim, 2) - 0.6479*pow(alphaLim, 3));
    volScalarField gammaC = dimensionedScalar("oneOverTime", dimless/dimTime, 1.0) 
        * exp(-6.1508 + 27.923*alphaLim - 25.6 * pow(alphaLim, 2) + 3.697 * pow(alphaLim, 3));
	
    return min(muMax_, mup_ * 
        pow(1.0 - 0.5 * alphaLim * (k0 + kinf*sqrt(strainRate(U) / gammaC)) / 
            (1.0 + sqrt(strainRate(U) / gammaC)), -2.0)
    );
}

Foam::tmp<Foam::volScalarField> Foam::viscosityModels::QuemadaC::dMuDalpha() const
{
    const volScalarField alphaLim( min(max(alpha(), scalar(0)), scalar(1)) );

    const volScalarField k0 = exp(3.874 - 10.41*alphaLim + 13.8*pow(alphaLim, 2) - 6.738*pow(alphaLim,3));
	const volScalarField kinf = exp(1.3435 - 2.803*alphaLim + 2.711*pow(alphaLim, 2) - 0.6479*pow(alphaLim, 3));
    const volScalarField gammaC = dimensionedScalar("oneOverTime", dimless/dimTime, 1.0) * exp(-6.1508 + 27.923*alphaLim - 25.6 * pow(alphaLim, 2) + 3.697 * pow(alphaLim, 3));
    const volScalarField dk0_dalpha = k0 * (-10.41 + 2*13.8*alphaLim - 3*6.738*pow(alphaLim, 2));
    const volScalarField dkinf_dalpha = kinf * (-2.803 + 2*2.711*alphaLim - 3*0.6479*pow(alphaLim, 2));
    const volScalarField dgammaC_dalpha = gammaC * (27.923 - 2*25.6 * alphaLim + 3*3.697 * pow(alphaLim, 2));

    const volScalarField gamma = strainRate(U_);
    const volScalarField num = k0 + kinf*sqrt(gamma / gammaC);
    const volScalarField denum = 1 + sqrt(gamma/gammaC);

    return tmp<volScalarField>(
        new volScalarField(
            mup_ * pow(1. - 0.5*alphaLim*num/denum, -3) * (
                num/denum * (1 + alphaLim*sqrt(gamma) * dgammaC_dalpha / (2*denum*pow(gammaC, 1.5))) +
                alphaLim/denum * (dk0_dalpha + dkinf_dalpha*sqrt(gamma/gammaC) - 0.5*kinf*sqrt(gamma)*dgammaC_dalpha / pow(gammaC, 1.5))
            )
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::viscosityModels::QuemadaC::dMuDgamma() const
{
    const volScalarField alphaLim( min(max(alpha(), scalar(0)), scalar(1)) );
    const volScalarField k0 = exp(3.874 - 10.41*alphaLim + 13.8*pow(alphaLim, 2) - 6.738*pow(alphaLim,3));
	const volScalarField kinf = exp(1.3435 - 2.803*alphaLim + 2.711*pow(alphaLim, 2) - 0.6479*pow(alphaLim, 3));
    const volScalarField gammaC = dimensionedScalar("oneOverTime", dimless/dimTime, 1.0) * exp(-6.1508 + 27.923*alphaLim - 25.6 * pow(alphaLim, 2) + 3.697 * pow(alphaLim, 3));

    const volScalarField gamma = max(strainRate(U_), dimensionedScalar("small", dimless/dimTime, 1e-3));
    const volScalarField num = k0 + kinf*sqrt(gamma / gammaC);
    const volScalarField denum = 1 + sqrt(gamma/gammaC);

    return tmp<volScalarField>(
        new volScalarField(
            mup_ * pow(1. - 0.5*alphaLim*num/denum, -3) * 
                0.5 * alphaLim * (kinf - k0) / (denum * (gamma + sqrt(gamma*gammaC)))
        )
    );
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
