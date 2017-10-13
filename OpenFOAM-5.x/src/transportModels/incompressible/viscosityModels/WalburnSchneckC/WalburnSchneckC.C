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

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::WalburnSchneckC::calcNu() const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            C1_
	    *exp
	    (scalar(0.0608)
	    *scalar(100.0)*min(max(VF(), scalar(0.25)), scalar(1))
	    )
	    *exp
	    (scalar(14.585)*TPMA_ 
	    / sqr(scalar(100.0)*min(max(VF(), scalar(0.25)), scalar(1)))
	    )
	    /( rho1_*min(max(VF(), scalar(0.25)), scalar(1))
	    + rho2_*(scalar(1.0)-min(max(VF(), scalar(0.25)), scalar(1)))
	     )
	    *pow
            (
                max
                (
                    dimensionedScalar("one", dimTime, 1.0)*strainRate(),
                    dimensionedScalar("VSMALL", dimless, VSMALL)
                ),
                (scalar(1.0) 
		- scalar(0.00499) 
		*scalar(100.0)*min(max(VF(), scalar(0.25)), scalar(1))
		- scalar(1.0)
		)
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
    nuMin_(WalburnSchneckCCoeffs_.lookup("nuMin")),
    nuMax_(WalburnSchneckCCoeffs_.lookup("nuMax")),
    rho1_(WalburnSchneckCCoeffs_.lookup("rho1")),
    rho2_(WalburnSchneckCCoeffs_.lookup("rho2")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        calcNu()
    )
{}


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
    WalburnSchneckCCoeffs_.lookup("nuMin") >> nuMin_;
    WalburnSchneckCCoeffs_.lookup("nuMax") >> nuMax_;
    WalburnSchneckCCoeffs_.lookup("rho1") >> rho1_;
    WalburnSchneckCCoeffs_.lookup("rho2") >> rho2_;

    return true;
}


// ************************************************************************* //
