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

#include "CassonC.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(CassonC, 0);

    addToRunTimeSelectionTable
    (
        viscosityModelC,
        CassonC,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::CassonC::calcNu() const
{
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
    tmp<volScalarField> sr(strainRate());
    return (min(nuMax_, sqr( (B_ / A_) 
	    * (scalar(1.0) / pow(scalar(1.0) 
		   - min(max(VF(), scalar(0)), scalar(1)),
		   ( A_ / scalar(2.0))) - scalar(1.0) )
	   * scalar(1.0) /pow(max(sr(), dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)),
		( scalar(1.0) / scalar(2.0)))

	    + 
	    sqrt( mup_ / 
	    pow(scalar(1.0) 
	 	- min(max(VF(), scalar(0)), scalar(1)),A_ ))
	 	            ) 
	         /( rho1_*min(max(VF(), scalar(0)), scalar(1))
		 + rho2_*(scalar(1.0)-min(max(VF(), scalar(0)), scalar(1)))
		  )
	        )
	    ) ;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::CassonC::CassonC
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& alpha1
)
:

    viscosityModelC(name, viscosityProperties, U, phi, alpha1),
    CassonCCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    mup_(CassonCCoeffs_.lookup("mup")),
    nuMax_(CassonCCoeffs_.lookup("nuMax")),
    A_(CassonCCoeffs_.lookup("A")),
    B_(CassonCCoeffs_.lookup("B")),
    rho1_(CassonCCoeffs_.lookup("rho1")),
    rho2_(CassonCCoeffs_.lookup("rho2")),
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

bool Foam::viscosityModels::CassonC::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModelC::read(viscosityProperties);

    CassonCCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    CassonCCoeffs_.lookup("mup") >> mup_;
    CassonCCoeffs_.lookup("nuMax") >> nuMax_;
    CassonCCoeffs_.lookup("A") >> A_;
    CassonCCoeffs_.lookup("B") >> B_;
    CassonCCoeffs_.lookup("rho1") >> rho1_;
    CassonCCoeffs_.lookup("rho2") >> rho2_;

    return true;
}


// ************************************************************************* //
