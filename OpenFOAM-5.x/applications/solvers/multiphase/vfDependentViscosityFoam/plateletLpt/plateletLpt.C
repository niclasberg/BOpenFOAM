/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "plateletLpt.H"
#include "vfDependentViscosityTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(plateletLpt, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        plateletLpt,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::plateletLpt::plateletLpt
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    g_
    (
        IOobject
        (
            "g",
            time_.constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector("g", dimAcceleration, Zero)
    ),
    transportModel_
    (
        mesh_.lookupObject<vfDependentViscosityTwoPhaseMixture>("transportProperties")
    ),
    rho_
    (
        mesh_.lookupObject<volScalarField>(dict.lookupOrDefault<word>("rho", "rho"))
    ),
    U_
    (
        mesh_.lookupObject<volVectorField>(dict.lookupOrDefault<word>("U", "U"))
    ),
    mu_(
		"mu",
		transportModel_.mu()
	),
	tau_(
		"tau",
		mu_ * dev(twoSymm(fvc::grad(U_)))
	),
    plateletCloudName_
    (
        dict.lookupOrDefault<word>("cloudName", "plateletCloud")
    ),
    plateletCloud_
    (
        plateletCloudName_,
        rho_,
        U_,
        mu_,
        g_,
		tau_
    )
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::plateletLpt::~plateletLpt()
{
    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::plateletLpt::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::plateletLpt::execute()
{
	// Update dynamic viscosity and shear stress field
	mu_ = rho_*transportModel_.nu();
	tau_ = mu_ * dev(twoSymm(fvc::grad(U_)));

	// Update cloud
    plateletCloud_.evolve();

    return true;
}


bool Foam::functionObjects::plateletLpt::write()
{
	plateletCloud_.write();
    return true;
}


// ************************************************************************* //
