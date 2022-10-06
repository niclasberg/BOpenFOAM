/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "powerLawHemolysis.H"
#include "surfaceFields.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{
    defineTypeNameAndDebug(powerLawHemolysis, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        powerLawHemolysis,
        dictionary
    );

}
}


// ******** Constructors ******* //


Foam::functionObjects::powerLawHemolysis::powerLawHemolysis
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "H")),
    nCorr_(0),
    fvOptions_(mesh_),
    s_
    (
        IOobject
        (
            fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    read(dict);
}

// ******** Destructor ******** //

Foam::functionObjects::powerLawHemolysis::~powerLawHemolysis()
{}

// ******** Member Functions ******** //

bool Foam::functionObjects::powerLawHemolysis::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    schemesField_ = dict.lookupOrDefault<word>("schemesField", fieldName_);

    
    // GW set by default - Gonna change this later
    C_ = dict.lookupOrDefault("C", 3.63E-7);
    alpha_ = dict.lookupOrDefault("alpha", 2.416);
    beta_ = dict.lookupOrDefault("beta", 0.785);
    rhoValue_ = dict.lookupOrDefault("rho", 1000);

    dict.readIfPresent("nCorr", nCorr_);

    if (dict.found("fvOptions"))
    {
        fvOptions_.reset(dict.subDict("fvOptions"));
    }

    return true;
}


bool Foam::functionObjects::powerLawHemolysis::execute()
{
    Info << type() << " write:" << endl;

    const surfaceScalarField& phi =
    mesh_.lookupObject<surfaceScalarField>(phiName_);

    // Calculate the diffusivity

    word divScheme("div(phi," + schemesField_ + ")");

    //Under relaxation coefficient
    scalar relaxCoeff = 0.0;
    if (mesh_.relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.equationRelaxationFactor(schemesField_);
    }

    if (phi.dimensions() == dimMass/dimTime) // Compressible
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        typedef compressible::turbulenceModel cmpModel;
        const cmpModel& model = mesh_.lookupObject<cmpModel>
        (
            turbulenceModel::propertiesName
        );
        volSymmTensorField tau = model.nu()*symm(fvc::grad(U));
        volScalarField tauScalar = sqrt(-0.5*(sqr(tr(tau))-tr(tau & tau)));
        for (label i = 0; i <= nCorr_; i++)
        {
            // I'm using model A_E from Yu et al. (2017)
            fvScalarMatrix sEqn
            (
                  fvm::ddt(rho, s_) 
                + fvm::div(phi, s_, divScheme)
                ==
                  pow(tauScalar,-alpha_/beta_)*pow(C_,1/beta_)
                  + fvm::Sp(0-pow(tauScalar,-alpha_/beta_)*pow(C_,1/beta_),s_)
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            sEqn.solve(mesh_.solverDict(schemesField_));

        }

    }
    else if (phi.dimensions() == dimVolume/dimTime) // Incompressible
    {
        typedef incompressible::turbulenceModel icoModel;
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        const icoModel& model = mesh_.lookupObject<icoModel>
        (
            turbulenceModel::propertiesName
        );
        volSymmTensorField tau = rhoValue_*model.nu()*symm(fvc::grad(U));
        volScalarField tauScalar = dimensionedScalar("one", dimDensity/dimPressure, 1) * 2 * sqrt(-0.5*(sqr(tr(tau))-tr(tau & tau)));
        Info << max(tauScalar) << endl;
        Info << min(tauScalar) << endl;
	//Info << max(pow(dimensionedScalar("one", dimensionSet(0, 0, -beta_, 0, 0,0, 0),C_),1/beta_)) << endl;
	Info << max(pow(tauScalar,-alpha_/beta_)*pow(dimensionedScalar("one", dimensionSet(0, 0, -beta_, 0, 0,0, 0),C_),1/beta_)) << endl;
	Info << min(pow(tauScalar,-alpha_/beta_)*pow(dimensionedScalar("one", dimensionSet(0, 0, -beta_, 0, 0,0, 0),C_),1/beta_)) << endl;
        //Info << tauScalar.dimensions() << endl;
        for (label i = 0; i <= nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(s_)
              + fvm::div(phi, s_, divScheme)
              ==
                pow(tauScalar,alpha_/beta_)*pow(dimensionedScalar("one", dimensionSet(0, 0, -beta_, 0, 0,0, 0),C_),1/beta_) // dimensionedScalar("one", dimless/dimTime,1)
              + fvm::Sp(-pow(tauScalar,alpha_/beta_)*pow(dimensionedScalar("one", dimensionSet(0, 0, -beta_, 0, 0,0, 0),C_),1/beta_),s_)// & dimensionedScalar("one", dimless/dimTime,1)  
            );

            sEqn.relax(relaxCoeff);

            fvOptions_.constrain(sEqn);

            sEqn.solve(mesh_.solverDict(schemesField_));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or " 
            << dimVolume/dimTime << exit(FatalError);
    }

    Info << endl;

    return true;

}
    
bool Foam::functionObjects::powerLawHemolysis::write()
{
    return true;
}

// ************************************************* //
