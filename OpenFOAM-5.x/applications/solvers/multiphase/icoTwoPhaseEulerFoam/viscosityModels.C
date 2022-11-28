#include "viscosityModels.H"
#include "fvc.H"
#include "fvm.H"
#include "twoPhaseModel.H"
#include "addToRunTimeSelectionTable.H"

// Create entries for the runtime selection table
namespace Foam
{
    defineTypeNameAndDebug(ViscosityModel, 0);
    defineRunTimeSelectionTable(ViscosityModel, viscModel);
    
    defineTypeNameAndDebug(Rbc, 0);
    addToRunTimeSelectionTable(ViscosityModel, Rbc, viscModel);

    defineTypeNameAndDebug(Plasma, 0);
    addToRunTimeSelectionTable(ViscosityModel, Plasma, viscModel);
}

// ViscosityModel base class methods
Foam::ViscosityModel::ViscosityModel(const dictionary & dict, PhaseModel & phaseModel)
:
    phaseModel_(phaseModel)
{

}

Foam::autoPtr<Foam::ViscosityModel> Foam::ViscosityModel::New(const dictionary& dict, PhaseModel & phaseModel)
{
    const word modelType(dict.lookup("viscosityModel"));

    viscModelConstructorTable::iterator cstrIter = viscModelConstructorTablePtr_->find(modelType);

    if (cstrIter == viscModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "ViscosityModel::New(...)"
        )   << "Unknown viscosityModel type "
            << modelType << nl << nl
            << "Valid viscosityModels are : " << endl
            << viscModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ViscosityModel>(cstrIter()(dict, phaseModel));
}

Foam::tmp<Foam::volScalarField> Foam::ViscosityModel::shearRate() const
{
    return sqrt(2.0)*mag(symm(fvc::grad(phaseModel_.U())));
}   

// Plasma model
Foam::Plasma::Plasma(const Foam::dictionary & dict, Foam::PhaseModel & phaseModel)
:
    ViscosityModel(dict, phaseModel),
    mu0_("mu0", dimMass/(dimLength*dimTime), dict)
{

}

Foam::tmp<Foam::volScalarField> Foam::Plasma::mu() const
{
    return tmp<volScalarField>(
		new volScalarField(	
			IOobject
			(
				IOobject::groupName("mu", phase().name()),
				phase().mesh().time().timeName(),
				phase().mesh(),
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			phase().mesh(),
			mu0_
		)
	);
}

Foam::tmp<Foam::fvVectorMatrix> Foam::Plasma::divDevReff(Foam::volVectorField & U)
{
    return (
        - fvm::laplacian(phase().alpha() * mu0_/phase().rho(), U)
        - fvc::div(phase().alpha() * (mu0_/phase().rho() * dev2(T(fvc::grad(U)))/* + lambda0_*fvc::div(U)*/), "div(Rc)")
    );
}

Foam::tmp<Foam::fvVectorMatrix> Foam::Plasma::divDevRhoReff(Foam::volVectorField & U)
{
    return (
        - fvm::laplacian(phase().alpha() * mu0_, U)
        - fvc::div(phase().alpha() * (mu0_ * dev2(T(fvc::grad(U)))/* + lambda0_*fvc::div(U)*/), "div(Rc)")
    );
}

// RBC model
Foam::Rbc::Rbc(const Foam::dictionary & dict, Foam::PhaseModel & phaseModel)
:
    ViscosityModel(dict, phaseModel),
    k_("k", dimTime, dict),
    muAlpha_
    (
        IOobject
        (
            IOobject::groupName("muAlpha", phase().name()),
            phaseModel.mesh().time().timeName(),
            phaseModel.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseModel.mesh(),
        dimensionedScalar("muAlpha", dimMass/(dimLength*dimTime), 0)
    )
{
    this->correct();
}

Foam::tmp<Foam::volScalarField> Foam::Rbc::mu() const
{
    volScalarField limitedAlpha = max(scalar(0), min(scalar(1), phase().alpha()));
    dimensionedScalar dimVisc("dimVisc", dimMass/(dimLength*dimTime), 1);
    volScalarField mu0OverAlpha = dimVisc*(537.002*pow(limitedAlpha, 2) + 55.006*limitedAlpha - 0.129);
    volScalarField muInfOverAlpha = dimVisc*(27.873*pow(limitedAlpha, 2) - 21.218*limitedAlpha + 14.439);

    return tmp<volScalarField>(
		new volScalarField(	
			IOobject
			(
				IOobject::groupName("mu", phase().name()),
				phase().mesh().time().timeName(),
				phase().mesh(),
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			muInfOverAlpha + (mu0OverAlpha - muInfOverAlpha) * (1 + log(scalar(1) + k_*this->shearRate())) / (scalar(1) + k_*this->shearRate())
		)
	);
}

void Foam::Rbc::correct() 
{ 
    // Evaluate viscosity
    volScalarField limitedAlpha = max(scalar(0), min(scalar(1), phase().alpha()));
    dimensionedScalar dimVisc("dimVisc", dimMass/(dimLength*dimTime), 1);
    volScalarField mu0 = dimVisc*(537.002*pow(limitedAlpha, 3) + 55.006*pow(limitedAlpha, 2) - 0.129*limitedAlpha);
    volScalarField muInf = dimVisc*(27.873*pow(limitedAlpha, 3) - 21.218*pow(limitedAlpha, 2) + 14.439*limitedAlpha);

    muAlpha_ = muInf + (mu0 - muInf) * (1 + log(1 + k_*this->shearRate())) / (1 + k_*this->shearRate());
    //beta2_ = beta20_ * limitedAlpha * (1 + limitedAlpha);
}

Foam::tmp<Foam::fvVectorMatrix> Foam::Rbc::divDevRhoReff(Foam::volVectorField & U)
{
    return (
        - fvm::laplacian(muAlpha_, U)
        - fvc::div(muAlpha_ * dev2(T(fvc::grad(U))), "div(Rc)")
    );
}

Foam::tmp<Foam::fvVectorMatrix> Foam::Rbc::divDevReff(Foam::volVectorField & U)
{
    return (
        - fvm::laplacian(muAlpha_/phase().rho(), U)
        - fvc::div(muAlpha_/phase().rho() * dev2(T(fvc::grad(U))), "div(Rc)")
    );
}