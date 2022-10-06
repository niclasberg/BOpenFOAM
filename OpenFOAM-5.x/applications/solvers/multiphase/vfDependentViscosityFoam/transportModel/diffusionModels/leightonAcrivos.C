#include "leightonAcrivos.H"
#include "volFields.H"
#include "fvMatrix.H"
#include "fvm.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(leightonAcrivos, 0);
	addToRunTimeSelectionTable
    (
        diffusionModel,
        leightonAcrivos,
        diffusion
    );
}	

Foam::leightonAcrivos::leightonAcrivos
(
    const dictionary& diffusionProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const viscosityModelC & viscosityModel)
:
    diffusionModel(U, phi, viscosityModel),
    Db_("Db", dimLength*dimLength/dimTime, diffusionProperties),
    a_("a", dimLength, diffusionProperties),
    Kgamma_("Kgamma", dimless, diffusionProperties),
    Kalpha_("Kalpha", dimless, diffusionProperties),
    Kmu_("Kmu", dimless, diffusionProperties),
    Dalpha_(
        IOobject
        (
            "Dalpha",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("Dalpha", dimLength*dimLength/dimTime, 0)
    ),
    Dgamma_(
        IOobject
        (
            "Dgamma",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedVector("Dgamma", dimLength/dimTime, Zero)
    )
{
    
}

void Foam::leightonAcrivos::correct(const Foam::volScalarField & alpha)
{
    const volScalarField gamma = shearRate();
    Dalpha_ = Db_ + a_*a_*alpha*gamma *(Kalpha_ + Kmu_* alpha * viscosityModel().dMuDalpha() / viscosityModel().mu());
    Dgamma_ = a_*a_*alpha*(Kgamma_ + Kmu_ * gamma * viscosityModel().dMuDgamma() / viscosityModel().mu()) * fvc::grad(gamma);
}

Foam::tmp<Foam::volVectorField> Foam::leightonAcrivos::flux(const Foam::volScalarField & alpha) const
{
	return tmp<volVectorField>(
        new volVectorField(
            IOobject
            (
                "flux",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Dalpha_ * fvc::grad(alpha) + Dgamma_ * alpha
        )
    );
}

//- Divergence of flux
Foam::tmp<Foam::fvScalarMatrix> Foam::leightonAcrivos::divFlux(const Foam::volScalarField & alpha) const
{
    // Surface field for the alpha diffusivity
    surfaceScalarField Dalphaf = fvc::interpolate(Dalpha_);
    surfaceScalarField Dgammaf = fvc::interpolate(Dgamma_) & U_.mesh().Sf();

	return fvm::laplacian(Dalphaf, alpha) + fvm::div(Dgammaf, alpha);
}

Foam::scalarField Foam::leightonAcrivos::noFluxBoundaryWeights(const Foam::fvPatch & patch) const 
{
    const scalarField projectedDgamma = Dgamma_.boundaryField()[patch.index()] & patch.nf();
    return projectedDgamma / (projectedDgamma + Dalpha_.boundaryField()[patch.index()] * patch.deltaCoeffs());
}