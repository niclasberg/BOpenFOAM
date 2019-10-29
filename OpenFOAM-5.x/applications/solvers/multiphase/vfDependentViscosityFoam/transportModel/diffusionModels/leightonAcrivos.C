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
    Kc_("Kc", dimless, diffusionProperties),
    Kmu_("Kmu", dimless, diffusionProperties)
{
    
}

Foam::tmp<Foam::volVectorField> Foam::leightonAcrivos::flux(const Foam::volScalarField & alpha) const
{
    const volScalarField gamma = shearRate();

	return -(Db_ + a_*a_ * alpha * gamma * 
                (Kc_ + Kmu_*alpha*viscosityModel().dMuDalpha() / this->viscosityModel().mu())) * fvc::grad(alpha)
           - a_*a_ * gamma * alpha * Kc_ * fvc::grad(gamma);
}

//- Divergence of flux
Foam::tmp<Foam::fvScalarMatrix> Foam::leightonAcrivos::divFlux(const Foam::volScalarField & alpha) const
{
    const volScalarField gamma = shearRate();

    // Surface field for the alpha diffusivity
    surfaceScalarField dAlpha = fvc::interpolate(Db_ + a_*a_*alpha*gamma *(Kc_ + Kmu_* viscosityModel().dMuDalpha() / viscosityModel().mu())) * U_.mesh().magSf();
    surfaceScalarField dShear = fvc::interpolate(a_*a_*alpha*Kc_*fvc::grad(gamma)) & U_.mesh().Sf();

	return fvm::laplacian(dAlpha, alpha) + fvm::div(dShear, alpha);
}
