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
    const surfaceScalarField& phi)
: 
diffusionModel(U, phi),
a_(diffusionProperties.lookup("a")),
Kgamma_(diffusionProperties.lookup("Kgamma")),
Kalpha_(diffusionProperties.lookup("Kalpha")),
Kmu_(diffusionProperties.lookup("Kmu")), 
mu_(U.mesh().lookup<volScalarField>("mu"))
{

}

Foam::tmp<Foam::volVectorField> Foam::leightonAcrivos::flux(const Foam::volScalarField & alpha) const
{
	return -pow(a, 2) * alpha * (
				shearRate()*(Kalpha_*fvc::grad(alpha) + Kmu_*alpha*fvc::grad(mu_) / mu_) + 
				Kgamma_*alpha*fvc::grad(shearRate()));
}

//- Divergence of flux
Foam::tmp<Foam::fvScalarMatrix> Foam::leightonAcrivos::divFlux(const Foam::volScalarField & alpha) const
{
	return fvm::laplacian(this->D(alpha), alpha);
}
