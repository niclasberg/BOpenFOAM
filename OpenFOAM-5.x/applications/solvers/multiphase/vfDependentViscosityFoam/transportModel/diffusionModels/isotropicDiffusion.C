#include "isotropicDiffusion.H"
#include "volFields.H"
#include "fvMatrix.H"
#include "fvm.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(isotropicDiffusion, 0);
}	

Foam::isotropicDiffusion::isotropicDiffusion
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:  diffusionModel(U, phi)
{

}

Foam::tmp<Foam::volVectorField> Foam::isotropicDiffusion::flux(const Foam::volScalarField & alpha) const
{
	return this->D(alpha) * fvc::grad(alpha);
}

//- Divergence of flux
Foam::tmp<Foam::fvScalarMatrix> Foam::isotropicDiffusion::divFlux(const Foam::volScalarField & alpha) const
{
	return fvm::laplacian(this->D(alpha), alpha);
}
