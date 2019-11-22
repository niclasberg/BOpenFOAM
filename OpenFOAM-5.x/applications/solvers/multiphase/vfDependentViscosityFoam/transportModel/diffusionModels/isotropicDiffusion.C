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
    const surfaceScalarField& phi,
    const viscosityModelC & viscosityModel
)
: diffusionModel(U, phi, viscosityModel)
{

}

Foam::tmp<Foam::volVectorField> Foam::isotropicDiffusion::flux(const Foam::volScalarField & alpha) const
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
            this->D(alpha) * fvc::grad(alpha)
        )
    );
}

//- Divergence of flux
Foam::tmp<Foam::fvScalarMatrix> Foam::isotropicDiffusion::divFlux(const Foam::volScalarField & alpha) const
{
	return fvm::laplacian(this->D(alpha), alpha);
}

Foam::scalarField Foam::isotropicDiffusion::noFluxBoundaryWeights(const Foam::fvPatch & patch) const 
{
    // Return weight corresponding to a zero gradient condition
    return Foam::scalarField(patch.size(), 1.);
}