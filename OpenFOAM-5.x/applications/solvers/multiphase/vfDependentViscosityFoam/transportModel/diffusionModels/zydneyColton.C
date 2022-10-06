#include "zydneyColton.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(zydneyColton, 0);
	addToRunTimeSelectionTable
    (
        diffusionModel,
        zydneyColton,
        diffusion
    );
}	

Foam::zydneyColton::zydneyColton
(
    const dictionary& diffusionProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const viscosityModelC & viscosityModel)
: 
    isotropicDiffusion(U, phi, viscosityModel),
    D0_(diffusionProperties.lookup("D0")),
    a_(diffusionProperties.lookup("a")),
    k_(diffusionProperties.lookup("k")),
    n_(diffusionProperties.lookup("n"))
{

}

Foam::tmp<Foam::volScalarField> Foam::zydneyColton::D(const volScalarField & alpha) const
{
	return tmp<volScalarField>(
		new volScalarField(
			"Dab", 
			D0_ + k_ * pow(a_, 2) * alpha * this->shearRate() * pow(1. - alpha, n_)
		)
	);
}

