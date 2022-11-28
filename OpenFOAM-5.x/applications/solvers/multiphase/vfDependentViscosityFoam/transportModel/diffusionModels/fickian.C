#include "fickian.H"
#include "volFields.H"
#include "fvMatrix.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(fickian, 0);
	addToRunTimeSelectionTable
    (
        diffusionModel,
        fickian,
        diffusion
    );
}	

Foam::fickian::fickian
(
    const dictionary& diffusionProperties,
    const volVectorField& U,
    const surfaceScalarField& phi, 
	const viscosityModelC & viscosityModel)
: 
	isotropicDiffusion(U, phi, viscosityModel),
	Dab_("Dab", dimLength*dimLength/dimTime, diffusionProperties)
{

}

Foam::tmp<Foam::volScalarField> Foam::fickian::D(const volScalarField & alpha) const
{
	return tmp<volScalarField>(
		new volScalarField(	
			IOobject
			(
				"Dab",
				U_.mesh().time().timeName(),
				U_.mesh(),
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			U_.mesh(),
			Dab_
		)
	);
}
