#ifndef isotropicDiffusion_H
#define isotropicDiffusion_H

#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class isotropicDiffusion : public diffusionModel 
{
protected:
	dimensionedScalar D_;
	
public:
    TypeName("isotropicDiffusion");

    isotropicDiffusion(const volVectorField&, const surfaceScalarField&, const viscosityModelC & viscosityModel);
    virtual ~isotropicDiffusion() {}

    //- Correct the diffusion model
    virtual void correct(const Foam::volScalarField & alpha) { }

	virtual tmp<volScalarField> D(const volScalarField &) const = 0;

	//- Get the flux
	virtual tmp<volVectorField> flux(const volScalarField &) const;

	//- Divergence of flux
	virtual tmp<fvScalarMatrix> divFlux(const volScalarField &) const;

	scalarField noFluxBoundaryWeights(const fvPatch & patch) const override;
};

} // End namespace Foam


#endif

// ************************************************************************* //
