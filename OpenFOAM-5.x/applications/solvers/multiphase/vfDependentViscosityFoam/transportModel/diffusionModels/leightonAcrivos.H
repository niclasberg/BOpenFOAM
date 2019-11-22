#ifndef leightonAcrivos_H
#define leightonAcrivos_H

#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class leightonAcrivos : public diffusionModel 
{

protected:
	dimensionedScalar Db_;			// Brownian diffusion
	dimensionedScalar a_;			// Particle radius
	dimensionedScalar Kgamma_;
	dimensionedScalar Kalpha_;
	dimensionedScalar Kmu_;
	volScalarField Dalpha_;
	volVectorField Dgamma_;
	
public:
    TypeName("leightonAcrivos");

    leightonAcrivos(const dictionary&, const volVectorField&, const surfaceScalarField&, const viscosityModelC &);
    virtual ~leightonAcrivos() {}

    //- Correct the diffusion model
    void correct(const Foam::volScalarField & alpha) override;

	//- Get the flux
	virtual tmp<volVectorField> flux(const volScalarField &) const;

	//- Divergence of flux
	virtual tmp<fvScalarMatrix> divFlux(const volScalarField &) const;

	scalarField noFluxBoundaryWeights(const fvPatch & patch) const override;
};

} // End namespace Foam


#endif

// ************************************************************************* //

