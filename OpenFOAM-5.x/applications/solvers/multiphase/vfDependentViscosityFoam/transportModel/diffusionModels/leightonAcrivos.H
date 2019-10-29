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
	dimensionedScalar Kc_;		// 
	dimensionedScalar Kmu_;
	
public:
    TypeName("leightonAcrivos");

    leightonAcrivos(const dictionary&, const volVectorField&, const surfaceScalarField&, const viscosityModelC &);
    virtual ~leightonAcrivos() {}

    //- Correct the diffusion model
    virtual void correct() { }

	//- Get the flux
	virtual tmp<volVectorField> flux(const volScalarField &) const;

	//- Divergence of flux
	virtual tmp<fvScalarMatrix> divFlux(const volScalarField &) const;
};

} // End namespace Foam


#endif

// ************************************************************************* //

