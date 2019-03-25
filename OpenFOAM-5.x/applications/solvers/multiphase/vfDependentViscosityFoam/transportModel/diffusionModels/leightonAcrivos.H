#ifndef leightonAcrivos_H
#define leightonAcrivos_H

#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class leightonAcrivos : public diffusionModel 
{

protected:
	dimensionedScalar a_;
	dimensionedScalar Kgamma_;
	dimensionedScalar Kalpha_;
	dimensionedScalar Kmu_;
	const volScalarField & mu_;
public:
    TypeName("leightonAcrivos");

    leightonAcrivos(const dictionary&, const volVectorField&, const surfaceScalarField&);
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

