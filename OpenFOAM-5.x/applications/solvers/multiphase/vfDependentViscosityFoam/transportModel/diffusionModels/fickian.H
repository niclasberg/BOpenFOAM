#ifndef fickianDiffusion_H
#define fickianDiffusion_H

#include "isotropicDiffusion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fickian : public isotropicDiffusion
{

protected:
	dimensionedScalar Dab_;
public:
    TypeName("fickian");

    fickian(const dictionary&, const volVectorField&, const surfaceScalarField&, const viscosityModelC & viscosityModel);
    virtual ~fickian() {}

	virtual tmp<volScalarField> D(const volScalarField &) const;
};

} // End namespace Foam


#endif

// ************************************************************************* //
