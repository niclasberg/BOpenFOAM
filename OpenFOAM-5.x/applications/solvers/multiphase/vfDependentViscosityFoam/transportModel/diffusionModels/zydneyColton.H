#ifndef zydneyColton_H
#define zydneyColton_H

#include "isotropicDiffusion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class zydneyColton : public isotropicDiffusion 
{
	dimensionedScalar D0_;
	dimensionedScalar a_;
	dimensionedScalar k_;
	dimensionedScalar n_;

public:
    TypeName("zydneyColton");

    zydneyColton(const dictionary&, const volVectorField&, const surfaceScalarField&, const viscosityModelC & viscosityModel);
    virtual ~zydneyColton() {}

	virtual tmp<volScalarField> D(const volScalarField & alpha) const;
};

} // End namespace Foam


#endif

// ************************************************************************* //
