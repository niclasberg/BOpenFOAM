#ifndef LEIGHTON_ACRIVOS_SOURCES_H_
#define LEIGHTON_ACRIVOS_SOURCES_H_

#include "volFields.H"
#include "fvMeshFunctionObject.H"
#include "writeLocalObjects.H"

namespace Foam 
{

namespace functionObjects 
{

class leightonAcrivosSources 
: 
    public fvMeshFunctionObject,
    public writeLocalObjects
{
public:
    static const constexpr char * alphaTermArrayName = "LAalphaTerm";
    static const constexpr char * gammaTermArrayName = "LAgammaTerm";
    static const constexpr char * muTermArrayName = "LAmuTerm";

    //- Runtime type information
    TypeName("leightonAcrivosSources");

    leightonAcrivosSources
    (
        const word & name,
        const Time& runTime,
        const dictionary& dict
    );


    //- Destructor
    virtual ~leightonAcrivosSources() { }

    virtual bool read(const dictionary&);
    virtual bool execute();
    virtual bool write();

private:
    word uName_;
    dimensionedScalar Db_;			// Brownian diffusion
	dimensionedScalar a_;			// Particle radius
	dimensionedScalar Kgamma_;
	dimensionedScalar Kalpha_;
	dimensionedScalar Kmu_;
    /*volVectorField vorticity_;
    volVectorField ddtVorticity_;
    volVectorField convectiveTerm_;
    volVectorField vorticityDotNablaU_;
    volVectorField vorticityDivU_;
    volVectorField baroclinicTerm_;
    volVectorField vorticitiyDiffusiveFlux_;*/
};

}

}

#endif