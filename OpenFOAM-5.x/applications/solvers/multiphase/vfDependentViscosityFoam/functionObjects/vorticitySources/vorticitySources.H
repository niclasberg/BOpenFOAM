#ifndef VORTICITY_SOURCES_H_
#define VORTICITY_SOURCES_H_

#include "volFields.H"
#include "fvMeshFunctionObject.H"
#include "writeLocalObjects.H"

namespace Foam 
{

namespace functionObjects 
{

class vorticitySources 
: 
    public fvMeshFunctionObject,
    public writeLocalObjects
{
public:
    static const constexpr char * vorticityArrayName = "vorticity";
    static const constexpr char * ddtVorticityArrayName = "ddt(vorticity)";
    static const constexpr char * convectiveTermName = "dot(U,nabla)vorticity";
    static const constexpr char * vortexStretchingTermName = "dot(vorticity,nabla)U";
    static const constexpr char * strainRateName = "strainRate";
    static const constexpr char * divTermName = "div(U)vorticity";
    static const constexpr char * baroclinicTermName = "cross(nabla_rho,nabla_p)_rhoSqr";
    static const constexpr char * diffusiveFluxTermName = "curl(div(tau)_rho)";

    //- Runtime type information
    TypeName("vorticitySources");

    vorticitySources
    (
        const word & name,
        const Time& runTime,
        const dictionary& dict
    );


    //- Destructor
    virtual ~vorticitySources() { }

    virtual bool read(const dictionary&);
    virtual bool execute();
    virtual bool write();

private:
    word uName_;
    word rhoName_;
    word pName_;
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