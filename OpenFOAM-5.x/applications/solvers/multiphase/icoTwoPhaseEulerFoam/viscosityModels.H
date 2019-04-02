#ifndef viscosityModel_H_
#define viscosityModel_H_

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrix.H"
#include "runTimeSelectionTables.H"

namespace Foam {

class PhaseModel;

class ViscosityModel {
public:
    ViscosityModel(const dictionary & dict, PhaseModel & phaseModel);
    virtual ~ViscosityModel() { }

    //- Runtime type information
    TypeName("ViscosityModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ViscosityModel,
        viscModel,
        (
            const dictionary& dict,
            PhaseModel & phaseModel
        ),
        (dict, phaseModel)
    );

    static autoPtr<ViscosityModel> New
    (
        const dictionary& dict,
        PhaseModel & phaseModel
    );

    virtual void correct() = 0;
    virtual tmp<volScalarField> mu() const = 0;
    virtual tmp<fvVectorMatrix> divDevReff(volVectorField & U) = 0;
    virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField & U) = 0;

    PhaseModel & phase() { return phaseModel_; }
    const PhaseModel & phase() const  { return phaseModel_; }
    
    tmp<volScalarField> shearRate() const;

private:
    PhaseModel & phaseModel_;
};

class Plasma : public ViscosityModel {
public:
    Plasma(const dictionary & dict, PhaseModel & phaseModel);

    TypeName("plasma");

    virtual void correct() { }
    
    virtual tmp<volScalarField> mu() const;
    virtual tmp<fvVectorMatrix> divDevReff(volVectorField & U);
    virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField & U);
private:
    dimensionedScalar mu0_;
    //dimensionedScalar lambda0_;
};

class Rbc  : public ViscosityModel {
public:
    Rbc(const dictionary & dict, PhaseModel & phaseModel);  

    TypeName("rbc");

    virtual void correct();
    virtual tmp<volScalarField> mu() const;
    virtual tmp<fvVectorMatrix> divDevReff(volVectorField & U);
    virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField & U);
private:
    //dimensionedScalar beta20_;
    dimensionedScalar k_;
    volScalarField muAlpha_;
    //volScalarField beta2_;
};

} //Namespace foam

#endif