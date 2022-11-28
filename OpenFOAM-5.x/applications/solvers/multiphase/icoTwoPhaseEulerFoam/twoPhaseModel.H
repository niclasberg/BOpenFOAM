#ifndef twoPhaseModel_H_
#define twoPhaseModel_H_

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "viscosityModels.H"

//class DragModel;
//class LiftModel;

namespace Foam
{

class PhaseModel 
{
public:
    PhaseModel(const fvMesh& mesh, const dictionary& dict, const word &);
    PhaseModel(const fvMesh& mesh, const dictionary& dict, const word &, const volScalarField &);

    const word & name() const { return name_; }
    const fvMesh & mesh() const { return mesh_; }

    volVectorField & U() { return U_; }
    const volVectorField & U() const { return U_; }

    volScalarField & alpha() { return alpha_; }
    const volScalarField & alpha() const { return alpha_; }

    surfaceScalarField & alphaPhi() { return alphaPhi_; }
    const surfaceScalarField & alphaPhi() const { return alphaPhi_; }

    surfaceScalarField & phi() { return phiPtr_(); }
    const surfaceScalarField & phi() const { return phiPtr_(); }

    dimensionedScalar & rho() { return rho0_; }
    const dimensionedScalar & rho() const { return rho0_; }

    ViscosityModel & viscosityModel() { return viscosityModel_(); }
    const ViscosityModel & viscosityModel() const { return viscosityModel_(); }

    void correctInflowOutflow(surfaceScalarField& alphaPhi) const;

private:
    void readPhi();

    word name_;
    dimensionedScalar rho0_;
    const fvMesh & mesh_;
    volScalarField alpha_;
    volVectorField U_;
    surfaceScalarField alphaPhi_;
    autoPtr<surfaceScalarField> phiPtr_;
    autoPtr<ViscosityModel> viscosityModel_;
};

class twoPhaseModel 
: 
    public IOdictionary
{
public:
    twoPhaseModel(const fvMesh& mesh);

    const fvMesh & mesh() const { return mesh_; }
    tmp<volScalarField> rho() const;

    surfaceScalarField & phi() { return phi_; }
    const surfaceScalarField & phi() const { return phi_; }

    tmp<volVectorField> U() const 
    { 
        return phase1_.alpha()*phase1_.U() + phase2_.alpha()*phase2_.U();;
    }

    tmp<volScalarField> dragForce() const;

    PhaseModel & phase1() { return phase1_; }
    PhaseModel & phase2() { return phase2_; }

    dimensionedScalar d() const { return d_; }

    //tmp<volScalarField> dragForceCoefficient() const;
    //tmp<volScalarField> liftForceCoefficient() const;

    void correct();
    void solve();
private:
    tmp<surfaceScalarField> calcPhi() const;

    const fvMesh & mesh_;
    PhaseModel phase1_;
    PhaseModel phase2_;
    dimensionedScalar d_;
    surfaceScalarField phi_;
};

} // Namespace Foam

#endif