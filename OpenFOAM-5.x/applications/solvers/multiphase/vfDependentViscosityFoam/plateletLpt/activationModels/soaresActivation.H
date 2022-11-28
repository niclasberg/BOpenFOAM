#ifndef SOARESACTIVATION_H_
#define SOARESACTIVATION_H_
#include "activationModel.H"

namespace Foam {

template<class CloudType>
class soaresActivation : public activationModel<CloudType> {
public:
    TypeName("soares");

    soaresActivation(const dictionary & dict);
    virtual void calc(typename CloudType::parcelType &, scalar dt) const;

    virtual autoPtr<activationModel<CloudType> > clone() const
    {
        return autoPtr<activationModel<CloudType> >(new soaresActivation(*this));
    }

private:
    scalar Sr_;
    scalar C_;
    scalar alpha_;
    scalar beta_;
    scalar Cr_;
    scalar gamma_;
    scalar delta_;
    scalar pasMin_;
};

template<class CloudType>
soaresActivation<CloudType>::soaresActivation(const dictionary & dict)
: 
    activationModel<CloudType>(dict),
    Sr_(dict.lookupOrDefault<scalar>("Sr", 1.5701e-7)),
    C_(dict.lookupOrDefault<scalar>("C", 1.4854e-7)),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 1.4854)),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.4401)),
    Cr_(dict.lookupOrDefault<scalar>("Cr", 1.3889e-4)),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 0.5720)),
    delta_(dict.lookupOrDefault<scalar>("delta", 1.5125)),
    pasMin_(dict.lookupOrDefault<scalar>("pasMin", 1.0e-12))
{

}

template<class CloudType>
void soaresActivation<CloudType>::calc(typename CloudType::parcelType & p, scalar dt) const 
{
    const scalar PaToDyneInvCmSqr = 10;

    // Compute stress magnitude
    scalar tau_mag = PaToDyneInvCmSqr * this->stressNorm(p.tau());
    scalar tau_mag_last = PaToDyneInvCmSqr * this->stressNorm(p.tau_last());
    scalar dTauDt = abs(tau_mag - tau_mag_last) / dt;

    // Limit pas
    p.pas() = max(p.pas(), pasMin_);

    // Pas terms
    const scalar S = Sr_ * p.pas() * p.stressHistory();
    const scalar F = pow(C_, 1./beta_) * beta_ * pow(p.pas(), (beta_ - 1.)/beta_) * pow(tau_mag, alpha_/beta_);
    const scalar G = pow(Cr_, 1./delta_) * delta_ * pow(p.pas(), (delta_-1.)/delta_) * pow(p.stressRateHistory(), gamma_/delta_);

    // Update PAS
    p.pas() += dt*(S + F + G) * (1-p.pas());

    // Update remaining terms
    // Stress history is the time integral of the shear stress
    p.stressHistory() += 0.5 * dt * (tau_mag + tau_mag_last);

    // Stress rate history is the time average of the derivative of stress
    // Here we use a running average computation
    p.stressRateHistory() = (p.age() * p.stressRateHistory() + dTauDt * dt) / (p.age() + dt);
}

}

#endif