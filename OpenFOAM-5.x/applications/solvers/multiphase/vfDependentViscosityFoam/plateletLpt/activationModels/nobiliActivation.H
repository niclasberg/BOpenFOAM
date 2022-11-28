#ifndef NOBILIACTIVATION_H_
#define NOBILIACTIVATION_H_
#include "activationModel.H"

namespace Foam {

template<class CloudType>
class nobiliActivation : public activationModel<CloudType> {
public:
    TypeName("nobili");
    
    nobiliActivation(const dictionary & dict);
    virtual void calc(typename CloudType::parcelType &, scalar dt) const;


    virtual autoPtr<activationModel<CloudType> > clone() const
    {
        return autoPtr<activationModel<CloudType> >(new nobiliActivation(*this));
    }

private:
    scalar a_;
    scalar b_;
    scalar c_;
};

template<class CloudType>
nobiliActivation<CloudType>::nobiliActivation(const dictionary & dict)
:
    activationModel<CloudType>(dict),
    a_(dict.lookupOrDefault<scalar>("a", 1.3198)), 
    b_(dict.lookupOrDefault<scalar>("b", 0.6256)), 
    c_(dict.lookupOrDefault<scalar>("c", 1.0e-5))
{

}

template<class CloudType>
void nobiliActivation<CloudType>::calc(typename CloudType::parcelType & p, scalar dt)  const 
{	
    const scalar PaToDyneInvCmSqr = 10;

    // Compute stress magnitude
    scalar tau_mag = PaToDyneInvCmSqr * this->stressNorm(p.tau());
    scalar tau_mag_last = PaToDyneInvCmSqr * this->stressNorm(p.tau_last());

	// Update mechanical dose
	scalar stressHistoryLast = p.stressHistory();
	p.stressHistory() += (dt / 2.0) * (Foam::pow(tau_mag, b_/a_) + Foam::pow(tau_mag_last, b_ / a_));

	// Update PAS
	p.pas() += (c_ * a_ * dt / 2.0) * 
		(Foam::pow(stressHistoryLast, a_-1) * Foam::pow(tau_mag_last, b_/a_) + 
		 Foam::pow(p.stressHistory(), a_-1) * Foam::pow(tau_mag, b_/a_));
}

}

#endif