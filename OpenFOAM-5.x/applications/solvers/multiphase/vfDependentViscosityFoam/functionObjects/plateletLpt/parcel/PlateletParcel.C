#include "PlateletParcel.H"
#include "meshTools.H"

// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::PlateletParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::setCellValues(td, dt, cellI);

	// Store last value
	tau_last_ = tau_;

	// Interpolate shear stress tetIndices tetIs = this->currentTetIndices();
    tau_ = td.tauInterp().interpolate
    (
        this->coordinates(),
        this->currentTetIndices()
    );
}

template<class ParcelType>
template<class TrackData>
void Foam::PlateletParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
	// Run parent calc method
	ParcelType::calc(td, dt, cellI);

	//typedef typename TrackData::cloudType cloudType;
	//typedef typename cloudType::activationType activationType;

	//activationType & activation = td.cloud().activation();

	// Compute stress magnitude 
	// OpenFoam computes the tensor magnitude as mag(T) = (T:T)^(1/2)
	//   which gives a factor sqrt(2) error for unidirectional flow fields 
	//   (e.g. the Couette flow)
	scalar tau_mag = mag(tau_) / sqrt(2.0);
	scalar tau_mag_last = mag(tau_last_) / sqrt(2.0);
	scalar b_over_a = td.cloud().constProps().b() / td.cloud().constProps().a();
	scalar a_minus_one = td.cloud().constProps().a() - 1;
	
	// Update mechanical dose
	scalar doseLast = mechanicalDose_;
	mechanicalDose_ += (dt / 2.0) * (Foam::pow(tau_mag, b_over_a) + Foam::pow(tau_mag_last, b_over_a));

	// Update PAS
	pas_ += (td.cloud().constProps().c()*td.cloud().constProps().a() * dt / 2.0) * 
		(Foam::pow(doseLast, a_minus_one)        * Foam::pow(tau_mag_last, b_over_a) + 
		 Foam::pow(mechanicalDose_, a_minus_one) * Foam::pow(tau_mag, b_over_a));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const PlateletParcel<ParcelType>& p
)
:
    ParcelType(p),
    mechanicalDose_(p.mechanicalDose_),
    pas_(p.pas_),
    tau_(p.tau_),
	tau_last_(p.tau_last_)
{}

template<class ParcelType>
Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const PlateletParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    mechanicalDose_(p.mechanicalDose_),
    pas_(p.pas_),
    tau_(p.tau_),
	tau_last_(p.tau_last_)
{}

