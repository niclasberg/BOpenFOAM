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

	// Update activation state
	td.cloud().activation().calc(*this, dt);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const PlateletParcel<ParcelType>& p
)
:
    ParcelType(p),
    stressHistory_(p.stressHistory_),
    pas_(p.pas_),
    stressRateHistory_(p.stressRateHistory_),
    tau_(p.tau_),
	tau_last_(p.tau_last_)
{
	
}

template<class ParcelType>
Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const PlateletParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    stressHistory_(p.stressHistory_),
    pas_(p.pas_),
    stressRateHistory_(p.stressRateHistory_),
    tau_(p.tau_),
	tau_last_(p.tau_last_)
{}

