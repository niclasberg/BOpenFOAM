template<class ParcelType>
inline Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const polyMesh& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(owner, position, cellI, tetFaceI, tetPtI),
    stressHistory_(0.0),
	pas_(0.0),
    stressRateHistory_(0.0),
    tau_(symmTensor::zero),
    tau_last_(symmTensor::zero)
{}

template<class ParcelType>
inline Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const polyMesh& owner,
    const vector& position,
    const label cellI
)
:
    ParcelType(owner, position, cellI),
    stressHistory_(0.0),
	pas_(0.0),
    stressRateHistory_(0.0),
    tau_(symmTensor::zero),
    tau_last_(symmTensor::zero)
{}


template<class ParcelType>
inline Foam::PlateletParcel<ParcelType>::PlateletParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    stressHistory_(0.0),
    pas_(0.0),
    stressRateHistory_(0.0),
	tau_(symmTensor::zero),
    tau_last_(symmTensor::zero)
{}

	
template<class ParcelType>
template<class CloudType>
inline Foam::PlateletParcel<ParcelType>::TrackingData<CloudType>::TrackingData
(
    CloudType& cloud,
    trackPart part
)
:
    ParcelType::template TrackingData<CloudType>(cloud, part),
    tauInterp_
    (
        interpolation<symmTensor>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.tau()
        )
    )
{}


