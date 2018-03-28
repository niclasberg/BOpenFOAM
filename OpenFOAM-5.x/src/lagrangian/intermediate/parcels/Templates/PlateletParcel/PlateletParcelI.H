// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::PlateletParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties(),
    a_(this->dict_, 0.0),
    b_(this->dict_, 0.0),
	c_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::PlateletParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    a_(cp.a_),
    b_(cp.b_),
	c_(cp.c_)
{}


template<class ParcelType>
inline Foam::PlateletParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    a_(this->dict_, "a"),
    b_(this->dict_, "b"),
	c_(this->dict_, "c")
{}


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
    mechanicalDose_(0.0),
	pas_(0.0),
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
    mechanicalDose_(0.0),
	pas_(0.0),
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
    mechanicalDose_(0.0),
    pas_(0.0),
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


template<class ParcelType>
template<class CloudType>
inline const Foam::interpolation<Foam::symmTensor>&
Foam::PlateletParcel<ParcelType>::TrackingData<CloudType>::tauInterp() const
{
    return tauInterp_();
}

// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar 
Foam::PlateletParcel<ParcelType>::constantProperties::a() const
{
	return a_.value();;
}

template<class ParcelType>
inline Foam::scalar 
Foam::PlateletParcel<ParcelType>::constantProperties::b() const
{
	return b_.value();;
}

template<class ParcelType>
inline Foam::scalar 
Foam::PlateletParcel<ParcelType>::constantProperties::c() const
{
	return c_.value();;
}
