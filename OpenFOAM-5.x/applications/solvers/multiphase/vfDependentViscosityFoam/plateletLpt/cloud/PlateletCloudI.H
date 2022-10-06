template<class CloudType>
inline const Foam::PlateletCloud<CloudType>&
Foam::PlateletCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}

template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::PlateletCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::PlateletCloud<CloudType>::constProps()
{
    return constProps_;
}
