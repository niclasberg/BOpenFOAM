#ifndef ACTIVATION_MODEL_C_
#define ACTIVATION_MODEL_C_

namespace Foam {

// ActivationModel
template<class CloudType>
autoPtr<activationModel<CloudType> > activationModel<CloudType>::New(
    const dictionary& dict ) 
{
    const word modelType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "activationModel::New(const dictionary&)"
        )   << "Unknown activaitonModel type "
            << modelType << nl << nl
            << "Valid activationModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<activationModel>(cstrIter()(dict));
}

template<class CloudType>
activationModel<CloudType>::activationModel(const dictionary & dict)
:
    stressNormType_(wordToStressNorm(dict.lookup("stressNorm")))
{

}

template<class CloudType>
typename activationModel<CloudType>::StressNorm activationModel<CloudType>::wordToStressNorm(const word & stressNormName)
{
    StressNorm sNorm = Frobenius; 
    if(stressNormName == "frobenius") {
        sNorm = Frobenius;
    } else if(stressNormName == "vonMises") {
        sNorm = vonMises;
    } else {
        FatalErrorIn
        (
            "activationModel::activationModel(const dictionary&)"
        )   << "Unknown stress norm type "
            << stressNormName << nl << nl
            << "Valid activationModels are : " << endl
            << "frobenius, vonMises"
            << exit(FatalError);
    }
    return sNorm;
}

template<class CloudType>
scalar activationModel<CloudType>::stressNorm(const symmTensor & tau) const
{
    switch(stressNormType_) {
    case Frobenius:
        return mag(tau) / sqrt(2.0);
    case vonMises:
        return sqrt(
                (1./6.) * (
                    (tau.xx() - tau.yy())*(tau.xx() - tau.yy()) +
                    (tau.yy() - tau.zz())*(tau.yy() - tau.zz()) + 
                    (tau.xx() - tau.zz())*(tau.xx() - tau.zz())
                ) + 
                tau.xy()*tau.xy() + tau.xz()*tau.xz() + tau.yz()*tau.yz()
        );
    }
    return 0;
}


}


#endif